"""Real LIANA+ ligand--receptor communication inference.

The capability in this module is deliberately independent from the legacy
``cell_cell_communication`` implementation.  It calls LIANA+'s
``mt.rank_aggregate`` directly for every donor/condition stratum and keeps
the resulting rank table intact.  Donors are never pooled into one AnnData
object, and condition comparisons operate on donor-level observations only.

LIANA+ and its resource files are optional.  A missing package, missing local
resource, or invalid metadata is an explicit failure; there is no local
fallback implementation that could be mistaken for a LIANA result.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import re
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


LIANA_AUDIT_REQUIREMENTS = (
    "liana_resource_provenance_check",
    "donor_condition_partition_check",
    "liana_rank_output_check",
    "all_ligand_receptor_pairs_retained",
    "donor_level_condition_comparison_check",
    "spatial_data_use_check",
)


class LianaCommunicationDependencyError(ImportError):
    """Raised when the real LIANA+ package cannot be imported."""


class LianaCommunicationInputError(ValueError):
    """Raised when the communication contract is incomplete or invalid."""


_SPECIES_ALIASES = {
    "human": {"human", "homo_sapiens", "homo sapiens", "9606"},
    "mouse": {"mouse", "mus_musculus", "mus musculus", "10090"},
    "rat": {"rat", "rattus_norvegicus", "rattus norvegicus", "10116"},
    "zebrafish": {"zebrafish", "danio_rerio", "danio rerio", "7955"},
}
_RESOURCE_SPECIES_COLUMNS = (
    "species",
    "organism",
    "species_a",
    "species_b",
    "ligand_species",
    "receptor_species",
    "source_species",
    "target_species",
)


def _normalise_species(value: Any) -> str:
    return re.sub(r"[\s\-]+", "_", str(value).strip().lower())


def _species_compatible(expected: Any, observed: Any) -> bool:
    expected_norm = _normalise_species(expected)
    observed_norm = _normalise_species(observed)
    if expected_norm == observed_norm:
        return True
    for aliases in _SPECIES_ALIASES.values():
        normalised = {_normalise_species(item) for item in aliases}
        if expected_norm in normalised and observed_norm in normalised:
            return True
    return False


def _pair_set(frame: pd.DataFrame) -> set[Tuple[str, str]]:
    return {
        (str(ligand).strip(), str(receptor).strip())
        for ligand, receptor in frame[["ligand", "receptor"]].itertuples(index=False, name=None)
    }


def _complex_components(value: Any) -> set[str]:
    """Return conservative gene tokens for a simple LIANA complex label."""

    label = str(value).strip()
    if not label:
        return set()
    return {token for token in re.split(r"[_+:|]", label) if token}


def _pair_is_gene_compatible(pair: Tuple[str, str], gene_names: set[str]) -> bool:
    return all(
        component in gene_names
        for member in pair
        for component in _complex_components(member)
    )


def _missing_pair_records(
    pairs: Iterable[Tuple[str, str]],
    gene_names: set[str],
) -> list[Dict[str, Any]]:
    """Describe resource pairs that cannot be evaluated in this assay."""

    records: list[Dict[str, Any]] = []
    for ligand, receptor in sorted(pairs):
        missing_genes = sorted(
            {
                component
                for member in (ligand, receptor)
                for component in _complex_components(member)
                if component not in gene_names
            }
        )
        records.append(
            {
                "ligand": str(ligand),
                "receptor": str(receptor),
                "missing_genes": missing_genes,
            }
        )
    return records


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _as_scdata(payload: Any) -> SCData:
    if isinstance(payload, SCData):
        return payload.copy()
    if isinstance(payload, Mapping) and "X" in payload and "obs" in payload:
        return SCData.from_dict(dict(payload))
    if hasattr(payload, "X") and hasattr(payload, "obs") and hasattr(payload, "var"):
        return SCData.from_anndata(payload)
    raise LianaCommunicationInputError(
        "LIANA communication requires an AnnData-like or SCData artifact"
    )


def _as_dense(value: Any) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    elif hasattr(value, "A"):
        value = value.A
    return np.asarray(value)


def _resolve_column(
    frame: pd.DataFrame,
    explicit: Any,
    candidates: Sequence[str],
    label: str,
) -> str:
    if explicit is not None:
        name = str(explicit)
        if name not in frame.columns:
            raise LianaCommunicationInputError(
                f"{label} metadata column {name!r} is absent; available columns: "
                f"{list(frame.columns)}"
            )
        return name
    for name in candidates:
        if name in frame.columns:
            return name
    raise LianaCommunicationInputError(
        f"{label} metadata is required; provide {label.lower()}_col explicitly"
    )


def _resource_path(parameters: Mapping[str, Any]) -> Path:
    raw = parameters.get("lr_resource_path", parameters.get("resource_path"))
    if raw is None:
        raw = parameters.get("ligand_receptor_resource_path")
    if raw is None:
        raise LianaCommunicationInputError(
            "an explicit local ligand/receptor resource path is required "
            "(lr_resource_path)"
        )
    value = str(raw).strip()
    if not value or "://" in value:
        raise LianaCommunicationInputError(
            "LIANA ligand/receptor resources must be a local file path; "
            "remote URLs and named-resource downloads are not allowed"
        )
    path = Path(value).expanduser()
    if not path.is_file():
        raise FileNotFoundError(f"ligand/receptor resource does not exist: {path}")
    return path.resolve()


def _read_resource(path: Path) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Read and validate a local LIANA resource without downloading anything."""

    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t")
    elif suffix == ".csv":
        frame = pd.read_csv(path)
    elif suffix in {".json", ".jsonl"}:
        frame = pd.read_json(path, lines=suffix == ".jsonl")
    elif suffix in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
    else:
        raise LianaCommunicationInputError(
            f"unsupported local ligand/receptor resource format {suffix!r}; "
            "use CSV, TSV, JSON, JSONL, or Parquet"
        )
    if not isinstance(frame, pd.DataFrame):
        raise LianaCommunicationInputError("ligand/receptor resource must decode to a table")

    # Accept case variants, while preserving every additional LIANA resource
    # column (complexes, evidence, and source annotations) for the actual API.
    columns = {str(col).strip().lower(): col for col in frame.columns}
    ligand_col = columns.get("ligand")
    receptor_col = columns.get("receptor")
    if ligand_col is None or receptor_col is None:
        raise LianaCommunicationInputError(
            "local ligand/receptor resource must contain 'ligand' and 'receptor' columns"
        )
    frame = frame.copy()
    frame = frame.rename(columns={ligand_col: "ligand", receptor_col: "receptor"})
    frame["ligand"] = frame["ligand"].astype("string").str.strip()
    frame["receptor"] = frame["receptor"].astype("string").str.strip()
    frame = frame.loc[
        frame["ligand"].notna()
        & frame["receptor"].notna()
        & frame["ligand"].ne("")
        & frame["receptor"].ne("")
    ].copy()
    frame = frame.drop_duplicates(subset=["ligand", "receptor"], keep="first").reset_index(drop=True)
    if frame.empty:
        raise LianaCommunicationInputError("local ligand/receptor resource contains no valid pairs")
    # A two-column custom resource is valid, but when an external resource
    # declares species we must validate it rather than trusting a provenance
    # string supplied in the task.  Keep these values in metadata so the
    # independent auditor can repeat the check from the immutable file.
    species_columns = [
        str(column)
        for column in frame.columns
        if str(column).strip().lower() in _RESOURCE_SPECIES_COLUMNS
    ]
    species_values = sorted({
        str(value).strip()
        for column in species_columns
        for value in frame[column].dropna().tolist()
        if str(value).strip()
    })
    return frame, {
        "sha256": digest,
        "path": str(path),
        "n_pairs": int(len(frame)),
        "species_columns": species_columns,
        "species_values": species_values,
    }


def _load_liana() -> Any:
    try:
        import liana as li
    except ImportError as exc:  # pragma: no cover - depends on optional environment
        raise LianaCommunicationDependencyError(
            "LIANA+ is required for liana_communication; install the pinned "
            "local LIANA dependency before running this capability"
        ) from exc
    if not hasattr(li, "mt") or not hasattr(li.mt, "rank_aggregate"):
        raise LianaCommunicationDependencyError(
            "installed LIANA package does not expose liana.mt.rank_aggregate"
        )
    return li


def _frame_from_liana(value: Any) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        frame = value.copy()
    elif isinstance(value, Mapping):
        # LIANA currently returns a DataFrame for inplace=False.  Handle the
        # documented dict-shaped return without selecting a fake algorithm.
        if "liana_res" in value:
            frame = value["liana_res"]
        elif len(value) == 1:
            frame = next(iter(value.values()))
        else:
            raise LianaCommunicationInputError(
                "LIANA rank_aggregate returned an ambiguous mapping"
            )
        if not isinstance(frame, pd.DataFrame):
            frame = pd.DataFrame(frame)
    else:
        raise LianaCommunicationInputError(
            "LIANA rank_aggregate did not return a DataFrame"
        )
    if frame.empty:
        raise LianaCommunicationInputError(
            "LIANA rank_aggregate returned no ligand/receptor rows for a donor-condition group"
        )
    frame = frame.reset_index(drop=True)
    # LIANA's resource output can expose complexes only.  Keep the original
    # columns and provide stable pair columns for downstream grouping.
    if "ligand" not in frame and "ligand_complex" in frame:
        frame["ligand"] = frame["ligand_complex"]
    if "receptor" not in frame and "receptor_complex" in frame:
        frame["receptor"] = frame["receptor_complex"]
    required = {"source", "target", "ligand", "receptor", "magnitude_rank", "specificity_rank"}
    missing = required - set(frame.columns)
    if missing:
        raise LianaCommunicationInputError(
            "LIANA rank_aggregate output is missing required columns: "
            + ", ".join(sorted(missing))
        )
    for col in ("magnitude_rank", "specificity_rank"):
        values = pd.to_numeric(frame[col], errors="coerce")
        numeric = values.to_numpy(dtype=float)
        finite = np.isfinite(numeric)
        if ((finite) & ((numeric < 0) | (numeric > 1))).any():
            raise LianaCommunicationInputError(
                f"LIANA rank output column {col!r} contains a finite value outside [0, 1]"
            )
        frame[col] = values.astype(float)
    for col in ("source", "target", "ligand", "receptor"):
        values = frame[col]
        if values.isna().any() or values.astype(str).str.strip().eq("").any():
            raise LianaCommunicationInputError(
                f"LIANA rank output column {col!r} contains empty labels"
            )
    return frame


def _prepare_adata(data: SCData, gene_names: Sequence[str]) -> Any:
    """Create the independent AnnData passed to one LIANA invocation."""

    adata = data.to_anndata()
    # LIANA resolves genes via var_names.  SCData stores the stable symbol in
    # a column so make the conversion explicit and auditable.
    adata.var_names = pd.Index([str(value) for value in gene_names])
    if not adata.var_names.is_unique:
        raise LianaCommunicationInputError("gene names must be unique for LIANA inference")
    if not adata.obs_names.is_unique:
        adata.obs_names = pd.Index([f"cell_{i}" for i in range(adata.n_obs)])
    return adata


def _output_uri(contract: TaskContract, index: int, study_id: str) -> Optional[str]:
    if len(contract.expected_outputs) > 2:
        raise LianaCommunicationInputError(
            "liana_communication accepts at most two expected output artifacts"
        )
    if len(contract.expected_outputs) > index:
        uri = ArtifactURI.parse(contract.expected_outputs[index]).to_string()
    elif index == 0:
        uri = f"table://{study_id}/liana_communication/v1"
    else:
        return None
    parsed = ArtifactURI.parse(uri)
    if parsed.study_id != study_id or parsed.scheme != "table":
        raise LianaCommunicationInputError(
            "LIANA communication outputs must be table:// URIs in the input study"
        )
    return uri


def _validate_resource_species(
    species: str,
    resource_meta: Mapping[str, Any],
    parameters: Mapping[str, Any],
) -> Optional[str]:
    declared = parameters.get("lr_resource_species", parameters.get("resource_species"))
    if declared is not None and not str(declared).strip():
        raise LianaCommunicationInputError("lr_resource_species must not be empty")
    if declared is not None and not _species_compatible(species, declared):
        raise LianaCommunicationInputError(
            f"ligand/receptor resource species {declared!r} does not match requested species {species!r}"
        )
    observed = [str(value) for value in resource_meta.get("species_values", [])]
    incompatible = [value for value in observed if not _species_compatible(species, value)]
    if incompatible:
        raise LianaCommunicationInputError(
            f"ligand/receptor resource species values {incompatible!r} do not match requested species {species!r}"
        )
    if declared is not None:
        return str(declared).strip()
    if observed:
        # The resource may contain one alias per row.  Preserve the complete
        # set in provenance elsewhere; this value is only a compact summary.
        return ",".join(sorted(set(observed)))
    return None


def _validate_external_resource_hash(
    parameters: Mapping[str, Any],
    observed_sha256: str,
) -> None:
    pinned = parameters.get("external_resource_sha256")
    if pinned is not None and not isinstance(pinned, Mapping):
        raise LianaCommunicationInputError(
            "external_resource_sha256 must map the local resource parameter to a SHA-256 digest"
        )
    expected: list[str] = []
    if isinstance(pinned, Mapping):
        for key in (
            "lr_resource_path",
            "resource_path",
            "ligand_receptor_resource_path",
        ):
            if key in pinned:
                expected.append(str(pinned[key]))
        if not expected and len(pinned) == 1:
            # A caller may use the canonical resource alias as the only key;
            # do not guess which unrelated external file was intended.
            only_key, only_value = next(iter(pinned.items()))
            if str(only_key).lower() in {"lr_resource", "ligand_receptor_resource"}:
                expected.append(str(only_value))
            else:
                raise LianaCommunicationInputError(
                    "external_resource_sha256 does not contain the ligand/receptor resource key"
                )
    for value in expected:
        digest = value.strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or digest != observed_sha256.lower():
            raise LianaCommunicationInputError(
                "ligand/receptor resource SHA-256 does not match the pinned external resource hash"
            )
    direct = parameters.get("lr_resource_sha256", parameters.get("resource_sha256"))
    if direct is not None:
        digest = str(direct).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or digest != observed_sha256.lower():
            raise LianaCommunicationInputError(
                "ligand/receptor resource SHA-256 does not match the declared resource_sha256"
            )


def _safe_rank_comparison_pvalue(
    values_a: Sequence[float],
    values_b: Sequence[float],
    *,
    paired: bool,
) -> float:
    """Return a p-value only for an estimable donor-level rank test.

    SciPy can return zero or a misleading finite value for a zero-variance
    sample while emitting only a warning.  A communication comparison must
    retain that interaction as ``not_estimated`` instead of presenting a
    degenerate test as evidence.
    """

    a = np.asarray(values_a, dtype=float)
    b = np.asarray(values_b, dtype=float)
    if a.ndim != 1 or b.ndim != 1 or len(a) < 2 or len(b) < 2:
        return float("nan")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return float("nan")
    if paired:
        if len(a) != len(b):
            return float("nan")
        difference = a - b
        if len(difference) < 2 or np.isclose(np.ptp(difference), 0.0, atol=1e-12, rtol=0):
            return float("nan")
        value = stats.ttest_rel(a, b, nan_policy="omit").pvalue
    else:
        if np.isclose(np.ptp(a), 0.0, atol=1e-12, rtol=0) or np.isclose(np.ptp(b), 0.0, atol=1e-12, rtol=0):
            return float("nan")
        value = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue
    value = float(value)
    return value if np.isfinite(value) and 0 <= value <= 1 else float("nan")


def _validate_liana_coverage(
    frame: pd.DataFrame,
    resource_pairs: set[Tuple[str, str]],
    evaluable_pairs: set[Tuple[str, str]],
    cell_types: set[str],
    group_key: str,
) -> Tuple[set[Tuple[str, str]], set[Tuple[str, str]], set[Tuple[str, str]]]:
    observed_pairs = _pair_set(frame)
    missing_evaluable = evaluable_pairs - observed_pairs
    missing_non_evaluable = (resource_pairs - evaluable_pairs) - observed_pairs
    unexpected = observed_pairs - resource_pairs
    if missing_evaluable or unexpected:
        details = []
        if missing_evaluable:
            details.append(f"missing evaluable {sorted(missing_evaluable)!r}")
        if unexpected:
            details.append(f"unexpected {sorted(unexpected)!r}")
        raise LianaCommunicationInputError(
            f"LIANA output for donor-condition group {group_key!r} does not cover the local resource: "
            + "; ".join(details)
        )
    observed_types = set(frame["source"].astype(str)) | set(frame["target"].astype(str))
    unknown_types = observed_types - cell_types
    if unknown_types:
        raise LianaCommunicationInputError(
            f"LIANA output for donor-condition group {group_key!r} contains cell types not present in that group: "
            f"{sorted(unknown_types)!r}"
        )
    return observed_pairs, missing_non_evaluable, unexpected


class LianaCommunicationCapability(BaseCapability):
    """Run LIANA+'s RRA rank aggregate independently per donor-condition."""

    contract_operations = [
        "validate_local_ligand_receptor_resource",
        "validate_donor_condition_partition",
        "convert_scdata_to_anndata_per_donor_condition",
        "run_liana_rank_aggregate",
        "retain_all_ligand_receptor_pairs",
        "materialize_donor_level_rank_table",
        "compare_conditions_at_donor_level",
        "record_spatial_data_use",
        "record_resource_provenance",
    ]
    required_audit_checks = LIANA_AUDIT_REQUIREMENTS

    def __init__(self):
        super().__init__(
            capability_name="liana_communication",
            implementation_id="liana_rank_aggregate_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["scRNA", "spatial"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            output_types=[ArtifactType.TABLE],
            suitable_for=["donor_level_cell_cell_communication", "ligand_receptor_rank_aggregation"],
        )

    def _condition_design(
        self,
        data: SCData,
        donor_col: str,
        condition_col: str,
        params: Mapping[str, Any],
    ) -> Dict[str, Any]:
        obs = data.obs
        if obs[[donor_col, condition_col]].isna().any().any():
            raise LianaCommunicationInputError(
                "donor and condition metadata must be complete for every cell"
            )
        donors = obs[donor_col].astype(str).str.strip()
        conditions = obs[condition_col].astype(str).str.strip()
        if donors.eq("").any() or conditions.eq("").any():
            raise LianaCommunicationInputError("donor and condition labels must be non-empty")
        groups = [(str(d), str(c)) for d, c in sorted(set(zip(donors, conditions)))]
        if not groups:
            raise LianaCommunicationInputError("no donor-condition groups were observed")

        # A condition comparison is optional, but when requested it must be a
        # real donor-level design.  Same-donor cross-condition observations
        # are paired only when explicitly requested.
        observed_conditions = list(dict.fromkeys(conditions.tolist()))
        cond_a = params.get("condition_a")
        cond_b = params.get("condition_b")
        if cond_a is not None:
            cond_a = str(cond_a).strip()
        if cond_b is not None:
            cond_b = str(cond_b).strip()
        if cond_a is None and cond_b is None and len(observed_conditions) == 2:
            cond_a, cond_b = observed_conditions
        elif (cond_a is None) != (cond_b is None):
            raise LianaCommunicationInputError(
                "condition_a and condition_b must be supplied together"
            )
        if cond_a is not None and cond_a == cond_b:
            raise LianaCommunicationInputError("condition_a and condition_b must differ")
        if cond_a is not None and (cond_a not in observed_conditions or cond_b not in observed_conditions):
            raise LianaCommunicationInputError(
                f"requested conditions {cond_a!r}, {cond_b!r} are not both observed"
            )
        donors_a = set(donors[conditions == cond_a]) if cond_a is not None else set()
        donors_b = set(donors[conditions == cond_b]) if cond_b is not None else set()
        paired = bool(params.get("paired", False))
        overlap = donors_a & donors_b
        if overlap and not paired:
            raise LianaCommunicationInputError(
                "the same donor occurs in both comparison conditions; set paired=True "
                "to request a donor-paired comparison"
            )
        min_donors = int(params.get("min_donors_for_comparison", 2))
        if min_donors < 2:
            raise LianaCommunicationInputError("min_donors_for_comparison must be at least 2")
        return {
            "donors": donors,
            "conditions": conditions,
            "groups": groups,
            "observed_conditions": observed_conditions,
            "condition_a": cond_a,
            "condition_b": cond_b,
            "donors_a": donors_a,
            "donors_b": donors_b,
            "paired": paired,
            "min_donors": min_donors,
        }

    @staticmethod
    def _comparison_table(
        rows: pd.DataFrame,
        design: Mapping[str, Any],
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        columns = [
            "source", "target", "ligand", "receptor",
            "condition_a", "condition_b", "n_donors_a", "n_donors_b",
            "comparison_test", "comparison_status",
            "mean_magnitude_rank_a", "mean_magnitude_rank_b",
            "mean_specificity_rank_a", "mean_specificity_rank_b",
            "comparison_p_value_magnitude", "comparison_p_value_specificity",
        ]
        empty = pd.DataFrame(columns=columns)
        cond_a = design.get("condition_a")
        cond_b = design.get("condition_b")
        if cond_a is None or cond_b is None:
            return empty, {
                "comparison_status": "not_requested",
                "comparison_test": None,
                "comparison_statistical_unit": "donor",
            }

        keys = ["source", "target", "ligand", "receptor"]
        result_rows = []
        paired = bool(design["paired"])
        min_donors = int(design["min_donors"])
        donors_a = set(design["donors_a"])
        donors_b = set(design["donors_b"])
        overlap = donors_a & donors_b
        if paired:
            usable = overlap
        else:
            usable = None

        for key, group in rows.groupby(keys, dropna=False, sort=True):
            group = group.copy()
            # A donor can have multiple source/target rows for one LR pair;
            # average only within that donor-condition before comparing.
            donor_rows = (
                group.groupby(["donor_id", "condition"], as_index=False)[
                    ["magnitude_rank", "specificity_rank"]
                ]
                .mean()
            )
            a = donor_rows[donor_rows["condition"] == cond_a].set_index("donor_id")
            b = donor_rows[donor_rows["condition"] == cond_b].set_index("donor_id")
            test_name = "paired_t_test" if paired else "welch_t_test"
            if paired:
                common = sorted(set(a.index) & set(b.index) & usable)
                a, b = a.loc[common], b.loc[common]
                n_a = n_b = len(common)
            else:
                common = []
                n_a, n_b = len(a), len(b)
            status = "tested" if (n_a >= min_donors and n_b >= min_donors) else "insufficient_donors"
            p_mag = np.nan
            p_spec = np.nan
            if status == "tested":
                if paired:
                    va = a.loc[common, "magnitude_rank"].to_numpy(float)
                    vb = b.loc[common, "magnitude_rank"].to_numpy(float)
                    sa = a.loc[common, "specificity_rank"].to_numpy(float)
                    sb = b.loc[common, "specificity_rank"].to_numpy(float)
                    p_mag = _safe_rank_comparison_pvalue(va, vb, paired=True)
                    p_spec = _safe_rank_comparison_pvalue(sa, sb, paired=True)
                else:
                    p_mag = _safe_rank_comparison_pvalue(
                        a["magnitude_rank"].to_numpy(float),
                        b["magnitude_rank"].to_numpy(float),
                        paired=False,
                    )
                    p_spec = _safe_rank_comparison_pvalue(
                        a["specificity_rank"].to_numpy(float),
                        b["specificity_rank"].to_numpy(float),
                        paired=False,
                    )
                if not np.isfinite([p_mag, p_spec]).all():
                    status = "not_estimated"
            result_rows.append({
                "source": key[0], "target": key[1], "ligand": key[2], "receptor": key[3],
                "condition_a": cond_a, "condition_b": cond_b,
                "n_donors_a": n_a, "n_donors_b": n_b,
                "comparison_test": test_name, "comparison_status": status,
                "mean_magnitude_rank_a": float(a["magnitude_rank"].mean()) if n_a else np.nan,
                "mean_magnitude_rank_b": float(b["magnitude_rank"].mean()) if n_b else np.nan,
                "mean_specificity_rank_a": float(a["specificity_rank"].mean()) if n_a else np.nan,
                "mean_specificity_rank_b": float(b["specificity_rank"].mean()) if n_b else np.nan,
                "comparison_p_value_magnitude": p_mag,
                "comparison_p_value_specificity": p_spec,
            })
        comparison = pd.DataFrame(result_rows, columns=columns)
        # Each score defines a separate, explicitly named testing family.
        for score in ("magnitude", "specificity"):
            p_column = f"comparison_p_value_{score}"
            q_column = f"comparison_fdr_{score}"
            comparison[q_column] = np.nan
            valid = comparison[p_column].notna()
            if valid.any():
                comparison.loc[valid, q_column] = stats.false_discovery_control(
                    comparison.loc[valid, p_column].to_numpy(), method="bh")
        status = "tested" if not comparison.empty and (comparison["comparison_status"] == "tested").any() else "insufficient_donors"
        return comparison, {
            "comparison_status": status,
            "comparison_test": "paired_t_test" if paired else "welch_t_test",
            "comparison_statistical_unit": "donor",
            "comparison_fdr_applied": True,
            "comparison_fdr_family": "all estimated interactions, separately for each rank score",
            "comparison_min_donors": min_donors,
        }

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        if len(contract.input_artifacts) != 1:
            raise LianaCommunicationInputError(
                "liana_communication requires exactly one input AnnData artifact"
            )
        in_uri = ArtifactURI.parse(contract.input_artifacts[0]).to_string()
        meta, payload = registry.get(in_uri)
        if meta.type not in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
            raise LianaCommunicationInputError(
                f"liana_communication requires AnnData/SpatialData input, received {meta.type}"
            )
        data = _as_scdata(payload)
        if data.n_obs == 0 or data.n_vars == 0:
            raise LianaCommunicationInputError("communication input must contain cells and genes")
        params = contract.parameters
        species = str(params.get("species", "")).strip()
        resource_version = str(params.get("lr_resource_version", params.get("resource_version", ""))).strip()
        resource_source = str(params.get("lr_resource_source", params.get("resource_source", ""))).strip()
        if not species or not resource_version or not resource_source:
            raise LianaCommunicationInputError(
                "species, lr_resource_version, and lr_resource_source are required explicitly"
            )
        resource_path = _resource_path(params)
        resource, resource_meta = _read_resource(resource_path)
        _validate_external_resource_hash(params, resource_meta["sha256"])
        resource_species = _validate_resource_species(species, resource_meta, params)
        if "return_all_lrs" in params and params["return_all_lrs"] is not True:
            raise LianaCommunicationInputError(
                "liana_communication always requires return_all_lrs=True; significance filtering is not allowed"
            )
        li = _load_liana()

        donor_col = _resolve_column(data.obs, params.get("donor_col"), ("donor_id", "donor", "mouse_id", "sample_id", "sample"), "donor")
        condition_col = _resolve_column(data.obs, params.get("condition_col"), ("condition", "experimental_condition", "group"), "condition")
        cell_type_col = str(params.get("cell_type_col", "cell_type"))
        if cell_type_col not in data.obs.columns:
            raise LianaCommunicationInputError(
                f"cell type metadata column {cell_type_col!r} is absent"
            )
        cell_types = data.obs[cell_type_col].astype(str).str.strip()
        if (
            data.obs[cell_type_col].isna().any()
            or cell_types.eq("").any()
        ):
            raise LianaCommunicationInputError(
                "cell type metadata must be complete and contain non-empty labels"
            )
        design = self._condition_design(data, donor_col, condition_col, params)
        min_cells = int(params.get("min_cells", 5))
        if min_cells < 1:
            raise LianaCommunicationInputError("min_cells must be positive")
        gene_names = (
            data.var["gene_name"].astype(str).tolist()
            if "gene_name" in data.var.columns
            else [str(value) for value in data.var.index.tolist()]
        )
        if len(gene_names) != data.n_vars or len(set(gene_names)) != len(gene_names):
            raise LianaCommunicationInputError("gene names must be present and unique")
        resource_pairs = _pair_set(resource)
        gene_name_set = {str(name).strip() for name in gene_names}
        usable_resource_pairs = {
            pair for pair in resource_pairs if _pair_is_gene_compatible(pair, gene_name_set)
        }
        if not usable_resource_pairs:
            raise LianaCommunicationInputError(
                "the local ligand/receptor resource has no pair whose ligand and receptor "
                "genes are present in the input assay"
            )

        requested_spatial = params.get("spatial_key")
        explicit_use_spatial = params.get("use_spatial")
        if requested_spatial is None and explicit_use_spatial is not False and "spatial" in data.obsm:
            requested_spatial = "spatial"
        if requested_spatial is not None and explicit_use_spatial is False:
            raise LianaCommunicationInputError(
                "spatial_key was supplied (or spatial coordinates are present), but use_spatial=False"
            )
        if meta.type == ArtifactType.SPATIAL_DATA and requested_spatial is None:
            raise LianaCommunicationInputError(
                "SpatialData input requires explicit spatial coordinates in obsm['spatial'] "
                "or a supplied spatial_key"
            )
        use_spatial = bool(explicit_use_spatial) if explicit_use_spatial is not None else requested_spatial is not None
        if use_spatial and requested_spatial is None:
            raise LianaCommunicationInputError(
                "use_spatial=True requires an explicit spatial_key or obsm['spatial']"
            )
        if requested_spatial is not None:
            if requested_spatial not in data.obsm:
                raise LianaCommunicationInputError(
                    f"spatial_key {requested_spatial!r} is not present in input obsm"
                )
            coords = _as_dense(data.obsm[requested_spatial])
            if coords.ndim != 2 or coords.shape[0] != data.n_obs or coords.shape[1] < 2 or not np.isfinite(coords).all():
                raise LianaCommunicationInputError(
                    "spatial coordinates must be finite with shape (n_cells, >=2)"
                )
        else:
            use_spatial = False

        rank_params: Dict[str, Any] = {
            "expr_prop": float(params.get("expr_prop", 0.05)),
            "min_cells": min_cells,
            "aggregate_method": str(params.get("aggregate_method", "rra")),
            "return_all_lrs": True,
            "use_raw": bool(params.get("use_raw", False)),
            "de_method": str(params.get("de_method", "t-test")),
            "n_perms": params.get("n_perms", 1000),
            "seed": int(params.get("random_seed", params.get("seed", 1337))),
            "n_jobs": int(params.get("n_jobs", 1)),
            "verbose": bool(params.get("verbose", False)),
        }
        if rank_params["aggregate_method"] not in {"rra", "mean"}:
            raise LianaCommunicationInputError("aggregate_method must be 'rra' or 'mean'")
        if rank_params["n_perms"] is not None:
            rank_params["n_perms"] = int(rank_params["n_perms"])
            if rank_params["n_perms"] < 0:
                raise LianaCommunicationInputError("n_perms must be non-negative or None")
        layer = params.get("layer")
        if layer is not None:
            if layer not in data.layers:
                raise LianaCommunicationInputError(f"requested expression layer {layer!r} is absent")
            rank_params["layer"] = str(layer)
        groupby_pairs = params.get("groupby_pairs")
        if groupby_pairs is not None:
            groupby_pairs = pd.DataFrame(groupby_pairs) if not isinstance(groupby_pairs, pd.DataFrame) else groupby_pairs.copy()
            if not {"source", "target"}.issubset(groupby_pairs.columns):
                raise LianaCommunicationInputError("groupby_pairs must contain source and target columns")
            rank_params["groupby_pairs"] = groupby_pairs
        if use_spatial:
            rank_params["spatial_key"] = str(requested_spatial)
            if params.get("spatial_kwargs") is not None:
                if not isinstance(params["spatial_kwargs"], Mapping):
                    raise LianaCommunicationInputError("spatial_kwargs must be a mapping")
                rank_params["spatial_kwargs"] = dict(params["spatial_kwargs"])

        result_frames = []
        group_cell_counts: Dict[str, int] = {}
        group_pair_coverage: Dict[str, int] = {}
        group_missing_non_evaluable: Dict[str, list[Dict[str, Any]]] = {}
        rank_aggregate = li.mt.rank_aggregate
        for donor, condition in design["groups"]:
            mask = (design["donors"] == donor) & (design["conditions"] == condition)
            group_data = data.subset_obs(mask.to_numpy())
            group_key = f"{donor}::{condition}"
            group_cell_counts[group_key] = int(group_data.n_obs)
            if group_data.n_obs < min_cells:
                raise LianaCommunicationInputError(
                    f"donor-condition group {group_key!r} has {group_data.n_obs} cells; "
                    f"min_cells={min_cells}"
                )
            group_cell_types = group_data.obs[cell_type_col].astype(str).str.strip()
            counts = group_cell_types.value_counts()
            if len(counts) < 2 or bool((counts < min_cells).any()):
                raise LianaCommunicationInputError(
                    f"donor-condition group {group_key!r} requires at least two cell types "
                    f"with {min_cells} cells each"
                )
            adata = _prepare_adata(group_data, gene_names)
            call_params = dict(rank_params)
            # Avoid passing a pandas object by reference into LIANA's mutable
            # preprocessing path on each independent donor invocation.
            if isinstance(call_params.get("groupby_pairs"), pd.DataFrame):
                call_params["groupby_pairs"] = call_params["groupby_pairs"].copy()
            raw_result = rank_aggregate(
                adata,
                groupby=cell_type_col,
                resource=resource.copy(),
                inplace=False,
                **call_params,
            )
            frame = _frame_from_liana(raw_result)
            # A resource can contain valid interactions whose genes are absent
            # from this assay.  Preserve the LIANA table for estimable rows,
            # while recording those omitted/unestimated pairs explicitly.
            rank_values = frame[["magnitude_rank", "specificity_rank"]].to_numpy(dtype=float)
            nonfinite_rows = ~np.isfinite(rank_values).all(axis=1)
            if nonfinite_rows.any():
                nonfinite_pairs = _pair_set(frame.loc[nonfinite_rows])
                invalid_evaluable = nonfinite_pairs & usable_resource_pairs
                if invalid_evaluable:
                    raise LianaCommunicationInputError(
                        f"LIANA returned non-estimable ranks for gene-compatible pairs in group {group_key!r}: "
                        f"{sorted(invalid_evaluable)!r}"
                    )
                frame = frame.loc[~nonfinite_rows].copy()
            covered_pairs, missing_non_evaluable, _ = _validate_liana_coverage(
                frame,
                resource_pairs,
                usable_resource_pairs,
                set(group_cell_types),
                group_key,
            )
            group_pair_coverage[group_key] = len(covered_pairs)
            group_missing_non_evaluable[group_key] = _missing_pair_records(
                missing_non_evaluable,
                gene_name_set,
            )
            if frame.empty:
                raise LianaCommunicationInputError(
                    f"LIANA returned no estimable ligand/receptor rows for donor-condition group {group_key!r}"
                )
            frame.insert(0, "donor_id", donor)
            frame.insert(1, "condition", condition)
            frame.insert(2, "donor_condition", group_key)
            frame.insert(3, "n_cells_group", int(group_data.n_obs))
            result_frames.append(frame)

        full_table = pd.concat(result_frames, ignore_index=True)
        full_table["species"] = species
        full_table["resource_version"] = resource_version
        full_table["resource_source"] = resource_source
        full_table["resource_sha256"] = resource_meta["sha256"]
        full_table["rank_method"] = rank_params["aggregate_method"]
        full_table["rank_is_fdr"] = False
        full_table["spatial_used"] = bool(use_spatial)
        comparison, comparison_metrics = self._comparison_table(full_table, design)

        study_id = ArtifactURI.parse(in_uri).study_id
        primary_uri = _output_uri(contract, 0, study_id)
        comparison_uri = _output_uri(contract, 1, study_id)
        provenance = {
            "resource_path": resource_meta["path"],
            "resource_sha256": resource_meta["sha256"],
            "resource_version": resource_version,
            "resource_source": resource_source,
            "species": species,
            "resource_species": resource_species,
            "resource_species_columns": list(resource_meta.get("species_columns", [])),
            "resource_species_values": list(resource_meta.get("species_values", [])),
        }
        effective_parameters = {
            "species": species,
            "donor_col": donor_col,
            "condition_col": condition_col,
            "cell_type_col": cell_type_col,
            "min_cells": min_cells,
            "condition_a": design["condition_a"],
            "condition_b": design["condition_b"],
            "paired": design["paired"],
            "min_donors_for_comparison": design["min_donors"],
            "spatial_key": requested_spatial if use_spatial else None,
            "resource": provenance,
            "rank_aggregate": {k: v for k, v in rank_params.items() if not isinstance(v, pd.DataFrame)},
        }
        if params.get("external_resource_sha256") is not None:
            effective_parameters["external_resource_sha256"] = dict(params["external_resource_sha256"])
        if params.get("lr_resource_sha256", params.get("resource_sha256")) is not None:
            effective_parameters["lr_resource_sha256"] = str(
                params.get("lr_resource_sha256", params.get("resource_sha256"))
            )
        if resource_species is not None:
            effective_parameters["lr_resource_species"] = resource_species
        summary_metrics: Dict[str, Any] = {
            "n_rank_rows": int(len(full_table)),
            "n_unique_interactions": int(full_table[["source", "target", "ligand", "receptor"]].drop_duplicates().shape[0]),
            "donor_condition_groups": list(group_cell_counts),
            "group_cell_counts": group_cell_counts,
            "resource_n_pairs": int(len(resource_pairs)),
            "resource_usable_n_pairs": int(len(usable_resource_pairs)),
            "output_pairs_per_group": group_pair_coverage,
            "resource_pair_coverage": {
                "expected_pairs": int(len(resource_pairs)),
                "expected_evaluable_pairs": int(len(usable_resource_pairs)),
                "covered_per_group": group_pair_coverage,
                "complete": all(value == len(resource_pairs) for value in group_pair_coverage.values()),
                "evaluable_complete": all(
                    value == len(usable_resource_pairs) for value in group_pair_coverage.values()
                ),
            },
            "group_missing_non_evaluable": group_missing_non_evaluable,
            "donors_by_condition": {condition: sorted(set(design["donors"][design["conditions"] == condition])) for condition in design["observed_conditions"]},
            "all_ligand_receptor_pairs_retained": all(
                value == len(resource_pairs) for value in group_pair_coverage.values()
            ),
            "all_evaluable_ligand_receptor_pairs_retained": all(
                value == len(usable_resource_pairs) for value in group_pair_coverage.values()
            ),
            "rank_columns": ["magnitude_rank", "specificity_rank"],
            "ranking_not_fdr": True,
            "spatial_used": bool(use_spatial),
            "spatial_key": requested_spatial if use_spatial else None,
            "resource_provenance": provenance,
            "required_audit_checks": list(self.required_audit_checks),
            **comparison_metrics,
        }
        software_versions = {"liana": _package_version("liana"), "anndata": _package_version("anndata")}
        registry.register(
            uri_str=primary_uri,
            payload=full_table,
            artifact_type=ArtifactType.TABLE,
            study_id=study_id,
            created_by_task=contract.task_id,
            operation="liana_rank_aggregate_donor_condition",
            parent_uris=[in_uri],
            parameters=effective_parameters,
            software_versions=software_versions,
            random_seed=rank_params["seed"],
            summary_metrics=summary_metrics,
        )
        output_artifacts = [primary_uri]
        if comparison_uri is not None:
            registry.register(
                uri_str=comparison_uri,
                payload=comparison,
                artifact_type=ArtifactType.TABLE,
                study_id=study_id,
                created_by_task=contract.task_id,
                operation="liana_rank_condition_comparison_donor_level",
                parent_uris=[primary_uri],
                parameters=effective_parameters,
                software_versions=software_versions,
                random_seed=rank_params["seed"],
                summary_metrics=summary_metrics,
            )
            output_artifacts.append(comparison_uri)

        operations = [
            "validate_local_ligand_receptor_resource",
            "validate_donor_condition_partition",
            "convert_scdata_to_anndata_per_donor_condition",
            "run_liana_rank_aggregate",
            "retain_all_ligand_receptor_pairs",
            "materialize_donor_level_rank_table",
        ]
        if design["condition_a"] is not None:
            operations.append("compare_conditions_at_donor_level")
        if use_spatial:
            operations.append("record_spatial_data_use")
        operations.append("record_resource_provenance")
        if contract.allowed_operations:
            # An explicitly narrowed contract must expose every operation it
            # authorizes.  Silent omission would make provenance misleading.
            missing_ops = [op for op in operations if op not in contract.allowed_operations]
            if missing_ops:
                raise LianaCommunicationInputError(
                    f"contract.allowed_operations omits operations actually required by LIANA: {missing_ops}"
                )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=output_artifacts,
            executed_operations=operations,
            metrics=summary_metrics,
        )


__all__ = [
    "LIANA_AUDIT_REQUIREMENTS",
    "LianaCommunicationCapability",
    "LianaCommunicationDependencyError",
    "LianaCommunicationInputError",
]
