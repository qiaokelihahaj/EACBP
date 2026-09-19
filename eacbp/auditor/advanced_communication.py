"""Independent integrity checks for LIANA donor ranks and comparisons.

The executor records provenance and donor strata, but this auditor reads the
local resource and input artifact again.  It therefore does not certify a
result merely because the executor set ``return_all_lrs=True`` or copied its
own summary metrics.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData


_SPECIES_ALIASES = {
    "human": {"human", "homo_sapiens", "homo sapiens", "9606"},
    "mouse": {"mouse", "mus_musculus", "mus musculus", "10090"},
    "rat": {"rat", "rattus_norvegicus", "rattus norvegicus", "10116"},
    "zebrafish": {"zebrafish", "danio_rerio", "danio rerio", "7955"},
}
_RESOURCE_SPECIES_COLUMNS = {
    "species",
    "organism",
    "species_a",
    "species_b",
    "ligand_species",
    "receptor_species",
    "source_species",
    "target_species",
}


def _normalise_species(value: Any) -> str:
    return re.sub(r"[\s\-]+", "_", str(value).strip().lower())


def _species_compatible(expected: Any, observed: Any) -> bool:
    left = _normalise_species(expected)
    right = _normalise_species(observed)
    if left == right:
        return True
    for aliases in _SPECIES_ALIASES.values():
        normalised = {_normalise_species(value) for value in aliases}
        if left in normalised and right in normalised:
            return True
    return False


def _pair_set(frame: pd.DataFrame) -> set[tuple[str, str]]:
    return {
        (str(ligand).strip(), str(receptor).strip())
        for ligand, receptor in frame[["ligand", "receptor"]].itertuples(index=False, name=None)
    }


def _complex_components(value: Any) -> set[str]:
    label = str(value).strip()
    return {token for token in re.split(r"[_+:|]", label) if token}


def _pair_is_gene_compatible(pair: tuple[str, str], gene_names: set[str]) -> bool:
    return all(
        component in gene_names
        for member in pair
        for component in _complex_components(member)
    )


def _missing_pair_records(
    pairs: set[tuple[str, str]],
    gene_names: set[str],
) -> list[dict[str, Any]]:
    """Describe non-evaluable resource pairs from the independent inputs."""

    records: list[dict[str, Any]] = []
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


def _canonical_missing_records(value: Any) -> list[tuple[str, str, tuple[str, ...]]] | None:
    """Normalize recorded missing-pair evidence for exact comparison."""

    if not isinstance(value, (list, tuple)):
        return None
    records: list[tuple[str, str, tuple[str, ...]]] = []
    for item in value:
        if not isinstance(item, Mapping):
            return None
        if "ligand" not in item or "receptor" not in item or "missing_genes" not in item:
            return None
        missing_genes = item["missing_genes"]
        if not isinstance(missing_genes, (list, tuple, set)):
            return None
        records.append(
            (
                str(item["ligand"]).strip(),
                str(item["receptor"]).strip(),
                tuple(sorted(str(value).strip() for value in missing_genes)),
            )
        )
    return sorted(records)


def _read_resource(path: Path) -> tuple[str, set[tuple[str, str]], list[str], list[str]]:
    """Read the local resource independently and return its immutable facts."""

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
        raise ValueError(f"unsupported local ligand/receptor resource format {suffix!r}")
    columns = {str(column).strip().lower(): column for column in frame.columns}
    ligand_col = columns.get("ligand")
    receptor_col = columns.get("receptor")
    if ligand_col is None or receptor_col is None:
        raise ValueError("local resource must contain ligand and receptor columns")
    frame = frame.rename(columns={ligand_col: "ligand", receptor_col: "receptor"}).copy()
    frame["ligand"] = frame["ligand"].astype("string").str.strip()
    frame["receptor"] = frame["receptor"].astype("string").str.strip()
    frame = frame.loc[
        frame["ligand"].notna()
        & frame["receptor"].notna()
        & frame["ligand"].ne("")
        & frame["receptor"].ne("")
    ].drop_duplicates(["ligand", "receptor"], keep="first")
    if frame.empty:
        raise ValueError("local resource contains no valid ligand/receptor pairs")
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
    return digest, _pair_set(frame), species_values, species_columns


def _safe_rank_comparison_pvalue(
    values_a: Sequence[float],
    values_b: Sequence[float],
    *,
    paired: bool,
) -> float:
    """Return a p-value only when the donor-level test is estimable.

    SciPy may return a finite-looking value while warning about a zero
    variance sample.  The capability treats that situation as
    ``not_estimated``; the independent auditor applies the same mathematical
    preconditions without importing the executor's helper.
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
        p_value = stats.ttest_rel(a, b, nan_policy="omit").pvalue
    else:
        if np.isclose(np.ptp(a), 0.0, atol=1e-12, rtol=0) or np.isclose(np.ptp(b), 0.0, atol=1e-12, rtol=0):
            return float("nan")
        p_value = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue
    p_value = float(p_value)
    return p_value if np.isfinite(p_value) and 0 <= p_value <= 1 else float("nan")


def _check_pinned_hashes(parameters: Mapping[str, Any], digest: str, errors: list[str]) -> None:
    pinned = parameters.get("external_resource_sha256")
    if pinned is not None and not isinstance(pinned, Mapping):
        errors.append("external_resource_sha256 is not a mapping")
    if isinstance(pinned, Mapping):
        found = False
        for key in ("lr_resource_path", "resource_path", "ligand_receptor_resource_path"):
            if key in pinned:
                found = True
                value = str(pinned[key]).strip().lower()
                if not re.fullmatch(r"[0-9a-f]{64}", value) or value != digest.lower():
                    errors.append("Pinned ligand/receptor resource hash differs from the local file")
        if not found and pinned:
            errors.append("external_resource_sha256 has no ligand/receptor resource key")
    direct = parameters.get("lr_resource_sha256", parameters.get("resource_sha256"))
    if direct is not None:
        value = str(direct).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", value) or value != digest.lower():
            errors.append("Declared ligand/receptor resource hash differs from the local file")


class LianaValidator(BaseAuditor):
    """Recompute provenance, coverage, grouping and donor-level statistics."""

    def __init__(self):
        super().__init__("liana_validator")

    def audit(self, contract, result, registry: ArtifactRegistry):
        report = ValidationReport(auditor_name=self.auditor_name, target_task_id=contract.task_id)
        if contract.capability != "liana_communication":
            return report
        errors: list[str] = []
        try:
            if len(contract.input_artifacts) != 1 or not result.output_artifacts:
                raise ValueError("LIANA result must have one input and at least one output artifact")
            input_meta, input_payload = registry.get(contract.input_artifacts[0])
            output_meta, output_payload = registry.get(result.output_artifacts[0])
            data = input_payload if isinstance(input_payload, SCData) else SCData.from_dict(input_payload)
            ranks = output_payload if isinstance(output_payload, pd.DataFrame) else pd.DataFrame(output_payload)
            parameters = output_meta.parameters
            summary_metrics = output_meta.summary_metrics
            resource = parameters.get("resource")
            if not isinstance(resource, Mapping):
                raise ValueError("LIANA output metadata has no resource provenance")
            resource_path = Path(str(resource["resource_path"])).expanduser().resolve()
            digest, resource_pairs, resource_species_values, resource_species_columns = _read_resource(resource_path)
            if str(resource.get("resource_sha256", "")).lower() != digest.lower():
                errors.append("Resource hash does not match local resource")
            if str(resource.get("resource_path", "")) != str(resource_path):
                errors.append("Resource path is not recorded canonically")
            if set(resource.get("resource_species_values", [])) != set(resource_species_values):
                errors.append("Recorded resource species values differ from local resource")
            if set(resource.get("resource_species_columns", [])) != set(resource_species_columns):
                errors.append("Recorded resource species columns differ from local resource")
            _check_pinned_hashes(contract.parameters, digest, errors)
            _check_pinned_hashes(parameters, digest, errors)

            expected_species = str(contract.parameters.get("species", parameters.get("species", ""))).strip()
            if not expected_species or str(parameters.get("species", "")).strip() != expected_species:
                errors.append("LIANA species provenance does not match the contract")
            for value in resource_species_values:
                if not _species_compatible(expected_species, value):
                    errors.append(f"Resource species {value!r} does not match {expected_species!r}")
            declared_resource_species = parameters.get("lr_resource_species", resource.get("resource_species"))
            if declared_resource_species and not _species_compatible(expected_species, declared_resource_species):
                errors.append("Declared resource species does not match the LIANA species")
            for key in ("resource_version", "resource_source"):
                if not str(resource.get(key, "")).strip():
                    errors.append(f"Missing resource provenance: {key}")

            required_rank_columns = {
                "donor_id", "condition", "donor_condition", "n_cells_group",
                "source", "target", "ligand", "receptor", "magnitude_rank",
                "specificity_rank", "species", "resource_version", "resource_source",
                "resource_sha256", "rank_is_fdr", "spatial_used",
            }
            missing = required_rank_columns - set(ranks.columns)
            if missing:
                errors.append(f"LIANA rank table is missing columns: {sorted(missing)}")
                raise ValueError("cannot continue rank integrity checks with missing columns")
            for key in ("species", "resource_version", "resource_source", "resource_sha256"):
                expected = resource.get(key)
                if not ranks[key].astype(str).eq(str(expected)).all():
                    errors.append(f"Rank table provenance differs from resource metadata: {key}")
            if not ranks.rank_is_fdr.eq(False).all():
                errors.append("Communication ranks must not be represented as FDR")
            for column in ("magnitude_rank", "specificity_rank"):
                values = pd.to_numeric(ranks[column], errors="coerce")
                numeric = values.to_numpy(dtype=float)
                if not np.isfinite(numeric).all() or not values.between(0, 1).all():
                    errors.append(f"Invalid rank: {column}")
            rank_aggregate = parameters.get("rank_aggregate", {})
            if not isinstance(rank_aggregate, Mapping) or rank_aggregate.get("return_all_lrs") is not True:
                errors.append("LIANA rank_aggregate was not configured with return_all_lrs=True")

            donor_col = str(parameters["donor_col"])
            condition_col = str(parameters["condition_col"])
            cell_type_col = str(parameters["cell_type_col"])
            for column in (donor_col, condition_col, cell_type_col):
                if column not in data.obs.columns:
                    raise ValueError(f"Input metadata column {column!r} is absent")
            donor_values = data.obs[donor_col].astype(str).str.strip()
            condition_values = data.obs[condition_col].astype(str).str.strip()
            cell_type_values = data.obs[cell_type_col].astype(str).str.strip()
            expected_groups = {
                (str(donor), str(condition)): int(size)
                for (donor, condition), size in data.obs.groupby([donor_col, condition_col], observed=True).size().items()
            }
            actual_groups = set(zip(ranks.donor_id.astype(str), ranks.condition.astype(str)))
            if actual_groups != set(expected_groups):
                errors.append("Rank table does not cover exactly the input donor-condition groups")
            key_columns = ["source", "target", "ligand", "receptor"]
            if ranks.duplicated(["donor_id", "condition", *key_columns]).any():
                errors.append("Duplicate donor interaction observations")
            gene_names = (
                data.var["gene_name"].astype(str).tolist()
                if "gene_name" in data.var.columns
                else [str(value) for value in data.var.index.tolist()]
            )
            gene_name_set = {str(value).strip() for value in gene_names}
            evaluable_pairs = {
                pair for pair in resource_pairs if _pair_is_gene_compatible(pair, gene_name_set)
            }
            if not evaluable_pairs:
                errors.append("No local resource pair is gene-compatible with the input assay")
            non_evaluable_pairs = resource_pairs - evaluable_pairs
            recorded_missing = summary_metrics.get("group_missing_non_evaluable")
            if non_evaluable_pairs and recorded_missing is None:
                errors.append(
                    "LIANA output metadata is missing group_missing_non_evaluable evidence"
                )
            elif recorded_missing is not None and not isinstance(recorded_missing, Mapping):
                errors.append("group_missing_non_evaluable metadata is not a mapping")
            for group, expected_count in expected_groups.items():
                donor, condition = group
                frame = ranks[
                    (ranks.donor_id.astype(str) == donor)
                    & (ranks.condition.astype(str) == condition)
                ]
                observed_pairs = _pair_set(frame) if not frame.empty else set()
                expected_missing = _missing_pair_records(
                    non_evaluable_pairs - observed_pairs,
                    gene_name_set,
                )
                if isinstance(recorded_missing, Mapping):
                    group_key = f"{donor}::{condition}"
                    if group_key not in recorded_missing:
                        errors.append(
                            f"group_missing_non_evaluable has no record for {group_key}"
                        )
                    else:
                        recorded = _canonical_missing_records(recorded_missing[group_key])
                        expected = _canonical_missing_records(expected_missing)
                        if recorded is None or recorded != expected:
                            errors.append(
                                f"group_missing_non_evaluable differs for {group_key}"
                            )
                if frame.empty:
                    continue
                if not evaluable_pairs.issubset(observed_pairs):
                    errors.append(f"LIANA output omits an evaluable LR pair for {donor}::{condition}")
                if not observed_pairs.issubset(resource_pairs):
                    errors.append(f"LIANA output contains an LR pair outside the local resource for {donor}::{condition}")
                input_mask = (donor_values == donor) & (condition_values == condition)
                expected_types = set(cell_type_values[input_mask])
                observed_types = set(frame.source.astype(str)) | set(frame.target.astype(str))
                if not observed_types.issubset(expected_types):
                    errors.append(f"Rank output contains cell types outside {donor}::{condition}")
                reported_counts = pd.to_numeric(frame.n_cells_group, errors="coerce")
                if not reported_counts.eq(expected_count).all():
                    errors.append(f"Incorrect donor-condition cell count for {donor}::{condition}")

            spatial_key = parameters.get("spatial_key")
            spatial_used = ranks.spatial_used.astype(bool)
            if spatial_key is None and spatial_used.any():
                errors.append("Spatial ranks are marked used without a spatial_key")
            if spatial_key is not None:
                if spatial_key not in data.obsm:
                    errors.append("Spatial inference has no input coordinates")
                else:
                    coords = data.obsm[spatial_key]
                    if hasattr(coords, "toarray"):
                        coords = coords.toarray()
                    coords = np.asarray(coords)
                    if coords.ndim != 2 or coords.shape[0] != data.n_obs or coords.shape[1] < 2 or not np.isfinite(coords).all():
                        errors.append("Spatial coordinates are not finite with shape (n_cells, >=2)")
                if not spatial_used.all():
                    errors.append("Spatial ranks are not consistently marked as spatial")
            if input_meta.type.value == "SpatialData" and spatial_key is None:
                errors.append("SpatialData input was audited without spatial coordinates")

            condition_a = parameters.get("condition_a")
            condition_b = parameters.get("condition_b")
            comparison_requested = condition_a is not None or condition_b is not None
            if comparison_requested and (condition_a is None or condition_b is None):
                errors.append("Condition comparison metadata has only one condition")
            if comparison_requested and len(result.output_artifacts) < 2:
                errors.append("Requested condition comparison has no second output artifact")
            if not comparison_requested and len(result.output_artifacts) > 1:
                errors.append("Unexpected condition comparison output without a requested comparison")
            if comparison_requested and len(result.output_artifacts) >= 2:
                _, comparison_payload = registry.get(result.output_artifacts[1])
                comparison = comparison_payload if isinstance(comparison_payload, pd.DataFrame) else pd.DataFrame(comparison_payload)
                if comparison.empty:
                    errors.append("Requested condition comparison is empty")
                required_comparison = {
                    *key_columns,
                    "n_donors_a", "n_donors_b", "comparison_status",
                    "mean_magnitude_rank_a", "mean_magnitude_rank_b",
                    "mean_specificity_rank_a", "mean_specificity_rank_b",
                    "comparison_p_value_magnitude", "comparison_p_value_specificity",
                    "comparison_fdr_magnitude", "comparison_fdr_specificity",
                }
                missing_comparison = required_comparison - set(comparison.columns)
                if missing_comparison:
                    errors.append(f"Condition comparison is missing columns: {sorted(missing_comparison)}")
                elif not comparison.empty:
                    if comparison.duplicated(key_columns).any():
                        errors.append("Duplicate condition comparisons")
                    rank_pairs = set(map(tuple, ranks[key_columns].drop_duplicates().to_numpy()))
                    comparison_pairs = set(map(tuple, comparison[key_columns].drop_duplicates().to_numpy()))
                    if comparison_pairs != rank_pairs:
                        errors.append("Condition comparison does not cover all interactions")
                    paired = bool(parameters.get("paired", False))
                    min_donors = int(parameters.get("min_donors_for_comparison", 2))
                    for _, row in comparison.iterrows():
                        selected = ranks
                        for key in key_columns:
                            selected = selected[selected[key] == row[key]]
                        a = selected[selected.condition == condition_a].set_index("donor_id")
                        b = selected[selected.condition == condition_b].set_index("donor_id")
                        if paired:
                            common = sorted(set(a.index) & set(b.index))
                            a, b = a.loc[common], b.loc[common]
                        elif len(set(a.index) & set(b.index)):
                            errors.append("Unpaired comparison contains overlapping donors")
                        if int(row.n_donors_a) != len(a) or int(row.n_donors_b) != len(b):
                            errors.append("Comparison donor counts are incorrect")
                        enough = min(len(a), len(b)) >= min_donors
                        p_values = {}
                        for score in ("magnitude", "specificity"):
                            rank_column = f"{score}_rank"
                            mean_a = a[rank_column].mean() if len(a) else np.nan
                            mean_b = b[rank_column].mean() if len(b) else np.nan
                            if not np.isclose(row[f"mean_{rank_column}_a"], mean_a, equal_nan=True):
                                errors.append("Reported mean differs from donor observations")
                            if not np.isclose(row[f"mean_{rank_column}_b"], mean_b, equal_nan=True):
                                errors.append("Reported mean differs from donor observations")
                            p_value = float("nan")
                            if enough:
                                p_value = _safe_rank_comparison_pvalue(
                                    a[rank_column].to_numpy(float),
                                    b[rank_column].to_numpy(float),
                                    paired=paired,
                                )
                            p_values[score] = p_value
                            if not np.isclose(row[f"comparison_p_value_{score}"], p_value, equal_nan=True):
                                errors.append("Condition p-value differs from donor-level test")
                        status = (
                            "insufficient_donors"
                            if not enough
                            else ("tested" if np.isfinite(list(p_values.values())).all() else "not_estimated")
                        )
                        if row.comparison_status != status:
                            errors.append("Invalid condition comparison status")
                    for score in ("magnitude", "specificity"):
                        p_values = pd.to_numeric(comparison[f"comparison_p_value_{score}"], errors="coerce")
                        q_values = pd.to_numeric(comparison[f"comparison_fdr_{score}"], errors="coerce")
                        valid = p_values.notna()
                        if (p_values[valid] < 0).any() or (p_values[valid] > 1).any():
                            errors.append("Condition p-value is outside [0, 1]")
                        if q_values[~valid].notna().any():
                            errors.append("Condition FDR is present for an unestimated interaction")
                        if valid.any():
                            expected_q = stats.false_discovery_control(
                                p_values[valid].to_numpy(float), method="bh"
                            )
                            if not np.allclose(q_values[valid], expected_q):
                                errors.append("Condition comparison FDR is inconsistent with BH")
        except Exception as exc:
            errors.append(f"{type(exc).__name__}: {exc}")
        report.add_check(
            "liana_integrity",
            not errors,
            ValidationSeverity.ERROR,
            "; ".join(errors)
            if errors
            else "Verified local resource provenance, LR coverage, donor strata, spatial semantics and donor-level comparisons.",
        )
        return report
