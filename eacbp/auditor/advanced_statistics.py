"""Independent audits for the donor-level statistical capabilities.

The executor records useful provenance, but this module reconstructs the
input contrast and design from the source artifact and checks the result
tables independently. In particular, a task cannot certify a donor LOO run
merely by reporting that all fits completed.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm, t

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.capabilities.advanced_statistics import _coerce_data
from eacbp.capabilities.sc_data import SCData


METHODS = {
    "pydeseq2_pseudobulk_v1",
    "decoupler_ulm_v2",
    "pydeseq2_leave_one_donor_out_v1",
}

# Shared with the computational auditor. The advanced validator still checks
# every nullable value's status and reason, so this only prevents a legitimate
# unestimated value from being rejected as a generic non-finite error.
NULLABLE_STATISTICS = {
    "base_mean", "log2_fold_change", "lfc_se", "ci_low", "ci_high", "statistic",
    "p_value", "fdr_q_value", "padj", "decoupler_p_value",
    "mean_decoupler_p_value_condition_a", "mean_decoupler_p_value_condition_b",
    "status_reason", "skip_reason", "covariates", "network_path",
    "activity_se", "activity_ci_low", "activity_ci_high",
    "activity_effect_condition_a_vs_b", "mean_activity_condition_a",
    "mean_activity_condition_b", "baseline_log2_fold_change",
    "baseline_fdr_q_value", "min_log2_fold_change", "max_log2_fold_change",
    "median_log2_fold_change", "estimated_coverage", "direction_consistency",
    "significance_retention",
}


@dataclass
class _AuditInput:
    counts: np.ndarray
    genes: list[str]
    obs: pd.DataFrame
    metadata: pd.DataFrame
    donor_col: str
    condition_col: str
    condition_a: str
    condition_b: str
    paired: bool
    covariates: list[str]
    donors_a: list[str]
    donors_b: list[str]
    formula: str
    design: pd.DataFrame
    design_rank: int
    residual_df: int


def _dense(value: Any) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    elif hasattr(value, "A"):
        value = value.A
    return np.asarray(value)


def _resolve_column(obs: pd.DataFrame, explicit: Any, candidates: tuple[str, ...], label: str) -> str:
    if explicit is not None:
        value = str(explicit)
        if value not in obs.columns:
            raise ValueError(f"{label} metadata column {value!r} is absent")
        return value
    for value in candidates:
        if value in obs.columns:
            return value
    raise ValueError(f"{label} metadata column is absent")


def _raw_counts(data: SCData, params: Mapping[str, Any]) -> np.ndarray:
    layer = str(params.get("counts_layer", params.get("raw_counts_layer", "counts")))
    if layer in (getattr(data, "layers", {}) or {}):
        value = data.layers[layer]
    elif bool(params.get("allow_x_as_counts", False)):
        value = data.X
    else:
        raise ValueError(f"raw integer counts are absent from layers[{layer!r}]")
    arr = _dense(value)
    if arr.ndim != 2 or tuple(arr.shape) != tuple(data.shape):
        raise ValueError(f"counts shape {getattr(arr, 'shape', None)} does not match input {data.shape}")
    try:
        finite = np.isfinite(arr)
        negative = np.any(arr < 0)
        integer = np.all(np.equal(arr, np.floor(arr)))
    except (TypeError, ValueError) as exc:
        raise ValueError("counts must be numeric") from exc
    if not finite.all():
        raise ValueError("counts contain NaN or infinite values")
    if negative:
        raise ValueError("counts contain negative values")
    if not integer:
        raise ValueError("counts must be integer-valued")
    if np.any(arr > np.iinfo(np.int64).max):
        raise ValueError("counts exceed int64 range")
    return np.asarray(arr, dtype=np.int64)


def _gene_names(data: SCData) -> list[str]:
    if "gene_name" in data.var.columns:
        genes = data.var["gene_name"].astype(str).tolist()
    else:
        genes = [str(v) for v in data.var.index.tolist()]
    if len(genes) != data.n_vars or any(not g or g.lower() in {"nan", "none"} for g in genes):
        raise ValueError("gene metadata do not define a complete gene universe")
    if len(set(genes)) != len(genes):
        raise ValueError("gene names are not unique")
    return genes


def _covariates(obs: pd.DataFrame, params: Mapping[str, Any]) -> list[str]:
    raw = params.get("covariates", params.get("covariate_cols", []))
    if raw is None:
        return []
    if isinstance(raw, str):
        raw = [raw]
    if not isinstance(raw, (list, tuple, set)):
        raise ValueError("covariates must be a string or sequence")
    values: list[str] = []
    for item in raw:
        name = str(item)
        if name not in obs.columns:
            raise ValueError(f"covariate metadata column {name!r} is absent")
        if name not in values:
            values.append(name)
    return values


def _formula(params: Mapping[str, Any], paired: bool, covariates: list[str]) -> str:
    explicit = params.get("design_formula", params.get("design"))
    if explicit is not None:
        formula = str(explicit).strip()
        if not formula:
            raise ValueError("design formula is empty")
        if not formula.startswith("~"):
            formula = "~" + formula
        cleaned = re.sub(r"\bC\s*\(\s*([^)]*)\s*\)", r"\1", formula)
        variables = set(re.findall(r"[A-Za-z_]\w*", cleaned))
        extra: list[str] = []
        if paired and "donor" not in variables:
            extra.append("donor")
        extra.extend(c for c in covariates if c not in variables)
        if extra:
            formula += " + " + " + ".join(extra)
        return formula
    terms = ["condition"]
    if paired:
        terms.append("donor")
    terms.extend(covariates)
    return "~" + " + ".join(terms)


def _reconstruct_input(
    data: SCData,
    params: Mapping[str, Any],
    result_metrics: Mapping[str, Any],
    method: str,
) -> _AuditInput:
    counts = _raw_counts(data, params)
    genes = _gene_names(data)
    obs = data.obs.copy(deep=True).reset_index(drop=True)
    if len(obs) != counts.shape[0]:
        raise ValueError("observation metadata length does not match counts")
    condition_col = _resolve_column(obs, params.get("condition_col"), ("condition",), "condition")
    donor_col = _resolve_column(
        obs,
        params.get("donor_col"),
        ("donor", "donor_id", "mouse_id", "sample_id", "sample"),
        "donor",
    )
    if obs[[condition_col, donor_col]].isna().any().any():
        raise ValueError("donor/condition metadata contain missing values")
    condition_values = obs[condition_col].astype(str)
    donor_values = obs[donor_col].astype(str)
    if donor_values.str.lower().isin({"nan", "none", "<na>", "nat"}).any():
        raise ValueError("donor metadata contain missing values")
    observed = list(pd.unique(condition_values))
    condition_a = str(params.get("condition_a", result_metrics.get("condition_a", "")))
    condition_b = str(params.get("condition_b", result_metrics.get("condition_b", "")))
    if not condition_a or not condition_b or condition_a == "None" or condition_b == "None":
        if len(observed) != 2:
            raise ValueError("condition_a and condition_b are required for this contrast")
        condition_a, condition_b = map(str, observed)
    if condition_a == condition_b or condition_a not in observed or condition_b not in observed:
        raise ValueError("requested condition contrast is not observed")
    covariates = _covariates(obs, params)
    paired = bool(params.get("paired", params.get("paired_design", False)))
    donors_a = sorted(set(donor_values.loc[condition_values == condition_a]))
    donors_b = sorted(set(donor_values.loc[condition_values == condition_b]))
    overlap = set(donors_a) & set(donors_b)
    if paired and set(donors_a) != set(donors_b):
        raise ValueError("paired design has unmatched donor sets")
    if not paired and overlap:
        raise ValueError("independent design has donors shared across conditions")

    selected = obs.loc[condition_values.isin([condition_a, condition_b])].copy()
    selected["__condition"] = condition_values.loc[selected.index].to_numpy()
    selected["__donor"] = donor_values.loc[selected.index].to_numpy()
    bulk_rows: list[dict[str, Any]] = []
    for (donor, condition), positions in selected.groupby(["__donor", "__condition"], sort=True, observed=True).groups.items():
        positions = list(positions)
        record: dict[str, Any] = {
            "sample_id": f"{donor}__{condition}",
            "donor": str(donor),
            "condition": str(condition),
            donor_col: str(donor),
            condition_col: str(condition),
        }
        for covariate in covariates:
            values = selected.loc[positions, covariate]
            if values.isna().any():
                raise ValueError(f"covariate {covariate!r} contains missing values")
            unique = pd.unique(values)
            if len(unique) != 1:
                raise ValueError(
                    f"covariate {covariate!r} is not constant within donor-condition pseudobulk group"
                )
            record[covariate] = unique[0]
        bulk_rows.append(record)
    metadata = pd.DataFrame(bulk_rows).set_index("sample_id", drop=True)
    if metadata.index.duplicated().any() or metadata.empty:
        raise ValueError("donor-condition pseudobulk metadata are empty or duplicated")

    formula = _formula(params, paired, covariates)
    try:
        from formulaic import model_matrix

        design = pd.DataFrame(model_matrix(formula, metadata), index=metadata.index)
    except Exception as exc:
        raise ValueError(f"unable to materialize design {formula!r}: {exc}") from exc
    values = design.to_numpy(dtype=float)
    if values.ndim != 2 or not values.size or not np.isfinite(values).all():
        raise ValueError("design matrix is empty or non-finite")
    rank = int(np.linalg.matrix_rank(values))
    if rank < design.shape[1]:
        raise ValueError(f"design matrix is rank deficient (rank {rank} < {design.shape[1]})")
    residual_df = int(design.shape[0] - rank)
    min_loo = max(3, int(params.get("min_loo_donors", 3)))
    underpowered_loo = method == "pydeseq2_leave_one_donor_out_v1" and min(len(donors_a), len(donors_b)) < min_loo
    if residual_df <= 0 and not underpowered_loo:
        raise ValueError(f"design has no residual degrees of freedom ({design.shape[0]} samples, rank {rank})")
    return _AuditInput(
        counts=counts,
        genes=genes,
        obs=obs,
        metadata=metadata,
        donor_col=donor_col,
        condition_col=condition_col,
        condition_a=condition_a,
        condition_b=condition_b,
        paired=paired,
        covariates=covariates,
        donors_a=donors_a,
        donors_b=donors_b,
        formula=formula,
        design=design,
        design_rank=rank,
        residual_df=residual_df,
    )


def _canonical_network(frame: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    network = frame.copy(deep=True)
    if not {"source", "target"}.issubset(network.columns):
        raise ValueError("network must contain source and target columns")
    if "weight" not in network.columns:
        network["weight"] = 1.0
    if network[["source", "target"]].isna().any().any():
        raise ValueError("network source and target cannot be missing")
    network["source"] = network["source"].astype(str)
    network["target"] = network["target"].astype(str)
    network["weight"] = pd.to_numeric(network["weight"], errors="coerce")
    if network["weight"].isna().any() or not np.isfinite(network["weight"].to_numpy(dtype=float)).all():
        raise ValueError("network weights must be finite")
    network = network.loc[(network["source"] != "") & (network["target"] != "")].copy()
    if network.empty:
        raise ValueError("network has no non-empty interactions")
    canonical = network.sort_values(["source", "target", "weight"], kind="mergesort").to_csv(index=False).encode()
    return network, hashlib.sha256(canonical).hexdigest()


def _read_network_file(path_value: Any) -> tuple[pd.DataFrame, str, str]:
    path_text = str(path_value)
    if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", path_text):
        raise ValueError("network resources must be local; remote URLs are not accepted")
    path = Path(path_text).expanduser()
    if not path.is_file():
        raise ValueError(f"local network file does not exist: {path}")
    suffix = path.suffix.lower()
    if suffix == ".csv":
        frame = pd.read_csv(path)
    elif suffix in {".tsv", ".txt"}:
        frame = pd.read_csv(path, sep="\t")
    elif suffix in {".parquet", ".pq"}:
        frame = pd.read_parquet(path)
    elif suffix == ".json":
        frame = pd.read_json(path)
    else:
        raise ValueError(f"unsupported local network format {suffix!r}")
    return frame, str(path.resolve()), hashlib.sha256(path.read_bytes()).hexdigest()


def _resource_hashes(params: Mapping[str, Any]) -> list[str]:
    hashes: list[str] = []
    if params.get("network_sha256") is not None:
        hashes.append(str(params["network_sha256"]).lower())
    external = params.get("external_resource_sha256")
    if isinstance(external, Mapping):
        for key in ("network_path", "network", str(params.get("network_path", ""))):
            if key and key in external:
                hashes.append(str(external[key]).lower())
    return hashes


def _network_for_audit(contract, registry) -> tuple[pd.DataFrame, dict[str, Any]]:
    params = contract.parameters
    network = params.get("network")
    source_kind = "parameters.network"
    resource_path: Optional[str] = None
    raw_hash: Optional[str] = None
    network_source = network_version = None
    if network is None and params.get("network_path") is not None:
        network, resource_path, raw_hash = _read_network_file(params["network_path"])
        source_kind = "local_file"
    if network is None and len(contract.input_artifacts) > 1:
        metadata, payload = registry.get(contract.input_artifacts[1])
        if isinstance(payload, pd.DataFrame):
            network = payload
        elif isinstance(payload, Mapping):
            network = pd.DataFrame(payload)
        if network is not None:
            source_kind = "table_artifact"
            network_source = metadata.summary_metrics.get("network_source")
            network_version = metadata.summary_metrics.get("network_version")
    if isinstance(network, (str, Path)):
        network, resource_path, raw_hash = _read_network_file(network)
        source_kind = "local_file"
    elif isinstance(network, Mapping):
        network = pd.DataFrame(network)
    if not isinstance(network, pd.DataFrame):
        raise ValueError("a local network DataFrame, path, or table artifact is required")
    network, canonical_hash = _canonical_network(network)
    digest = raw_hash or canonical_hash
    expected = _resource_hashes(params)
    if expected and any(value != digest.lower() for value in expected):
        raise ValueError("network hash mismatch")
    species = params.get("species")
    network_source = params.get("network_source", network_source)
    network_version = params.get("network_version", network_version)
    if species is None or not str(species).strip():
        raise ValueError("species is required")
    if network_source is None or not str(network_source).strip():
        raise ValueError("network_source is required")
    if network_version is None or not str(network_version).strip():
        raise ValueError("network_version is required")
    if "species" in network.columns:
        species_values = {str(value) for value in network["species"].dropna().unique()}
        if species_values and species_values != {str(species)}:
            raise ValueError("network species does not match requested species")
    return network, {
        "species": str(species),
        "network_source": str(network_source),
        "network_version": str(network_version),
        "network_sha256": digest,
        "external_resource_sha256": digest,
        "network_resource_kind": source_kind,
        "network_path": resource_path,
    }


def _finite_ci_errors(table: pd.DataFrame, alpha: float) -> list[str]:
    errors: list[str] = []
    if {"ci_low", "ci_high", "lfc_se", "log2_fold_change"}.issubset(table.columns):
        z = float(norm.ppf(1.0 - alpha / 2.0))
        frame = table[["ci_low", "ci_high", "lfc_se", "log2_fold_change"]].apply(pd.to_numeric, errors="coerce").dropna()
        if (frame["lfc_se"] < 0).any():
            errors.append("confidence interval standard errors must be non-negative")
        if not np.allclose(frame["ci_low"], frame["log2_fold_change"] - z * frame["lfc_se"], rtol=1e-6, atol=1e-8):
            errors.append("invalid lower confidence interval")
        if not np.allclose(frame["ci_high"], frame["log2_fold_change"] + z * frame["lfc_se"], rtol=1e-6, atol=1e-8):
            errors.append("invalid upper confidence interval")
    if {"activity_ci_low", "activity_ci_high", "activity_se", "activity_effect_condition_a_vs_b"}.issubset(table.columns):
        activity_columns = ["activity_ci_low", "activity_ci_high", "activity_se", "activity_effect_condition_a_vs_b"]
        frame = table[activity_columns].apply(pd.to_numeric, errors="coerce")
        has_residual_df = "residual_degrees_of_freedom" in table
        if not has_residual_df:
            errors.append("activity confidence intervals require residual degrees of freedom")
            residual_df = pd.Series(np.nan, index=table.index, dtype=float)
        else:
            residual_df = pd.to_numeric(table["residual_degrees_of_freedom"], errors="coerce")
        residual_values = residual_df.to_numpy(dtype=float)
        valid_residual_df = (
            has_residual_df
            and (
                len(residual_values) == 0
                or (np.isfinite(residual_values).all() and (residual_values > 0).all())
            )
        )
        if has_residual_df and not valid_residual_df:
            errors.append("activity confidence intervals require finite positive residual degrees of freedom")
        if (frame["activity_se"].dropna() < 0).any():
            errors.append("activity standard errors must be non-negative")
        finite = frame.notna().all(axis=1) & residual_df.notna()
        if valid_residual_df and finite.any():
            critical = t.ppf(1.0 - alpha / 2.0, residual_df.loc[finite].to_numpy(dtype=float))
            expected_low = frame.loc[finite, "activity_effect_condition_a_vs_b"] - critical * frame.loc[finite, "activity_se"]
            expected_high = frame.loc[finite, "activity_effect_condition_a_vs_b"] + critical * frame.loc[finite, "activity_se"]
            if not np.allclose(frame.loc[finite, "activity_ci_low"], expected_low, rtol=1e-6, atol=1e-8):
                errors.append("invalid activity lower confidence interval")
            if not np.allclose(frame.loc[finite, "activity_ci_high"], expected_high, rtol=1e-6, atol=1e-8):
                errors.append("invalid activity upper confidence interval")
    return errors


def _table_statistics_errors(table: pd.DataFrame, params: Mapping[str, Any]) -> list[str]:
    errors: list[str] = []
    if not isinstance(table, pd.DataFrame) or not len(table.columns):
        return ["result must be a non-empty structured table"]
    numeric = table.select_dtypes(include=[np.number])
    if np.isinf(numeric.to_numpy(dtype=float)).any():
        errors.append("result contains infinite numeric values")
    for column in ("p_value", "fdr_q_value", "padj", "decoupler_p_value"):
        if column in table:
            values = pd.to_numeric(table[column], errors="coerce").dropna()
            if not values.between(0, 1).all():
                errors.append(f"{column} is outside [0, 1]")
    if "fdr_q_value" in table and "p_value" in table:
        q = pd.to_numeric(table["fdr_q_value"], errors="coerce")
        p = pd.to_numeric(table["p_value"], errors="coerce")
        if (q.notna() & p.isna()).any():
            errors.append("adjusted significance has no raw p-value")
    if "status" in table:
        status = table["status"].astype(str)
        reason = table.get("status_reason", table.get("skip_reason", pd.Series("", index=table.index))).fillna("").astype(str)
        explainable = status.isin({"not_estimated", "estimated_fdr_missing", "skipped"})
        if (explainable & reason.str.strip().eq("")).any():
            errors.append("unestimated or skipped rows require an explicit reason")
        if "p_value" in table:
            p = pd.to_numeric(table["p_value"], errors="coerce")
            if (p.isna() & ~explainable).any():
                errors.append("missing p-values are not marked not_estimated")
        if "fdr_q_value" in table:
            q = pd.to_numeric(table["fdr_q_value"], errors="coerce")
            if (q.isna() & ~explainable).any():
                errors.append("missing FDR values are not marked not_estimated or estimated_fdr_missing")
    errors.extend(_finite_ci_errors(table, float(params.get("alpha", 0.05))))
    return errors


def _check_summary_against_fits(summary: pd.DataFrame, primary: pd.DataFrame, details: _AuditInput) -> list[str]:
    errors: list[str] = []
    required = {
        "gene", "baseline_log2_fold_change", "baseline_fdr_q_value", "n_estimated_loo",
        "estimated_coverage", "direction_consistency", "min_log2_fold_change",
        "max_log2_fold_change", "median_log2_fold_change", "significance_retention",
    }
    if not required.issubset(summary.columns):
        return [f"LOO summary is missing columns {sorted(required - set(summary.columns))}"]
    genes = set(details.genes)
    if set(summary["gene"].astype(str)) != genes or summary["gene"].astype(str).duplicated().any():
        errors.append("LOO summary does not contain exactly one row per input gene")
        return errors
    primary = primary.copy()
    if "left_out_donor" not in primary:
        return errors + ["LOO primary output lacks left_out_donor"]
    alpha = .05
    for gene in details.genes:
        base = primary.loc[(primary["gene"].astype(str) == gene) & (primary["left_out_donor"].astype(str) == "__full_model__")]
        fits = primary.loc[(primary["gene"].astype(str) == gene) & (primary["left_out_donor"].astype(str) != "__full_model__")]
        row = summary.loc[summary["gene"].astype(str) == gene].iloc[0]
        effects = pd.to_numeric(fits.get("log2_fold_change", pd.Series(dtype=float)), errors="coerce").dropna()
        q_values = pd.to_numeric(fits.get("fdr_q_value", pd.Series(dtype=float)), errors="coerce").dropna()
        baseline_lfc = pd.to_numeric(base.get("log2_fold_change", pd.Series(dtype=float)), errors="coerce").dropna()
        baseline_q = pd.to_numeric(base.get("fdr_q_value", pd.Series(dtype=float)), errors="coerce").dropna()
        if len(base) != 1:
            errors.append(f"gene {gene} lacks one full-model baseline row")
            continue
        expected_n = len(effects)
        if int(row["n_estimated_loo"]) != expected_n:
            errors.append(f"gene {gene} has incorrect n_estimated_loo")
        requested = int(pd.to_numeric(primary["n_requested_fits"], errors="coerce").dropna().iloc[0]) if "n_requested_fits" in primary and primary["n_requested_fits"].notna().any() else 0
        expected_coverage = expected_n / requested if requested else np.nan
        actual_coverage = pd.to_numeric(pd.Series([row["estimated_coverage"]]), errors="coerce").iloc[0]
        if np.isfinite(expected_coverage) and (not np.isfinite(actual_coverage) or not np.isclose(actual_coverage, expected_coverage)):
            errors.append(f"gene {gene} has incorrect estimated coverage")
        if len(effects) and len(baseline_lfc) and float(baseline_lfc.iloc[0]) != 0:
            expected_direction = float((np.sign(effects.to_numpy()) == np.sign(float(baseline_lfc.iloc[0]))).mean())
            actual_direction = pd.to_numeric(pd.Series([row["direction_consistency"]]), errors="coerce").iloc[0]
            if not np.isfinite(actual_direction) or not np.isclose(actual_direction, expected_direction):
                errors.append(f"gene {gene} has incorrect direction consistency")
        if len(effects):
            for column, expected_value in (("min_log2_fold_change", effects.min()), ("max_log2_fold_change", effects.max()), ("median_log2_fold_change", effects.median())):
                actual = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
                if not np.isfinite(actual) or not np.isclose(actual, expected_value):
                    errors.append(f"gene {gene} has incorrect {column}")
        if len(q_values) and len(baseline_q):
            expected_retention = float(((q_values.to_numpy() < alpha) == (float(baseline_q.iloc[0]) < alpha)).mean())
            actual_retention = pd.to_numeric(pd.Series([row["significance_retention"]]), errors="coerce").iloc[0]
            if not np.isfinite(actual_retention) or not np.isclose(actual_retention, expected_retention):
                errors.append(f"gene {gene} has incorrect significance retention")
    return errors


class AdvancedStatisticsValidator(BaseAuditor):
    """Independently audit input scale, design, provenance, and statistics."""

    def __init__(self):
        super().__init__("advanced_statistics_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(auditor_name=self.auditor_name, target_task_id=contract.task_id)
        if result.method_used not in METHODS:
            return report
        errors: list[str] = []
        method = result.method_used
        params = dict(contract.parameters)
        details: Optional[_AuditInput] = None
        tables: list[pd.DataFrame] = []

        try:
            _, payload = registry.get(contract.input_artifacts[0])
            data = payload if isinstance(payload, SCData) else _coerce_data(payload)
            details = _reconstruct_input(data, params, result.metrics, method)
            if method != "pydeseq2_leave_one_donor_out_v1":
                threshold = max(2, int(params.get("min_donors", 2)))
                if min(len(details.donors_a), len(details.donors_b)) < threshold:
                    errors.append(f"fewer than {threshold} donors per condition")
        except Exception as exc:
            errors.append(f"input reconstruction failed: {type(exc).__name__}: {exc}")

        try:
            if contract.expected_outputs and list(result.output_artifacts) != list(contract.expected_outputs):
                errors.append("result output artifacts do not satisfy expected_outputs")
            for uri in result.output_artifacts:
                _, payload = registry.get(uri)
                tables.append(payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload))
        except Exception as exc:
            errors.append(f"result table retrieval failed: {type(exc).__name__}: {exc}")

        input_errors = [e for e in errors if e.startswith("input reconstruction") or "fewer than" in e]
        report.add_check(
            "advanced_statistics_input_counts_metadata",
            not input_errors,
            ValidationSeverity.ERROR,
            "; ".join(input_errors) if input_errors else "Raw integer counts, metadata, donor pairing, and covariate constancy were independently reconstructed.",
        )

        design_errors: list[str] = []
        if details is None:
            design_errors.append("design cannot be audited because input reconstruction failed")
        else:
            if result.metrics.get("design_rank") is not None and int(result.metrics["design_rank"]) != details.design_rank:
                design_errors.append("reported design rank does not match the reconstructed design")
            if result.metrics.get("design") is not None and str(result.metrics["design"]) != details.formula:
                design_errors.append("reported design formula does not match parameters and metadata")
            for table in tables:
                if "design" in table:
                    values = {str(v) for v in table["design"].dropna().unique()}
                    if values and values != {details.formula}:
                        design_errors.append("result table design formula differs from reconstructed design")
                if "design_rank" in table:
                    ranks = pd.to_numeric(table["design_rank"], errors="coerce").dropna()
                    if len(ranks) and not np.all(ranks.astype(int) == details.design_rank):
                        design_errors.append("result table design rank differs from reconstructed design")
                if "residual_degrees_of_freedom" in table:
                    residual = pd.to_numeric(table["residual_degrees_of_freedom"], errors="coerce").dropna()
                    if len(residual) and not np.all(residual.astype(int) == details.residual_df):
                        design_errors.append("result table residual degrees of freedom differs from reconstructed design")
                if "covariates" in table:
                    expected_covariates = ",".join(details.covariates)
                    values = {str(v) for v in table["covariates"].dropna().unique()}
                    if values and values != {expected_covariates}:
                        design_errors.append("result table covariates differ from the requested design")
        report.add_check(
            "advanced_statistics_design_rank",
            not design_errors,
            ValidationSeverity.ERROR,
            "; ".join(dict.fromkeys(design_errors)) if design_errors else "Design formula, covariates, rank, and residual degrees of freedom are consistent.",
            {"design_rank": details.design_rank, "residual_degrees_of_freedom": details.residual_df} if details else {},
        )
        errors.extend(design_errors)

        resource_errors: list[str] = []
        if method == "decoupler_ulm_v2":
            try:
                network, provenance = _network_for_audit(contract, registry)
                if details is not None:
                    overlap = network.loc[network["target"].isin(details.genes)]
                    tmin = max(1, int(params.get("tmin", params.get("min_network_size", 1))))
                    sources = overlap.groupby("source", observed=True)["target"].nunique()
                    if overlap.empty or not (sources >= tmin).any():
                        resource_errors.append("local network has no source meeting the target overlap threshold")
                for table in tables:
                    for key, expected in provenance.items():
                        if key not in {"network_path", "network_resource_kind"} and key in table:
                            values = {str(v) for v in table[key].dropna().unique()}
                            if values and values != {str(expected)}:
                                resource_errors.append(f"result network provenance {key} differs from local resource")
                        elif key in {"species", "network_source", "network_version", "network_sha256", "external_resource_sha256"}:
                            resource_errors.append(f"result table lacks network provenance {key}")
            except Exception as exc:
                resource_errors.append(f"network provenance audit failed: {type(exc).__name__}: {exc}")
        report.add_check(
            "advanced_statistics_resource_provenance",
            not resource_errors,
            ValidationSeverity.ERROR,
            "; ".join(dict.fromkeys(resource_errors)) if resource_errors else ("Local network species, source, version, and hash are independently verified." if method == "decoupler_ulm_v2" else "No external network resource is used by this method."),
        )
        errors.extend(resource_errors)

        result_errors: list[str] = []
        for table in tables:
            result_errors.extend(_table_statistics_errors(table, params))
        if method == "pydeseq2_pseudobulk_v1" and details is not None:
            gene_tables = [table for table in tables if "gene" in table]
            if not gene_tables:
                result_errors.append("DEG output lacks a gene column")
            else:
                genes = gene_tables[0]["gene"].astype(str)
                if set(genes) != set(details.genes) or genes.duplicated().any():
                    result_errors.append("DEG output does not preserve one row for every input gene")
        if method == "decoupler_ulm_v2":
            for table in tables:
                if "source" not in table or table["source"].duplicated().any():
                    result_errors.append("functional activity output must contain one row per network source")
                if "fdr_q_value" not in table or "p_value" not in table:
                    result_errors.append("functional activity output lacks raw p-values or FDR")
        report.add_check(
            "advanced_statistics_result_statistics",
            not result_errors,
            ValidationSeverity.ERROR,
            "; ".join(dict.fromkeys(result_errors)) if result_errors else "Result status, probability ranges, intervals, and retained features are consistent.",
        )
        errors.extend(result_errors)

        loo_errors: list[str] = []
        if method == "pydeseq2_leave_one_donor_out_v1":
            if details is None:
                loo_errors.append("LOO coverage cannot be audited because input reconstruction failed")
            elif not tables:
                loo_errors.append("LOO produced no output table")
            else:
                primary = tables[0]
                threshold = max(3, int(params.get("min_loo_donors", 3)))
                n_a, n_b = len(details.donors_a), len(details.donors_b)
                skipped = bool(result.metrics.get("skipped")) or ("status" in primary and primary["status"].astype(str).eq("skipped").all())
                if skipped:
                    if n_a >= threshold and n_b >= threshold:
                        loo_errors.append("LOO was skipped despite meeting the donor threshold")
                    if not bool(result.metrics.get("skipped")):
                        loo_errors.append("skipped LOO table is not marked skipped in result metrics")
                    for table in tables:
                        if "status" not in table or not table["status"].astype(str).eq("skipped").all():
                            loo_errors.append("every expected LOO output must carry skipped status")
                        if "skip_reason" not in table or table["skip_reason"].fillna("").astype(str).str.strip().eq("").any():
                            loo_errors.append("skipped LOO outputs require a reason")
                    if float(result.metrics.get("complete_fit_coverage", np.nan)) != 0.0:
                        loo_errors.append("skipped LOO must report zero complete-fit coverage")
                else:
                    if n_a < threshold or n_b < threshold:
                        loo_errors.append("LOO ran below its minimum donor threshold")
                    if bool(result.metrics.get("skipped")):
                        loo_errors.append("LOO result metrics claim skipped despite fitted rows")
                    if "left_out_donor" not in primary or "gene" not in primary:
                        loo_errors.append("LOO primary output lacks donor or gene identifiers")
                    else:
                        genes = set(details.genes)
                        baseline = primary.loc[primary["left_out_donor"].astype(str) == "__full_model__"]
                        if set(baseline["gene"].astype(str)) != genes or baseline["gene"].astype(str).duplicated().any():
                            loo_errors.append("LOO full-model baseline does not preserve the full gene universe")
                        donor_ids = set(details.donors_a) | set(details.donors_b)
                        observed = set(primary.loc[primary["left_out_donor"].astype(str) != "__full_model__", "left_out_donor"].astype(str))
                        if not observed.issubset(donor_ids):
                            loo_errors.append("LOO output contains an unexpected donor")
                        missing_donors = donor_ids - observed
                        n_successful = sum(
                            len(primary.loc[primary["left_out_donor"].astype(str) == donor, "gene"].astype(str).unique()) == len(genes)
                            for donor in observed
                        )
                        claimed_successful = int(result.metrics.get("n_successful_fits", -1))
                        if claimed_successful != n_successful:
                            loo_errors.append("reported successful LOO fits do not match complete donor rows")
                        if missing_donors and int(result.metrics.get("n_failed_fits", -1)) != len(missing_donors):
                            loo_errors.append("missing donor LOO fits are not recorded as failures")
                        requested = len(donor_ids)
                        expected_coverage = n_successful / requested if requested else 0.0
                        coverage = float(result.metrics.get("complete_fit_coverage", np.nan))
                        if not np.isclose(coverage, expected_coverage):
                            loo_errors.append("complete-fit coverage does not match donor-level fits")
                        for donor in observed:
                            donor_rows = primary.loc[primary["left_out_donor"].astype(str) == donor]
                            if set(donor_rows["gene"].astype(str)) != genes or donor_rows["gene"].astype(str).duplicated().any():
                                loo_errors.append(f"LOO donor {donor} does not have a complete gene table")
                            expected_removed = 2 if details.paired else 1
                            if "removed_sample_count" not in donor_rows or set(pd.to_numeric(donor_rows["removed_sample_count"], errors="coerce").dropna().astype(int)) != {expected_removed}:
                                loo_errors.append(f"LOO donor {donor} does not remove the expected complete donor sample")
                        if "paired_leaveout_removes_complete_donor" in primary and details.paired and not primary["paired_leaveout_removes_complete_donor"].astype(bool).all():
                            loo_errors.append("paired LOO does not declare complete-donor removal")
                        summary = next((table for table in tables[1:] if {"direction_consistency", "estimated_coverage"}.issubset(table.columns)), None)
                        if summary is None and {"direction_consistency", "estimated_coverage"}.issubset(primary.columns):
                            summary = primary.drop_duplicates("gene")
                        if len(contract.expected_outputs) >= 2 and summary is None:
                            loo_errors.append("expected LOO summary output is absent")
                        if summary is not None:
                            loo_errors.extend(_check_summary_against_fits(summary, primary, details))
                    if result.metrics.get("scientific_robustness_claim_supported") is not False:
                        loo_errors.append("LOO must not claim scientific robustness from fit completion")
                    for table in tables:
                        if "scientific_robustness_claim_supported" in table and table["scientific_robustness_claim_supported"].astype(bool).any():
                            loo_errors.append("LOO output contains an unsupported robustness claim")
        report.add_check(
            "advanced_statistics_loo_coverage",
            not loo_errors,
            ValidationSeverity.ERROR,
            "; ".join(dict.fromkeys(loo_errors)) if loo_errors else ("Every requested donor LOO fit removes the complete donor and has independently checked effect sensitivity." if method == "pydeseq2_leave_one_donor_out_v1" else "Not applicable to this method."),
        )
        errors.extend(loo_errors)

        report.add_check(
            "advanced_statistics_integrity",
            not errors,
            ValidationSeverity.ERROR,
            "; ".join(dict.fromkeys(errors)) if errors else "Independently checked raw counts, donor design, resource provenance, result status, and sensitivity coverage.",
        )
        return report
