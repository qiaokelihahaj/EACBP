"""Advanced donor-level statistical capabilities.

This module contains the optional, library-backed statistics used by the
standard analysis profile.  It deliberately keeps the optional imports lazy:
the baseline EACBP installation remains importable when ``pydeseq2`` and
``decoupler`` are not installed.

The three public capabilities are:

``PyDESeq2PseudobulkCapability``
    Aggregates integer raw counts by donor and condition and fits a real
    PyDESeq2 negative-binomial model.  The model supports paired donor
    blocking, covariates, an explicit formula and an explicit contrast.

``DecouplerFunctionalAnalysisCapability``
    Runs the current decoupler ULM API on donor-level log-CPM expression and
    compares inferred pathway or TF activities at the donor level.  A prior
    network is always supplied by the caller from a local path, table
    artifact, or in-memory DataFrame.  No decoupler resource downloader is
    called here.

``PyDESeq2LeaveOneDonorOutCapability``
    Re-fits the PyDESeq2 model after removing one complete biological donor at
    a time.  A paired leave-one-out replicate removes both condition samples
    belonging to that donor.  Fewer than three donors per condition (or pairs)
    yields an explicit skip artifact instead of a robustness claim.

The implementation is intentionally independent of the legacy Welch DEG
capability.  Shared helpers are kept private to avoid changing the behavior
of existing workflows.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.numerics import benjamini_hochberg
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


_LN2 = float(np.log(2.0))
_DEFAULT_MIN_DONORS = 2
_DEFAULT_LOO_MIN_DONORS = 3


class AdvancedStatisticsInputError(ValueError):
    """Raised when a statistical prerequisite is not satisfied."""


class AdvancedStatisticsDependencyError(ImportError):
    """Raised when an optional scientific dependency is unavailable."""


@dataclass
class _PseudobulkInput:
    """Prepared donor-condition counts and design information."""

    counts: pd.DataFrame
    metadata: pd.DataFrame
    condition_a: str
    condition_b: str
    condition_col: str
    donor_col: str
    counts_source: str
    design: str
    paired: bool
    covariates: list[str]
    donor_ids: list[str]
    donor_ids_a: list[str]
    donor_ids_b: list[str]
    design_rank: Optional[int] = None
    design_columns: tuple[str, ...] = ()


def _software_version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unavailable"


def _require_pydeseq2():
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:  # pragma: no cover - exercised without extra
        raise AdvancedStatisticsDependencyError(
            "PyDESeq2 is required for this capability; install eacbp[advanced-statistics]"
        ) from exc
    return DeseqDataSet, DeseqStats


def _require_decoupler():
    try:
        import decoupler as dc
    except ImportError as exc:  # pragma: no cover - exercised without extra
        raise AdvancedStatisticsDependencyError(
            "decoupler is required for this capability; install eacbp[advanced-statistics]"
        ) from exc
    return dc


def _as_dense(value: Any) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    elif hasattr(value, "A"):
        value = value.A
    return np.asarray(value)


def _coerce_data(payload: Any) -> SCData:
    if isinstance(payload, SCData):
        return payload.copy()
    if isinstance(payload, Mapping) and "X" in payload and "obs" in payload:
        return SCData.from_dict(dict(payload))
    if hasattr(payload, "X") and hasattr(payload, "obs") and hasattr(payload, "var"):
        return SCData.from_anndata(payload)
    raise TypeError(
        "advanced statistical capabilities require an SCData, AnnData-like payload, "
        "or SCData dictionary"
    )


def _resolve_column(obs: pd.DataFrame, explicit: Any, candidates: Sequence[str], label: str) -> str:
    if explicit is not None:
        name = str(explicit)
        if name not in obs.columns:
            raise AdvancedStatisticsInputError(
                f"{label} metadata column {name!r} is not present; available columns: {list(obs.columns)}"
            )
        return name
    for name in candidates:
        if name in obs.columns:
            return name
    raise AdvancedStatisticsInputError(
        f"{label} metadata is required; provide {label.lower()}_col explicitly"
    )


def _resolve_conditions(obs: pd.DataFrame, params: Mapping[str, Any], condition_col: str) -> tuple[str, str]:
    values = obs[condition_col]
    if values.isna().any():
        raise AdvancedStatisticsInputError(
            f"condition metadata column {condition_col!r} contains missing values"
        )
    observed = [str(v) for v in pd.unique(values)]
    if len(observed) < 2:
        raise AdvancedStatisticsInputError(
            f"at least two observed conditions are required; found {observed}"
        )
    requested_a = params.get("condition_a")
    requested_b = params.get("condition_b")
    if requested_a is None and requested_b is None:
        if len(observed) != 2:
            raise AdvancedStatisticsInputError(
                "condition_a and condition_b are required when more than two conditions are observed"
            )
        return observed[0], observed[1]
    if requested_a is None or requested_b is None:
        raise AdvancedStatisticsInputError("condition_a and condition_b must be supplied together")
    cond_a, cond_b = str(requested_a), str(requested_b)
    if cond_a == cond_b:
        raise AdvancedStatisticsInputError("condition_a and condition_b must be different")
    if cond_a not in observed or cond_b not in observed:
        raise AdvancedStatisticsInputError(
            f"requested conditions {cond_a!r}, {cond_b!r} are not both observed: {observed}"
        )
    return cond_a, cond_b


def _resolve_covariates(obs: pd.DataFrame, params: Mapping[str, Any]) -> list[str]:
    raw = params.get("covariates", params.get("covariate_cols", []))
    if raw is None:
        return []
    if isinstance(raw, str):
        raw = [raw]
    if not isinstance(raw, Iterable):
        raise AdvancedStatisticsInputError("covariates must be a string or a sequence of metadata columns")
    result: list[str] = []
    for value in raw:
        name = str(value)
        if name in result:
            continue
        if name not in obs.columns:
            raise AdvancedStatisticsInputError(f"covariate metadata column {name!r} is not present")
        result.append(name)
    return result


def _validate_counts(data: SCData, params: Mapping[str, Any]) -> tuple[Any, str]:
    layers = getattr(data, "layers", {}) or {}
    layer_name = params.get("counts_layer", params.get("raw_counts_layer", "counts"))
    counts_source = f"layers.{layer_name}"
    if layer_name not in layers:
        if bool(params.get("allow_x_as_counts", False)):
            # This opt-in is useful for AnnData files whose X is documented as
            # raw counts.  It is still labelled explicitly in provenance.
            value = data.X
            counts_source = "X_explicitly_declared_raw_counts"
        else:
            raise AdvancedStatisticsInputError(
                f"raw integer counts are required in layers[{layer_name!r}]; "
                "set allow_x_as_counts=True only when X is documented as raw counts"
            )
    else:
        value = layers[layer_name]
    arr = _as_dense(value)
    if arr.ndim != 2 or arr.shape != data.shape:
        raise AdvancedStatisticsInputError(
            f"counts matrix must have shape {data.shape}; received {getattr(arr, 'shape', None)}"
        )
    if arr.dtype.kind in "fc" and not np.isfinite(arr).all():
        raise AdvancedStatisticsInputError("raw counts contain NaN or infinite values")
    if arr.dtype.kind not in "biufc":
        raise AdvancedStatisticsInputError(f"raw counts have unsupported dtype {arr.dtype}")
    if np.any(arr < 0):
        raise AdvancedStatisticsInputError("raw counts must be non-negative")
    if not np.all(np.equal(arr, np.floor(arr))):
        raise AdvancedStatisticsInputError(
            "raw counts must be integer-valued; normalized/log-transformed values are not accepted"
        )
    max_int = np.iinfo(np.int64).max
    if np.any(arr > max_int):
        raise AdvancedStatisticsInputError("raw counts exceed int64 range")
    return np.asarray(arr, dtype=np.int64), counts_source


def _gene_names(data: SCData) -> list[str]:
    if "gene_name" in data.var.columns:
        genes = data.var["gene_name"].astype(str).tolist()
    else:
        genes = [str(v) for v in data.var.index.tolist()]
    if len(genes) != data.n_vars:
        raise AdvancedStatisticsInputError("gene metadata length does not match the count matrix")
    if any(not gene or gene.lower() in {"nan", "none"} for gene in genes):
        raise AdvancedStatisticsInputError("gene names must be non-empty")
    if len(set(genes)) != len(genes):
        duplicates = sorted({gene for gene in genes if genes.count(gene) > 1})[:10]
        raise AdvancedStatisticsInputError(f"gene names must be unique; duplicates include {duplicates}")
    return genes


def _aggregate_pseudobulk(
    counts: np.ndarray,
    obs: pd.DataFrame,
    genes: Sequence[str],
    condition_col: str,
    donor_col: str,
    condition_a: str,
    condition_b: str,
    covariates: Sequence[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    condition_values = obs[condition_col].astype(str)
    raw_donors = obs[donor_col]
    if raw_donors.isna().any():
        raise AdvancedStatisticsInputError(f"donor metadata column {donor_col!r} contains missing values")
    donor_values = raw_donors.astype(str)
    mask = condition_values.isin([condition_a, condition_b]).to_numpy()
    if not mask.any():
        raise AdvancedStatisticsInputError("no cells belong to the requested condition contrast")
    # Reset the selected observation index because groupby(...).groups returns
    # labels from the frame index, while ``selected_counts`` is addressed by
    # positional rows.  Cell IDs are often strings and must never be passed as
    # NumPy positions accidentally.
    selected_obs = obs.loc[mask].copy().reset_index(drop=True)
    selected_obs["__condition"] = condition_values.loc[mask].to_numpy()
    selected_obs["__donor"] = donor_values.loc[mask].to_numpy()
    if selected_obs["__donor"].isna().any() or selected_obs["__donor"].str.lower().isin({"nan", "none", "<na>", "nat"}).any():
        raise AdvancedStatisticsInputError(f"donor metadata column {donor_col!r} contains missing values")
    selected_counts = counts[mask]

    rows: list[np.ndarray] = []
    metadata_rows: list[dict[str, Any]] = []
    # groupby(sort=True) gives deterministic pseudobulk order independent of
    # cell order while retaining the original counts exactly.
    for (donor, condition), positions in selected_obs.groupby(
        ["__donor", "__condition"], sort=True, observed=True
    ).groups.items():
        pos = np.asarray(list(positions), dtype=int)
        rows.append(selected_counts[pos].sum(axis=0, dtype=np.int64))
        record: dict[str, Any] = {"donor": str(donor), "condition": str(condition)}
        # Store the caller's names too, so a user formula written with custom
        # metadata column names remains inspectable in the metadata artifact.
        record[donor_col] = str(donor)
        record[condition_col] = str(condition)
        for covariate in covariates:
            values = selected_obs.iloc[pos][covariate]
            if values.isna().any():
                raise AdvancedStatisticsInputError(
                    f"covariate {covariate!r} contains missing values in donor-condition group "
                    f"({donor!r}, {condition!r})"
                )
            unique = pd.unique(values)
            if len(unique) != 1:
                raise AdvancedStatisticsInputError(
                    f"covariate {covariate!r} must be constant within each donor-condition group; "
                    f"group ({donor!r}, {condition!r}) has {len(unique)} values"
                )
            record[covariate] = unique[0]
        record["sample_id"] = f"{donor}__{condition}"
        metadata_rows.append(record)
    if not rows:
        raise AdvancedStatisticsInputError("the requested contrast produced no donor-condition samples")
    metadata = pd.DataFrame(metadata_rows)
    if metadata["sample_id"].duplicated().any():
        raise AdvancedStatisticsInputError("donor-condition pseudobulk sample identifiers are not unique")
    metadata = metadata.set_index("sample_id", drop=True)
    counts_df = pd.DataFrame(np.vstack(rows), index=metadata.index, columns=list(genes), dtype=np.int64)
    return counts_df, metadata


def _formula_variables(formula: str) -> set[str]:
    # This is only used to determine whether an automatically requested donor
    # block/covariate is already present.  PyDESeq2/formulaic remains the
    # authoritative parser and rank checker.
    cleaned = re.sub(r"\bC\s*\(\s*([^)]*)\s*\)", r"\1", formula)
    return set(re.findall(r"[A-Za-z_]\w*", cleaned))


def _make_design(
    metadata: pd.DataFrame,
    params: Mapping[str, Any],
    paired: bool,
    covariates: Sequence[str],
) -> str:
    explicit = params.get("design_formula", params.get("design"))
    if explicit is not None:
        design = str(explicit).strip()
        if not design:
            raise AdvancedStatisticsInputError("design_formula cannot be empty")
        if not design.startswith("~"):
            design = "~" + design
        variables = _formula_variables(design)
        terms: list[str] = []
        if paired and "donor" not in variables:
            terms.append("donor")
        for covariate in covariates:
            if covariate not in variables:
                terms.append(covariate)
        if terms:
            design += " + " + " + ".join(terms)
        return design
    terms = ["condition"]
    if paired:
        terms.append("donor")
    terms.extend(covariates)
    return "~" + " + ".join(terms)


def _validate_sample_qualification(
    metadata: pd.DataFrame,
    condition_a: str,
    condition_b: str,
    paired: bool,
    min_donors: int,
) -> tuple[list[str], list[str], list[str]]:
    by_condition = {
        condition: set(metadata.loc[metadata["condition"] == condition, "donor"].astype(str))
        for condition in (condition_a, condition_b)
    }
    donors_a = sorted(by_condition[condition_a])
    donors_b = sorted(by_condition[condition_b])
    overlap = sorted(set(donors_a).intersection(donors_b))
    if paired:
        if len(overlap) < min_donors:
            raise AdvancedStatisticsInputError(
                f"paired design requires at least {min_donors} donors observed in both conditions; "
                f"found {len(overlap)} ({overlap})"
            )
        # A paired model must not silently discard unmatched donors.  If a
        # caller explicitly asks for paired analysis, unmatched donor samples
        # would change the intended pairing and are rejected.
        if set(donors_a) != set(donors_b):
            raise AdvancedStatisticsInputError(
                "paired=True requires every donor to have both conditions; "
                f"condition A donors={donors_a}, condition B donors={donors_b}"
            )
        return donors_a, donors_b, overlap
    if overlap:
        raise AdvancedStatisticsInputError(
            "the same donor occurs in both conditions; set paired=True for a paired design"
        )
    if len(donors_a) < min_donors or len(donors_b) < min_donors:
        raise AdvancedStatisticsInputError(
            f"independent design requires at least {min_donors} donors per condition; "
            f"found {len(donors_a)} and {len(donors_b)}"
        )
    return donors_a, donors_b, []


def _prepare_pseudobulk(data: SCData, params: Mapping[str, Any]) -> _PseudobulkInput:
    obs = data.obs.copy(deep=True)
    condition_col = _resolve_column(obs, params.get("condition_col"), ("condition",), "Condition")
    donor_col = _resolve_column(
        obs,
        params.get("donor_col"),
        ("donor", "donor_id", "mouse_id", "sample_id", "sample"),
        "Donor",
    )
    condition_a, condition_b = _resolve_conditions(obs, params, condition_col)
    covariates = _resolve_covariates(obs, params)
    paired = bool(params.get("paired", params.get("paired_design", False)))
    counts, counts_source = _validate_counts(data, params)
    genes = _gene_names(data)
    counts_df, metadata = _aggregate_pseudobulk(
        counts,
        obs,
        genes,
        condition_col,
        donor_col,
        condition_a,
        condition_b,
        covariates,
    )
    donors_a, donors_b, paired_donors = _validate_sample_qualification(
        metadata,
        condition_a,
        condition_b,
        paired,
        max(2, int(params.get("min_donors", _DEFAULT_MIN_DONORS))),
    )
    design = _make_design(metadata, params, paired, covariates)
    # Normalize custom metadata columns into canonical columns while retaining
    # aliases.  This gives the default formula stable names and lets user
    # formulas reference either the canonical or original column.
    metadata = metadata.copy()
    metadata["condition"] = metadata["condition"].astype(str)
    metadata["donor"] = metadata["donor"].astype(str)
    return _PseudobulkInput(
        counts=counts_df,
        metadata=metadata,
        condition_a=condition_a,
        condition_b=condition_b,
        condition_col=condition_col,
        donor_col=donor_col,
        counts_source=counts_source,
        design=design,
        paired=paired,
        covariates=list(covariates),
        donor_ids=sorted(set(donors_a + donors_b)),
        donor_ids_a=donors_a,
        donor_ids_b=donors_b,
    )


def _validate_design_and_make_dds(
    prepared: _PseudobulkInput,
    params: Mapping[str, Any],
    *,
    counts: Optional[pd.DataFrame] = None,
    metadata: Optional[pd.DataFrame] = None,
):
    DeseqDataSet, _ = _require_pydeseq2()
    count_table = prepared.counts if counts is None else counts
    metadata_table = prepared.metadata if metadata is None else metadata
    n_cpus = params.get("n_cpus", params.get("n_processes", 1))
    n_cpus = None if n_cpus is None else max(1, int(n_cpus))
    min_replicates = max(2, int(params.get("min_replicates", 2)))
    try:
        dds = DeseqDataSet(
            counts=count_table,
            metadata=metadata_table,
            design=prepared.design,
            min_replicates=min_replicates,
            n_cpus=n_cpus,
            quiet=True,
            low_memory=bool(params.get("low_memory", False)),
            fit_type=str(params.get("fit_type", "parametric")),
            size_factors_fit_type=str(params.get("size_factors_fit_type", "ratio")),
        )
    except TypeError:
        # Keep compatibility with an older PyDESeq2 release in the optional
        # range if it lacks one of the newer keyword arguments.
        dds = DeseqDataSet(
            counts=count_table,
            metadata=metadata_table,
            design=prepared.design,
            min_replicates=min_replicates,
            n_cpus=n_cpus,
            quiet=True,
        )
    design_matrix = np.asarray(dds.obsm.get("design_matrix"), dtype=float)
    if design_matrix.ndim != 2 or design_matrix.shape[0] != len(metadata_table):
        raise AdvancedStatisticsInputError("PyDESeq2 did not produce a valid design matrix")
    if not np.isfinite(design_matrix).all():
        raise AdvancedStatisticsInputError("design matrix contains non-finite values")
    rank = int(np.linalg.matrix_rank(design_matrix))
    if rank < design_matrix.shape[1]:
        raise AdvancedStatisticsInputError(
            f"design matrix is rank deficient (rank {rank} < {design_matrix.shape[1]} columns); "
            "remove a redundant covariate or revise the formula"
        )
    if design_matrix.shape[0] <= rank:
        raise AdvancedStatisticsInputError(
            f"design has no residual degrees of freedom ({design_matrix.shape[0]} samples, rank {rank})"
        )
    prepared.design_rank = rank
    prepared.design_columns = tuple(str(v) for v in dds.obsm["design_matrix"].columns)
    return dds


def _resolve_contrast(prepared: _PseudobulkInput, params: Mapping[str, Any]) -> list[str] | np.ndarray:
    supplied = params.get("contrast")
    if supplied is None:
        return ["condition", prepared.condition_a, prepared.condition_b]
    if isinstance(supplied, np.ndarray):
        if supplied.ndim != 1:
            raise AdvancedStatisticsInputError("numeric contrast must be one-dimensional")
        return supplied.astype(float)
    if isinstance(supplied, (list, tuple)):
        if len(supplied) == 3 and all(isinstance(v, str) for v in supplied):
            if str(supplied[0]) not in {"condition", prepared.condition_col}:
                raise AdvancedStatisticsInputError(
                    "categorical contrast must target the condition factor"
                )
            if str(supplied[1]) not in {prepared.condition_a, prepared.condition_b} or str(supplied[2]) not in {
                prepared.condition_a,
                prepared.condition_b,
            }:
                raise AdvancedStatisticsInputError("contrast levels must be observed condition levels")
            if str(supplied[1]) == str(supplied[2]):
                raise AdvancedStatisticsInputError("contrast tested and reference levels must differ")
            return [str(v) for v in supplied]
        try:
            return np.asarray(supplied, dtype=float)
        except (TypeError, ValueError) as exc:
            raise AdvancedStatisticsInputError(
                "contrast must be ['factor', 'tested', 'reference'] or a numeric contrast vector"
            ) from exc
    raise AdvancedStatisticsInputError(
        "contrast must be ['factor', 'tested', 'reference'] or a numeric contrast vector"
    )


def _run_pydeseq2(
    prepared: _PseudobulkInput,
    params: Mapping[str, Any],
    *,
    counts: Optional[pd.DataFrame] = None,
    metadata: Optional[pd.DataFrame] = None,
) -> tuple[pd.DataFrame, Any]:
    _, DeseqStats = _require_pydeseq2()
    dds = _validate_design_and_make_dds(prepared, params, counts=counts, metadata=metadata)
    dds.deseq2()
    contrast = _resolve_contrast(prepared, params)
    n_cpus = params.get("n_cpus", params.get("n_processes", 1))
    n_cpus = None if n_cpus is None else max(1, int(n_cpus))
    stats_obj = DeseqStats(
        dds,
        contrast=contrast,
        alpha=float(params.get("alpha", 0.05)),
        cooks_filter=bool(params.get("cooks_filter", True)),
        independent_filter=bool(params.get("independent_filter", True)),
        quiet=True,
        n_cpus=n_cpus,
    )
    stats_obj.summary()
    result = stats_obj.results_df.copy()
    result.index = result.index.astype(str)
    # PyDESeq2 0.4/0.5 uses the columns below; retain a clear error for an
    # incompatible future API instead of falling back to a different test.
    required = {"baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"}
    missing = required.difference(result.columns)
    if missing:
        raise AdvancedStatisticsDependencyError(
            f"unsupported PyDESeq2 results API; missing columns {sorted(missing)}"
        )
    return result, stats_obj


def _result_table(
    result: pd.DataFrame,
    prepared: _PseudobulkInput,
    params: Mapping[str, Any],
    *,
    status: str = "estimated",
) -> pd.DataFrame:
    out = pd.DataFrame(index=result.index.astype(str))
    out["gene"] = out.index
    out["condition_a"] = prepared.condition_a
    out["condition_b"] = prepared.condition_b
    out["base_mean"] = pd.to_numeric(result["baseMean"], errors="coerce").to_numpy()
    out["log2_fold_change"] = pd.to_numeric(result["log2FoldChange"], errors="coerce").to_numpy()
    out["lfc_se"] = pd.to_numeric(result["lfcSE"], errors="coerce").to_numpy()
    alpha = float(params.get("alpha", 0.05))
    if not 0.0 < alpha < 1.0:
        raise AdvancedStatisticsInputError("alpha must be between zero and one")
    z_critical = float(stats.norm.ppf(1.0 - alpha / 2.0))
    out["ci_low"] = out["log2_fold_change"] - z_critical * out["lfc_se"]
    out["ci_high"] = out["log2_fold_change"] + z_critical * out["lfc_se"]
    out["statistic"] = pd.to_numeric(result["stat"], errors="coerce").to_numpy()
    out["p_value"] = pd.to_numeric(result["pvalue"], errors="coerce").to_numpy()
    out["fdr_q_value"] = pd.to_numeric(result["padj"], errors="coerce").to_numpy()
    out["padj"] = out["fdr_q_value"]
    out["effect_definition"] = "PyDESeq2 log2 fold change for condition_a versus condition_b"
    # Keep the legacy auditor's accepted pseudobulk label while exposing the
    # more precise donor-condition unit in a separate field.
    out["statistical_unit"] = "donor_pseudobulk"
    out["pseudobulk_unit"] = "donor_condition"
    out["donor_col"] = prepared.donor_col
    out["counts_source"] = prepared.counts_source
    out["design"] = prepared.design
    out["paired"] = prepared.paired
    out["covariates"] = ",".join(prepared.covariates)
    out["n_pseudobulk_samples"] = int(len(prepared.metadata))
    out["n_donors_condition_a"] = int(len(prepared.donor_ids_a))
    out["n_donors_condition_b"] = int(len(prepared.donor_ids_b))
    finite_lfc = out["log2_fold_change"].notna()
    finite_p = out["p_value"].notna()
    finite_q = out["fdr_q_value"].notna()
    row_status = np.full(len(out), status, dtype=object)
    row_reason = np.full(len(out), "", dtype=object)
    row_status[~(finite_lfc & finite_p)] = "not_estimated"
    row_reason[~finite_lfc] = "log2_fold_change_not_estimated_by_pydeseq2"
    row_reason[finite_lfc & ~finite_p] = "p_value_not_estimated_or_cooks_filtered_by_pydeseq2"
    row_status[finite_lfc & finite_p & ~finite_q] = "estimated_fdr_missing"
    row_reason[finite_lfc & finite_p & ~finite_q] = "adjusted_p_value_not_estimated_by_independent_filtering"
    out["status"] = row_status
    out["status_reason"] = row_reason
    # NaN p-values are preserved as not-estimated.  They are never converted
    # to one, because doing so would present an invented p-value to users.
    significance = out["fdr_q_value"].lt(float(params.get("alpha", 0.05)))
    out["significant_fdr05"] = pd.Series(
        pd.array(significance.where(out["fdr_q_value"].notna(), pd.NA), dtype="boolean"),
        index=out.index,
    )
    out = out.reset_index(drop=True)
    return out


def _fit_metrics(prepared: _PseudobulkInput, result_table: pd.DataFrame) -> dict[str, Any]:
    finite_p = result_table["p_value"].notna()
    finite_q = result_table["fdr_q_value"].notna()
    return {
        "n_genes": int(len(result_table)),
        "n_genes_with_p_value": int(finite_p.sum()),
        "n_genes_with_fdr": int(finite_q.sum()),
        "n_genes_not_estimated": int(result_table["status"].eq("not_estimated").sum()),
        "n_genes_fdr_missing": int(result_table["status"].eq("estimated_fdr_missing").sum()),
        "n_significant_fdr05": int(result_table["significant_fdr05"].eq(True).sum()),
        "condition_a": prepared.condition_a,
        "condition_b": prepared.condition_b,
        "condition_col": prepared.condition_col,
        "donor_col": prepared.donor_col,
        "paired": bool(prepared.paired),
        "covariates": list(prepared.covariates),
        "design": prepared.design,
        "design_rank": prepared.design_rank,
        "design_columns": list(prepared.design_columns),
        "n_pseudobulk_samples": int(len(prepared.metadata)),
        "n_donors_condition_a": int(len(prepared.donor_ids_a)),
        "n_donors_condition_b": int(len(prepared.donor_ids_b)),
        "counts_source": prepared.counts_source,
        "all_genes_retained": True,
        "statistical_unit": "donor_pseudobulk",
        "is_pseudobulk": True,
    }


def _output_uri(contract: TaskContract, input_uri: str, default_path: str) -> str:
    if contract.expected_outputs:
        return str(contract.expected_outputs[0])
    return default_path.format(study_id=ArtifactURI.parse(input_uri).study_id)


def _register_table(
    registry: ArtifactRegistry,
    uri: str,
    table: pd.DataFrame,
    input_uris: Sequence[str],
    contract: TaskContract,
    operation: str,
    parameters: Mapping[str, Any],
    summary_metrics: Mapping[str, Any],
) -> None:
    study_id = ArtifactURI.parse(uri).study_id
    registry.register(
        uri_str=uri,
        payload=table,
        artifact_type=ArtifactType.TABLE,
        study_id=study_id,
        created_by_task=contract.task_id,
        operation=operation,
        parent_uris=list(input_uris),
        parameters=dict(parameters),
        software_versions={
            "pydeseq2": _software_version("pydeseq2"),
            "decoupler": _software_version("decoupler"),
        },
        summary_metrics=dict(summary_metrics),
    )


class PyDESeq2PseudobulkCapability(BaseCapability):
    """Fit a PyDESeq2 donor-condition pseudobulk differential model."""

    CONTRACT_OPERATIONS = (
        "validate_raw_integer_counts",
        "validate_metadata",
        "validate_replicate_qualification",
        "aggregate_counts_by_donor_condition",
        "validate_design_rank",
        "fit_pydeseq2_negative_binomial",
        "wald_contrast",
        "benjamini_hochberg_fdr",
        "retain_all_genes",
        "skip_insufficient_replicates",
    )

    def __init__(self, implementation_id: str = "pydeseq2_pseudobulk_v1"):
        super().__init__(
            capability_name="deg",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
            suitable_for=["raw_count_differential_expression", "paired_donor_design", "covariate_adjustment"],
        )
        self.legacy_aliases = {
            "pydeseq2_deg_v1": self.implementation_id,
            "pydeseq2_donor_pseudobulk_v1": self.implementation_id,
        }
        self.contract_operations = list(self.CONTRACT_OPERATIONS)

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        if not contract.input_artifacts:
            raise AdvancedStatisticsInputError("DEG requires one input AnnData/SCData artifact")
        input_uri = contract.input_artifacts[0]
        _, payload = registry.get(input_uri)
        data = _coerce_data(payload)
        prepared = _prepare_pseudobulk(data, contract.parameters)
        result, _ = _run_pydeseq2(prepared, contract.parameters)
        table = _result_table(result, prepared, contract.parameters)
        metrics = _fit_metrics(prepared, table)
        out_uri = _output_uri(contract, input_uri, "table://{study_id}/pydeseq2_deg/v1")
        _register_table(
            registry,
            out_uri,
            table,
            [input_uri],
            contract,
            "differential_expression_pydeseq2_donor_pseudobulk",
            {
                "condition_a": prepared.condition_a,
                "condition_b": prepared.condition_b,
                "condition_col": prepared.condition_col,
                "donor_col": prepared.donor_col,
                "paired": prepared.paired,
                "covariates": prepared.covariates,
                "design": prepared.design,
                "contrast": _contrast_for_record(prepared, contract.parameters),
                "counts_source": prepared.counts_source,
            },
            metrics,
        )
        operations = [
            "validate_raw_integer_counts",
            "validate_metadata",
            "validate_replicate_qualification",
            "aggregate_counts_by_donor_condition",
            "validate_design_rank",
            "fit_pydeseq2_negative_binomial",
            "wald_contrast",
            "benjamini_hochberg_fdr",
            "retain_all_genes",
        ]
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[input_uri],
            output_artifacts=[out_uri],
            executed_operations=operations,
            metrics=metrics,
        )


def _contrast_for_record(prepared: _PseudobulkInput, params: Mapping[str, Any]) -> Any:
    contrast = _resolve_contrast(prepared, params)
    if isinstance(contrast, np.ndarray):
        return contrast.tolist()
    return list(contrast)


def _read_local_network(path_value: Any) -> tuple[pd.DataFrame, str, str]:
    path_text = str(path_value)
    if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", path_text):
        raise AdvancedStatisticsInputError(
            "network_path must point to a local file; remote network downloads are disabled"
        )
    path = Path(path_text).expanduser()
    if not path.is_file():
        raise AdvancedStatisticsInputError(f"local network file does not exist: {path}")
    suffix = path.suffix.lower()
    if suffix in {".csv"}:
        network = pd.read_csv(path)
    elif suffix in {".tsv", ".txt"}:
        network = pd.read_csv(path, sep="\t")
    elif suffix in {".parquet", ".pq"}:
        network = pd.read_parquet(path)
    elif suffix in {".json"}:
        network = pd.read_json(path)
    else:
        raise AdvancedStatisticsInputError(
            f"unsupported local network format {suffix!r}; use CSV, TSV, JSON, or Parquet"
        )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return network, str(path.resolve()), digest


def _network_from_inputs(
    contract: TaskContract,
    registry: ArtifactRegistry,
    loaded_payloads: Sequence[Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    params = contract.parameters
    network: Any = params.get("network")
    source_kind = "parameters.network"
    resource_hash: Optional[str] = None
    resource_path: Optional[str] = None
    if network is None and params.get("network_path") is not None:
        network, resource_path, resource_hash = _read_local_network(params.get("network_path"))
        source_kind = "local_file"
    if network is None:
        for payload in loaded_payloads[1:]:
            if isinstance(payload, pd.DataFrame) and {"source", "target"}.issubset(payload.columns):
                network = payload
                source_kind = "table_artifact"
                break
            if isinstance(payload, Mapping) and {"source", "target"}.issubset(payload):
                network = pd.DataFrame(payload)
                source_kind = "table_artifact"
                break
    if network is None:
        raise AdvancedStatisticsInputError(
            "a local prior network is required; provide network_path, parameters.network, "
            "or a network table artifact as an additional input"
        )
    if isinstance(network, (str, Path)):
        network, resource_path, resource_hash = _read_local_network(network)
        source_kind = "local_file"
    elif isinstance(network, Mapping):
        network = pd.DataFrame(network)
    elif not isinstance(network, pd.DataFrame):
        raise AdvancedStatisticsInputError("network must be a DataFrame, mapping, local path, or table artifact")
    network = network.copy(deep=True)
    required = {"source", "target"}
    if not required.issubset(network.columns):
        raise AdvancedStatisticsInputError("network must contain source and target columns")
    if "weight" not in network.columns:
        network["weight"] = 1.0
    if network[["source", "target"]].isna().any().any():
        raise AdvancedStatisticsInputError("network source and target values cannot be missing")
    network["source"] = network["source"].astype(str)
    network["target"] = network["target"].astype(str)
    network["weight"] = pd.to_numeric(network["weight"], errors="coerce")
    if network["weight"].isna().any() or not np.isfinite(network["weight"].to_numpy()).all():
        raise AdvancedStatisticsInputError("network weights must be finite numeric values")
    network = network.loc[(network["source"] != "") & (network["target"] != "")].copy()
    if network.empty:
        raise AdvancedStatisticsInputError("network has no non-empty source-target interactions")
    # The hash is deterministic for an in-memory network and can be verified
    # by the orchestrator when a resource hash is supplied in the contract.
    if resource_hash is None:
        canonical = network.sort_values(["source", "target", "weight"], kind="mergesort").to_csv(index=False).encode()
        resource_hash = hashlib.sha256(canonical).hexdigest()
    expected_hash = params.get("network_sha256")
    # The orchestrator may attach hashes for all *_path resources under the
    # generic external_resource_sha256 mapping.  Verify the network resource
    # itself here so a changed file cannot silently reuse a prior contract.
    external_hashes = params.get("external_resource_sha256")
    if isinstance(external_hashes, Mapping):
        for key in ("network_path", "network", str(params.get("network_path", "")), resource_path or ""):
            if key and key in external_hashes:
                expected_hash = external_hashes[key]
                break
    if expected_hash is not None and str(expected_hash).lower() != resource_hash.lower():
        raise AdvancedStatisticsInputError(
            f"network hash mismatch: network_sha256 does not match the supplied local network ({resource_hash})"
        )
    species = params.get("species")
    if species is None or not str(species).strip():
        raise AdvancedStatisticsInputError("species is required for functional analysis")
    network_source = params.get("network_source")
    network_version = params.get("network_version")
    # Table artifacts may carry provenance from their producer, but explicit
    # contract parameters take precedence and are required for bare tables.
    if network_source is None and source_kind == "table_artifact" and len(contract.input_artifacts) > 1:
        try:
            network_source = registry.get_metadata(contract.input_artifacts[1]).summary_metrics.get("network_source")
            network_version = registry.get_metadata(contract.input_artifacts[1]).summary_metrics.get("network_version")
        except Exception:
            pass
    if network_source is None or not str(network_source).strip():
        raise AdvancedStatisticsInputError("network_source is required and must identify the local prior resource")
    if network_version is None or not str(network_version).strip():
        raise AdvancedStatisticsInputError("network_version is required for functional analysis")
    if "species" in network.columns:
        species_values = {str(v) for v in network["species"].dropna().unique()}
        if species_values and species_values != {str(species)}:
            raise AdvancedStatisticsInputError(
                f"network species {sorted(species_values)} does not match requested species {species!r}"
            )
    provenance = {
        "species": str(species),
        "network_source": str(network_source),
        "network_version": str(network_version),
        "network_sha256": resource_hash,
        "external_resource_sha256": resource_hash,
        "network_resource_kind": source_kind,
        "network_path": resource_path,
        "network_interactions": int(len(network)),
    }
    return network, provenance


def _network_sources_for_overlap(network: pd.DataFrame, genes: Sequence[str], tmin: int) -> pd.DataFrame:
    gene_set = set(map(str, genes))
    overlap = network.loc[network["target"].isin(gene_set)].copy()
    if overlap.empty:
        raise AdvancedStatisticsInputError("the local network has no target genes in the expression matrix")
    source_counts = overlap.groupby("source", observed=True)["target"].nunique()
    keep = source_counts[source_counts >= int(tmin)].index
    overlap = overlap.loc[overlap["source"].isin(keep)].copy()
    if overlap.empty:
        raise AdvancedStatisticsInputError(
            f"no network sources have at least tmin={int(tmin)} target genes in the expression matrix"
        )
    return overlap


def _decoupler_result_frames(
    estimate: Any,
    pvalues: Any,
    index: pd.Index,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    estimates = estimate.copy() if isinstance(estimate, pd.DataFrame) else pd.DataFrame(estimate)
    pvals = pvalues.copy() if isinstance(pvalues, pd.DataFrame) else pd.DataFrame(pvalues)
    if len(estimates) != len(index) or len(pvals) != len(index):
        raise AdvancedStatisticsDependencyError("decoupler ULM returned an unexpected number of observations")
    estimates.index = index
    pvals.index = index
    estimates.columns = estimates.columns.astype(str)
    pvals.columns = pvals.columns.astype(str)
    return estimates, pvals


def _activity_design_matrix(prepared: _PseudobulkInput) -> tuple[pd.DataFrame, np.ndarray, tuple[str, ...], int]:
    """Construct the same formula/design used for donor DEG metadata.

    Formulaic is a PyDESeq2 dependency and is also the parser used by the
    current PyDESeq2 release.  Evaluating the formula once with condition A
    and once with condition B yields an adjusted A-versus-B contrast vector,
    including donor fixed effects, covariates, and interactions supplied by a
    caller.  This avoids silently reporting an unadjusted t-test when the DEG
    contract requested adjustment.
    """
    try:
        from formulaic import model_matrix
    except ImportError as exc:  # pragma: no cover - bundled with PyDESeq2
        raise AdvancedStatisticsDependencyError(
            "formulaic is required for design-adjusted functional activity analysis"
        ) from exc
    metadata = prepared.metadata.copy()
    try:
        design_matrix = pd.DataFrame(model_matrix(prepared.design, metadata), index=metadata.index)
    except Exception as exc:
        raise AdvancedStatisticsInputError(
            f"unable to materialize functional activity design {prepared.design!r}: {exc}"
        ) from exc
    if design_matrix.empty or not np.isfinite(design_matrix.to_numpy(dtype=float)).all():
        raise AdvancedStatisticsInputError("functional activity design matrix is empty or non-finite")
    matrix = design_matrix.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(matrix))
    if rank < matrix.shape[1]:
        raise AdvancedStatisticsInputError(
            f"functional activity design is rank deficient (rank {rank} < {matrix.shape[1]} columns)"
        )
    if matrix.shape[0] <= rank:
        raise AdvancedStatisticsInputError(
            f"functional activity design has no residual degrees of freedom ({matrix.shape[0]} samples, rank {rank})"
        )

    def _counterfactual(level: str) -> pd.DataFrame:
        altered = metadata.copy()
        # Preserve both levels when evaluating a counterfactual.  Formulaic
        # otherwise drops a categorical coefficient when every row in the
        # hypothetical frame has the same condition.
        levels = [prepared.condition_a, prepared.condition_b]
        altered["condition"] = pd.Categorical([level] * len(altered), categories=levels)
        if prepared.condition_col in altered.columns:
            altered[prepared.condition_col] = pd.Categorical([level] * len(altered), categories=levels)
        try:
            return pd.DataFrame(model_matrix(prepared.design, altered), index=metadata.index).reindex(columns=design_matrix.columns, fill_value=0.0)
        except Exception as exc:
            raise AdvancedStatisticsInputError(
                f"unable to construct condition contrast for level {level!r}: {exc}"
            ) from exc

    design_a = _counterfactual(prepared.condition_a).to_numpy(dtype=float)
    design_b = _counterfactual(prepared.condition_b).to_numpy(dtype=float)
    contrast_vector = np.nanmean(design_a - design_b, axis=0)
    if not np.isfinite(contrast_vector).all() or not np.any(np.abs(contrast_vector) > 0):
        raise AdvancedStatisticsInputError("functional activity contrast is empty or non-finite")
    return design_matrix, contrast_vector, tuple(str(v) for v in design_matrix.columns), rank


def _activity_test(
    estimates: pd.DataFrame,
    pvals: pd.DataFrame,
    metadata: pd.DataFrame,
    prepared: _PseudobulkInput,
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    design_matrix, contrast_vector, design_columns, design_rank = _activity_design_matrix(prepared)
    design_values = design_matrix.to_numpy(dtype=float)
    residual_df = int(design_values.shape[0] - design_rank)
    t_critical = float(stats.t.ppf(1.0 - float(alpha) / 2.0, residual_df))
    # The same design is used for every source; fit an OLS activity model
    # rather than silently ignoring paired donor effects or covariates.
    xtx_inv = np.linalg.pinv(design_values.T @ design_values)
    contrast_variance_factor = float(contrast_vector @ xtx_inv @ contrast_vector)
    records: list[dict[str, Any]] = []
    donor_rows: list[dict[str, Any]] = []
    for sample_id in estimates.index:
        row_meta = metadata.loc[sample_id]
        for source in estimates.columns:
            donor_rows.append(
                {
                    "sample_id": str(sample_id),
                    "donor": str(row_meta["donor"]),
                    "condition": str(row_meta["condition"]),
                    "source": str(source),
                    "activity": float(estimates.loc[sample_id, source]),
                    "decoupler_p_value": float(pvals.loc[sample_id, source])
                    if pd.notna(pvals.loc[sample_id, source])
                    else np.nan,
                }
            )
    for source in estimates.columns:
        y = estimates[source].to_numpy(dtype=float)
        if not np.isfinite(y).all():
            effect = standard_error = statistic = p_value = np.nan
            test_name = "not_estimated_nonfinite_decoupler_activity"
        else:
            coefficients, _, _, _ = np.linalg.lstsq(design_values, y, rcond=None)
            residuals = y - design_values @ coefficients
            residual_variance = float(np.sum(residuals**2) / residual_df)
            effect = float(contrast_vector @ coefficients)
            standard_error = float(np.sqrt(max(0.0, residual_variance * contrast_variance_factor)))
            if standard_error > 0 and np.isfinite(standard_error):
                statistic = float(effect / standard_error)
                p_value = float(2.0 * stats.t.sf(abs(statistic), residual_df))
            else:
                statistic = p_value = np.nan
            test_name = "ols_activity_design"
        mask_a = metadata["condition"].astype(str) == prepared.condition_a
        mask_b = metadata["condition"].astype(str) == prepared.condition_b
        values_a = estimates.loc[mask_a, source].to_numpy(dtype=float)
        values_b = estimates.loc[mask_b, source].to_numpy(dtype=float)
        records.append(
            {
                "source": str(source),
                "condition_a": prepared.condition_a,
                "condition_b": prepared.condition_b,
                "mean_activity_condition_a": float(np.nanmean(values_a)) if len(values_a) else np.nan,
                "mean_activity_condition_b": float(np.nanmean(values_b)) if len(values_b) else np.nan,
                "activity_effect_condition_a_vs_b": effect,
                "activity_se": standard_error,
                "activity_ci_low": effect - t_critical * standard_error
                if np.isfinite(effect) and np.isfinite(standard_error)
                else np.nan,
                "activity_ci_high": effect + t_critical * standard_error
                if np.isfinite(effect) and np.isfinite(standard_error)
                else np.nan,
                "statistic": statistic,
                "p_value": p_value,
                "n_donors_condition_a": int(len(values_a)),
                "n_donors_condition_b": int(len(values_b)),
                "paired": bool(prepared.paired),
                "test": test_name,
                "mean_decoupler_p_value_condition_a": float(np.nanmean(pvals.loc[mask_a, source]))
                if pvals.loc[mask_a, source].notna().any()
                else np.nan,
                "mean_decoupler_p_value_condition_b": float(np.nanmean(pvals.loc[mask_b, source]))
                if pvals.loc[mask_b, source].notna().any()
                else np.nan,
                "design": prepared.design,
                "design_rank": design_rank,
                "design_columns": ",".join(design_columns),
                "residual_degrees_of_freedom": residual_df,
                "status": "estimated" if np.isfinite(p_value) else "not_estimated",
                "status_reason": "" if np.isfinite(p_value) else test_name,
            }
        )
    summary = pd.DataFrame(records)
    if not summary.empty:
        # An unestimated raw p-value must not become an invented adjusted 1.0.
        estimated = np.isfinite(summary["p_value"].to_numpy(dtype=float))
        adjusted = np.full(len(summary), np.nan)
        adjusted[estimated] = benjamini_hochberg(summary.loc[estimated, "p_value"].to_numpy(dtype=float))
        summary["fdr_q_value"] = adjusted
        summary["padj"] = summary["fdr_q_value"]
        summary["significant_fdr05"] = pd.Series(
            pd.array(
                summary["fdr_q_value"].lt(0.05).where(summary["fdr_q_value"].notna(), pd.NA),
                dtype="boolean",
            ),
            index=summary.index,
        )
    return summary, pd.DataFrame(donor_rows)


class DecouplerFunctionalAnalysisCapability(BaseCapability):
    """Run decoupler ULM and compare donor-level pathway/TF activities."""

    CONTRACT_OPERATIONS = (
        "validate_raw_integer_counts",
        "validate_metadata",
        "validate_replicate_qualification",
        "validate_local_network_provenance",
        "aggregate_counts_by_donor_condition",
        "log1p_cpm_transform_for_decoupler",
        "run_decoupler_ulm",
        "build_formula_activity_design",
        "compare_donor_activities",
        "benjamini_hochberg_fdr",
        "retain_all_network_sources",
        "skip_insufficient_donors",
    )

    def __init__(self, implementation_id: str = "decoupler_ulm_v2"):
        super().__init__(
            capability_name="functional_activity",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
            output_types=[ArtifactType.TABLE],
            suitable_for=["donor_level_pathway_activity", "donor_level_tf_activity"],
        )
        self.legacy_aliases = {"decoupler_functional_v1": self.implementation_id}
        self.contract_operations = list(self.CONTRACT_OPERATIONS)

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        if not contract.input_artifacts:
            raise AdvancedStatisticsInputError("functional activity requires an expression input artifact")
        loaded: list[Any] = []
        for uri in contract.input_artifacts:
            loaded.append(registry.get(uri)[1])
        data = _coerce_data(loaded[0])
        prepared = _prepare_pseudobulk(data, contract.parameters)
        network, provenance = _network_from_inputs(contract, registry, loaded)
        tmin = max(1, int(contract.parameters.get("tmin", contract.parameters.get("min_network_size", 1))))
        network = _network_sources_for_overlap(network, list(prepared.counts.columns), tmin)
        totals = prepared.counts.sum(axis=1).to_numpy(dtype=float)
        if np.any(totals <= 0):
            raise AdvancedStatisticsInputError("a donor-condition pseudobulk has zero total counts")
        expr = np.log1p(prepared.counts.div(totals, axis=0) * 1_000_000.0)
        dc = _require_decoupler()
        # decoupler >=2 exposes mt.ulm; using it directly keeps this capability
        # tied to the current documented API instead of a removed run_ulm shim.
        try:
            estimate, pvalues = dc.mt.ulm(
                data=expr,
                net=network[["source", "target", "weight"]],
                tmin=tmin,
                verbose=False,
            )
        except Exception as exc:
            raise AdvancedStatisticsDependencyError(f"decoupler ULM failed: {exc}") from exc
        estimates, pvals = _decoupler_result_frames(estimate, pvalues, expr.index)
        summary, donor_activity = _activity_test(
            estimates,
            pvals,
            prepared.metadata,
            prepared,
            alpha=float(contract.parameters.get("alpha", 0.05)),
        )
        feature_kind = str(contract.parameters.get("analysis_kind", contract.parameters.get("network_kind", "functional")))
        summary["feature_type"] = feature_kind
        summary["decoupler_method"] = "ulm"
        for key, value in provenance.items():
            summary[key] = value
        summary["statistical_unit"] = "donor_pseudobulk_activity"
        summary["is_pseudobulk"] = True
        summary["activity_input"] = "log1p_CPM_from_layers_counts"
        out_uri = _output_uri(contract, contract.input_artifacts[0], "table://{study_id}/functional_activity/v1")
        metrics = {
            "n_sources": int(len(summary)),
            "n_sources_not_estimated": int(summary["status"].eq("not_estimated").sum()) if not summary.empty else 0,
            "n_significant_fdr05": int(summary["significant_fdr05"].eq(True).sum()) if not summary.empty else 0,
            "n_donors_condition_a": len(prepared.donor_ids_a),
            "n_donors_condition_b": len(prepared.donor_ids_b),
            "paired": prepared.paired,
            "condition_a": prepared.condition_a,
            "condition_b": prepared.condition_b,
            "network_source": provenance["network_source"],
            "network_version": provenance["network_version"],
            "network_sha256": provenance["network_sha256"],
            "species": provenance["species"],
            "feature_type": feature_kind,
            "design": prepared.design,
            "covariates": list(prepared.covariates),
            "statistical_unit": "donor_pseudobulk_activity",
            "is_pseudobulk": True,
            "p_values_source": "formulaic OLS donor activity comparison; decoupler ULM p-values retained per donor",
            "all_sources_retained": True,
        }
        _register_table(
            registry,
            out_uri,
            summary,
            contract.input_artifacts,
            contract,
            "functional_activity_decoupler_ulm_donor_comparison",
            {**provenance, "tmin": tmin, "analysis_kind": feature_kind},
            metrics,
        )
        outputs = [out_uri]
        # If the contract requests a second output, provide donor-level
        # activities as a transparent companion table.  Never create it when
        # it is not requested so expected_outputs remains authoritative.
        if len(contract.expected_outputs) >= 2:
            donor_uri = str(contract.expected_outputs[1])
            donor_table = donor_activity.copy()
            donor_table["feature_type"] = feature_kind
            donor_table["decoupler_method"] = "ulm"
            for key, value in provenance.items():
                donor_table[key] = value
            _register_table(
                registry,
                donor_uri,
                donor_table,
                contract.input_artifacts,
                contract,
                "functional_activity_decoupler_ulm_donor_scores",
                {**provenance, "tmin": tmin, "analysis_kind": feature_kind},
                {"n_rows": int(len(donor_table)), **metrics},
            )
            outputs.append(donor_uri)
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=list(contract.input_artifacts),
            output_artifacts=outputs,
            executed_operations=[
                "validate_raw_integer_counts",
                "validate_metadata",
                "validate_replicate_qualification",
                "validate_local_network_provenance",
                "aggregate_counts_by_donor_condition",
                "log1p_cpm_transform_for_decoupler",
                "run_decoupler_ulm",
                "build_formula_activity_design",
                "compare_donor_activities",
                "benjamini_hochberg_fdr",
                "retain_all_network_sources",
            ],
            metrics=metrics,
        )


def _skip_table(reason: str, prepared: Optional[_PseudobulkInput], params: Mapping[str, Any]) -> pd.DataFrame:
    values = {
        "gene": None,
        "left_out_donor": None,
        "status": "skipped",
        "skip_reason": reason,
        "condition_a": prepared.condition_a if prepared else params.get("condition_a"),
        "condition_b": prepared.condition_b if prepared else params.get("condition_b"),
        "p_value": np.nan,
        "fdr_q_value": np.nan,
    }
    return pd.DataFrame([values])


def _leave_one_out_summary(
    full_table: pd.DataFrame,
    fit_tables: Sequence[pd.DataFrame],
    n_requested_fits: int,
    alpha: float,
) -> pd.DataFrame:
    """Summarize effect and decision stability for every gene.

    ``direction_consistency`` is the fraction of *estimated* leave-one-out
    effects having the same non-zero sign as the full model.  It is NaN when
    the baseline effect is zero/unknown or no leave-one-out effect was
    estimated.  ``estimated_coverage`` uses all requested fits as its
    denominator, so a partially failed run cannot look complete after
    silently dropping unestimated values.
    """
    loo = pd.concat(list(fit_tables), ignore_index=True) if fit_tables else pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for _, baseline in full_table.iterrows():
        gene = str(baseline["gene"])
        if not loo.empty:
            gene_loo = loo.loc[loo["gene"].astype(str) == gene]
        else:
            gene_loo = pd.DataFrame()
        effects = pd.to_numeric(gene_loo.get("log2_fold_change", pd.Series(dtype=float)), errors="coerce").dropna()
        q_values = pd.to_numeric(gene_loo.get("fdr_q_value", pd.Series(dtype=float)), errors="coerce").dropna()
        n_estimated = int(len(effects))
        baseline_lfc = pd.to_numeric(pd.Series([baseline.get("log2_fold_change")]), errors="coerce").iloc[0]
        baseline_q = pd.to_numeric(pd.Series([baseline.get("fdr_q_value")]), errors="coerce").iloc[0]
        if n_estimated and pd.notna(baseline_lfc) and float(baseline_lfc) != 0.0:
            direction_consistency = float((np.sign(effects.to_numpy()) == np.sign(float(baseline_lfc))).mean())
        else:
            direction_consistency = np.nan
        if len(q_values) and pd.notna(baseline_q):
            baseline_significant = bool(float(baseline_q) < alpha)
            loo_significant = q_values.to_numpy(dtype=float) < alpha
            # Retention is relative to the full-model decision: among the
            # estimated fits, how often did the same significant/non-significant
            # decision remain?
            significance_retention = float((loo_significant == baseline_significant).mean())
            n_significant_loo = int(loo_significant.sum())
        else:
            baseline_significant = pd.NA
            significance_retention = np.nan
            n_significant_loo = 0
        rows.append(
            {
                "gene": gene,
                "baseline_log2_fold_change": float(baseline_lfc) if pd.notna(baseline_lfc) else np.nan,
                "baseline_fdr_q_value": float(baseline_q) if pd.notna(baseline_q) else np.nan,
                "baseline_significant_fdr05": baseline_significant,
                "n_estimated_loo": n_estimated,
                "estimated_coverage": float(n_estimated / n_requested_fits) if n_requested_fits else np.nan,
                "direction_consistency": direction_consistency,
                "min_log2_fold_change": float(effects.min()) if n_estimated else np.nan,
                "max_log2_fold_change": float(effects.max()) if n_estimated else np.nan,
                "median_log2_fold_change": float(effects.median()) if n_estimated else np.nan,
                "n_significant_loo": n_significant_loo,
                "significance_retention": significance_retention,
                "summary_status": "estimated" if n_estimated else "not_estimated",
            }
        )
    return pd.DataFrame(rows)


class PyDESeq2LeaveOneDonorOutCapability(BaseCapability):
    """Assess DEG sensitivity by re-fitting after each donor is removed."""

    CONTRACT_OPERATIONS = (
        "validate_raw_integer_counts",
        "validate_metadata",
        "validate_replicate_qualification",
        "validate_leave_one_donor_qualification",
        "skip_insufficient_donors",
        "fit_full_pydeseq2_model",
        "remove_complete_donor_for_each_fit",
        "refit_pydeseq2_per_donor",
        "retain_all_genes_per_fit",
        "record_fit_failures_without_claiming_robustness",
    )

    def __init__(self, implementation_id: str = "pydeseq2_leave_one_donor_out_v1"):
        super().__init__(
            capability_name="donor_sensitivity",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
            suitable_for=["donor_leave_one_out", "paired_donor_sensitivity"],
        )
        self.legacy_aliases = {"donor_loo_pydeseq2_v1": self.implementation_id}
        self.contract_operations = list(self.CONTRACT_OPERATIONS)

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        if not contract.input_artifacts:
            raise AdvancedStatisticsInputError("donor sensitivity requires one input AnnData/SCData artifact")
        input_uri = contract.input_artifacts[0]
        _, payload = registry.get(input_uri)
        data = _coerce_data(payload)
        # Prepare with the structural minimum so an underpowered LOO request
        # can produce an explicit, auditable skip artifact.  The ordinary
        # DEG/functional capabilities still require their configured minimum
        # donor count; LOO applies its separate ``min_loo_donors`` gate below.
        loo_params = dict(contract.parameters)
        loo_params["min_donors"] = 1
        prepared = _prepare_pseudobulk(data, loo_params)
        min_donors = max(3, int(contract.parameters.get("min_loo_donors", _DEFAULT_LOO_MIN_DONORS)))
        n_a, n_b = len(prepared.donor_ids_a), len(prepared.donor_ids_b)
        if n_a < min_donors or n_b < min_donors:
            reason = (
                f"skipped: leave-one-donor-out requires at least {min_donors} donors per condition; "
                f"found {n_a} and {n_b}"
            )
            table = _skip_table(reason, prepared, contract.parameters)
            out_uri = _output_uri(contract, input_uri, "table://{study_id}/donor_sensitivity/v1")
            metrics = {
                "skipped": True,
                "skip_reason": reason,
                "n_donors_condition_a": n_a,
                "n_donors_condition_b": n_b,
                "paired": prepared.paired,
                "n_successful_fits": 0,
                "n_requested_fits": 0,
                "complete_fit_coverage": 0.0,
                "all_requested_fits_completed": False,
                "scientific_robustness_claim_supported": False,
                "scientific_robustness_limit": "leave-one-donor-out is a sensitivity analysis, not independent validation",
            }
            _register_table(
                registry,
                out_uri,
                table,
                [input_uri],
                contract,
                "pydeseq2_leave_one_donor_out_skipped",
                {"skip_reason": reason, "paired": prepared.paired},
                metrics,
            )
            outputs = [out_uri]
            if len(contract.expected_outputs) >= 2:
                summary_uri = str(contract.expected_outputs[1])
                summary_table = table.copy()
                summary_table["complete_fit_coverage"] = 0.0
                summary_table["scientific_robustness_claim_supported"] = False
                _register_table(
                    registry,
                    summary_uri,
                    summary_table,
                    [input_uri],
                    contract,
                    "pydeseq2_leave_one_donor_out_summary_skipped",
                    {"skip_reason": reason, "paired": prepared.paired},
                    {"n_rows": int(len(summary_table)), **metrics},
                )
                outputs.append(summary_uri)
            return TaskResult(
                task_id=contract.task_id,
                status=TaskStatus.SUCCESS,
                capability=self.capability_name,
                method_used=self.implementation_id,
                input_artifacts=[input_uri],
                output_artifacts=outputs,
                executed_operations=["validate_leave_one_donor_qualification", "skip_insufficient_donors"],
                metrics=metrics,
            )

        # Fit the complete model first.  This is also used as the explicit
        # baseline row for every gene in the output.
        full_result, _ = _run_pydeseq2(prepared, contract.parameters)
        full_table = _result_table(full_result, prepared, contract.parameters)
        fit_tables: list[pd.DataFrame] = []
        successful_fits = 0
        failed_fits: dict[str, str] = {}
        # Paired designs remove both condition rows for one donor by masking on
        # the canonical donor column.  Independent designs remove one donor
        # sample at a time from its observed condition.
        donors_to_remove = sorted(set(prepared.donor_ids_a + prepared.donor_ids_b))
        for donor in donors_to_remove:
            keep = prepared.metadata["donor"].astype(str) != str(donor)
            sub_counts = prepared.counts.loc[keep].copy()
            sub_metadata = prepared.metadata.loc[keep].copy()
            try:
                sub_result, _ = _run_pydeseq2(
                    prepared,
                    contract.parameters,
                    counts=sub_counts,
                    metadata=sub_metadata,
                )
                # _run_pydeseq2 updates the design diagnostics on prepared;
                # every successful result contains the complete gene table.
                sub_table = _result_table(sub_result, prepared, contract.parameters)
                sub_table.insert(1, "left_out_donor", str(donor))
                sub_table["fit_status"] = "estimated"
                sub_table["removed_sample_count"] = int((~keep).sum())
                fit_tables.append(sub_table)
                successful_fits += 1
            except Exception as exc:
                failed_fits[str(donor)] = f"{type(exc).__name__}: {exc}"
        baseline = full_table.copy()
        baseline.insert(1, "left_out_donor", "__full_model__")
        baseline["fit_status"] = "full_model"
        baseline["removed_sample_count"] = 0
        output_table = pd.concat([baseline] + fit_tables, ignore_index=True)
        summary_table = _leave_one_out_summary(
            full_table,
            fit_tables,
            n_requested_fits=len(donors_to_remove),
            alpha=float(contract.parameters.get("alpha", 0.05)),
        )
        # Keep the detailed donor-by-gene table as the primary artifact, while
        # repeating the per-gene stability summaries on each row for callers
        # that request only one output URI.
        output_table = output_table.merge(summary_table, on="gene", how="left")
        output_table["n_requested_fits"] = len(donors_to_remove)
        output_table["n_successful_fits"] = successful_fits
        output_table["paired_leaveout_removes_complete_donor"] = bool(prepared.paired)
        output_table["complete_fit_coverage"] = float(successful_fits / len(donors_to_remove)) if donors_to_remove else np.nan
        output_table["all_requested_fits_completed"] = bool(successful_fits == len(donors_to_remove))
        # Fit coverage and effect consistency are descriptive sensitivity
        # metrics.  The capability never promotes them to a scientific
        # robustness claim; that requires independent evidence and review.
        output_table["scientific_robustness_claim_supported"] = False
        out_uri = _output_uri(contract, input_uri, "table://{study_id}/donor_sensitivity/v1")
        metrics = {
            "skipped": False,
            "skip_reason": None,
            "n_requested_fits": int(len(donors_to_remove)),
            "n_successful_fits": int(successful_fits),
            "n_failed_fits": int(len(failed_fits)),
            "failed_fits": failed_fits,
            "n_donors_condition_a": n_a,
            "n_donors_condition_b": n_b,
            "paired": prepared.paired,
            "paired_leaveout_removes_complete_donor": bool(prepared.paired),
            "complete_fit_coverage": float(successful_fits / len(donors_to_remove)) if donors_to_remove else np.nan,
            "all_requested_fits_completed": bool(successful_fits == len(donors_to_remove)),
            "scientific_robustness_claim_supported": False,
            "scientific_robustness_limit": "leave-one-donor-out is a sensitivity analysis, not independent validation",
            "n_genes_with_direction_consistency": int(summary_table["direction_consistency"].notna().sum()),
            "n_genes_with_significance_retention": int(summary_table["significance_retention"].notna().sum()),
            "n_genes": int(len(full_table)),
            "all_genes_retained_per_successful_fit": True,
        }
        _register_table(
            registry,
            out_uri,
            output_table,
            [input_uri],
            contract,
            "pydeseq2_leave_one_donor_out",
            {
                "paired": prepared.paired,
                "min_loo_donors": min_donors,
                "condition_a": prepared.condition_a,
                "condition_b": prepared.condition_b,
                "design": prepared.design,
            },
            metrics,
        )
        outputs = [out_uri]
        if len(contract.expected_outputs) >= 2:
            summary_uri = str(contract.expected_outputs[1])
            summary_output = summary_table.copy()
            summary_output["n_requested_fits"] = len(donors_to_remove)
            summary_output["n_successful_fits"] = successful_fits
            summary_output["paired_leaveout_removes_complete_donor"] = bool(prepared.paired)
            summary_output["complete_fit_coverage"] = float(successful_fits / len(donors_to_remove)) if donors_to_remove else np.nan
            summary_output["all_requested_fits_completed"] = bool(successful_fits == len(donors_to_remove))
            summary_output["scientific_robustness_claim_supported"] = False
            _register_table(
                registry,
                summary_uri,
                summary_output,
                [input_uri],
                contract,
                "pydeseq2_leave_one_donor_out_summary",
                {
                    "paired": prepared.paired,
                    "min_loo_donors": min_donors,
                    "condition_a": prepared.condition_a,
                    "condition_b": prepared.condition_b,
                    "design": prepared.design,
                },
                {"n_rows": int(len(summary_output)), **metrics},
            )
            outputs.append(summary_uri)
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[input_uri],
            output_artifacts=outputs,
            executed_operations=[
                "validate_raw_integer_counts",
                "validate_metadata",
                "validate_replicate_qualification",
                "validate_leave_one_donor_qualification",
                "fit_full_pydeseq2_model",
                "remove_complete_donor_for_each_fit",
                "refit_pydeseq2_per_donor",
                "retain_all_genes_per_fit",
                "record_fit_failures_without_claiming_robustness",
            ],
            metrics=metrics,
        )


# Public aliases make the intent discoverable to callers that use the longer
# capability names while preserving the stable registry class names above.
PyDESeq2DifferentialExpressionCapability = PyDESeq2PseudobulkCapability
AdvancedDifferentialExpressionCapability = PyDESeq2PseudobulkCapability
DecouplerActivityCapability = DecouplerFunctionalAnalysisCapability
FunctionalActivityCapability = DecouplerFunctionalAnalysisCapability
DonorSensitivityCapability = PyDESeq2LeaveOneDonorOutCapability
LeaveOneDonorOutCapability = PyDESeq2LeaveOneDonorOutCapability


__all__ = [
    "AdvancedStatisticsDependencyError",
    "AdvancedStatisticsInputError",
    "PyDESeq2PseudobulkCapability",
    "PyDESeq2DifferentialExpressionCapability",
    "AdvancedDifferentialExpressionCapability",
    "DecouplerFunctionalAnalysisCapability",
    "DecouplerActivityCapability",
    "FunctionalActivityCapability",
    "PyDESeq2LeaveOneDonorOutCapability",
    "DonorSensitivityCapability",
    "LeaveOneDonorOutCapability",
]
