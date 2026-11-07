#!/usr/bin/env python3
"""Describe donor sensitivity for a public Kang cell-type contrast.

This harness consumes the immutable, already accepted Kang subset and its full
PyDESeq2 result table. It re-fits the same contrast after removing each paired
donor, audits the resulting artifacts independently, and writes durable
diagnostics. The run is a sensitivity analysis, not independent validation or
a claim of novel discovery.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path, PurePath, PureWindowsPath
from typing import Any

PROTOCOL_ALPHA = 0.05
PROTOCOL_MIN_ABS_LFC = 1.0
EXPECTED_DONORS = 8
EXPECTED_GENES = 15_706
KNOWN_POSITIVE_MARKERS = ("ISG15", "IFI6", "LY6E", "MX1", "IFIT1", "IFIT2", "ISG20")
MARKER_REFERENCE_URL = "https://satijalab.org/seurat/articles/integration_introduction.html"
PROJECT_ROOT = Path(__file__).resolve().parents[1]

OUTPUT_NAMES = (
    "full_baseline.tsv",
    "loo_primary.tsv",
    "loo_summary.tsv",
    "gene_sensitivity_summary.tsv",
    "candidate_result_cards.tsv",
    "known_marker_result_cards.tsv",
    "audit.json",
    "contract_environment.json",
    "run_summary.json",
    "report.md",
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_external_path(work_dir: PurePath, *, project_root: PurePath) -> PurePath:
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute external path")
    if isinstance(work_dir, PureWindowsPath) != isinstance(project_root, PureWindowsPath):
        raise ValueError("--work-dir and repository root must use the same path style")
    try:
        work_dir.relative_to(project_root)
    except ValueError:
        pass
    else:
        raise ValueError("--work-dir must be outside the repository/workspace")
    if isinstance(work_dir, PureWindowsPath) and work_dir.drive.casefold() == "c:":
        raise ValueError("--work-dir must not use the C: drive")
    return work_dir


def _validate_external_work_dir(work_dir: Path, *, project_root: Path = PROJECT_ROOT) -> Path:
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute external path")
    resolved = work_dir.expanduser().resolve()
    _validate_external_path(resolved, project_root=project_root.resolve())
    return resolved


def _json_safe(value: Any) -> Any:
    import numpy as np
    import pandas as pd

    if hasattr(value, "model_dump"):
        return _json_safe(value.model_dump(mode="json"))
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if pd is not None and value is pd.NA:
        return None
    return str(value)


def _write_json(payload: Any, path: Path) -> None:
    path.write_text(
        json.dumps(_json_safe(payload), ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _sha256_assertions(
    *, data_sha256: str, manifest_data_sha256: str, manifest_file_sha256: str,
    source_summary: dict[str, Any], full_de_sha256: str
) -> dict[str, Any]:
    files = source_summary.get("files", {})
    expected_data = (files.get("data") or {}).get("sha256")
    expected_manifest = (files.get("manifest") or {}).get("sha256")
    expected_full_de = (files.get("table") or {}).get("sha256")
    checks = {
        "data_matches_manifest_declaration": data_sha256 == manifest_data_sha256,
        "data_matches_source_summary": data_sha256 == expected_data,
        "manifest_file_matches_source_summary": manifest_file_sha256 == expected_manifest,
        "full_de_matches_source_summary": full_de_sha256 == expected_full_de,
    }
    if not all(checks.values()):
        raise ValueError(
            "Input provenance SHA256 mismatch: "
            + json.dumps(
                {
                    "checks": checks,
                    "data_sha256": data_sha256,
                    "manifest_data_sha256": manifest_data_sha256,
                    "manifest_file_sha256": manifest_file_sha256,
                    "source_summary_manifest_sha256": expected_manifest,
                    "source_summary_data_sha256": expected_data,
                    "full_de_sha256": full_de_sha256,
                    "source_summary_full_de_sha256": expected_full_de,
                },
                sort_keys=True,
            )
        )
    return {
        "checks": checks,
        "data_sha256": data_sha256,
        "manifest_data_sha256": manifest_data_sha256,
        "manifest_file_sha256": manifest_file_sha256,
        "full_de_sha256": full_de_sha256,
    }


def _validate_inputs(
    *, manifest: Any, source_summary: dict[str, Any], adata: Any, full_table: Any
) -> dict[str, Any]:
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if manifest.dataset_kind != "real":
        raise ValueError("Only the declared real-data Kang subset is accepted")
    if source_summary.get("dataset_kind") != "real" or source_summary.get("study_id") != manifest.study_id:
        raise ValueError("Source summary dataset kind/study_id does not match the real-data manifest")
    if source_summary.get("overall_status") != "passed":
        raise ValueError("The accepted Kang input run must have overall_status='passed'")
    if source_summary.get("files", {}).get("data", {}).get("sha256") != manifest.sha256:
        raise ValueError("Source summary and manifest do not declare the same derived H5AD SHA256")
    source_record = source_summary.get("source", {})
    if (
        source_record.get("url") != manifest.source_url
        or source_record.get("citation") != manifest.source
        or source_record.get("license") != manifest.license
    ):
        raise ValueError("Source citation, URL, or license differs between the passed summary and manifest")
    lock_summary = source_summary.get("lockfile", {})
    expected_env = manifest.reference_environment
    if expected_env is None or lock_summary.get("sha256") != expected_env.lockfile_sha256:
        raise ValueError("Source summary lock hash does not match the manifest's fixed environment")
    if manifest.statistics.alpha != PROTOCOL_ALPHA:
        raise ValueError(f"This predeclared sensitivity protocol requires alpha={PROTOCOL_ALPHA}")
    if not manifest.donor.paired:
        raise ValueError("This Kang protocol requires the declared paired-donor design")
    if manifest.contrast.tested != "stim" or manifest.contrast.reference != "ctrl":
        raise ValueError("The declared Kang contrast must remain stim versus ctrl")
    if manifest.contrast.column not in adata.obs or manifest.donor.column not in adata.obs:
        raise ValueError("Subset is missing the manifest condition or donor column")
    if manifest.counts_layer not in adata.layers:
        raise ValueError(f"Subset lacks the manifest raw-count layer {manifest.counts_layer!r}")
    if bool(adata.uns.get("is_simulated", False)):
        raise ValueError("Subset AnnData is marked is_simulated")
    if adata.n_vars != EXPECTED_GENES:
        raise ValueError(f"Expected all {EXPECTED_GENES} genes, observed {adata.n_vars}")
    declared_subset = source_summary.get("subset", {})
    if int(declared_subset.get("n_genes", -1)) != int(adata.n_vars):
        raise ValueError("Source summary gene count does not match the subset H5AD")
    if int(declared_subset.get("n_cells", -1)) != int(adata.n_obs):
        raise ValueError("Source summary cell count does not match the subset H5AD")
    if declared_subset.get("all_genes_and_matching_cells_retained") is not True:
        raise ValueError("Source summary does not attest to retaining all matching cells and genes")

    conditions = adata.obs[manifest.contrast.column].astype("string")
    donors = adata.obs[manifest.donor.column].astype("string")
    if conditions.isna().any() or donors.isna().any():
        raise ValueError("Subset has missing condition or donor metadata")
    if set(conditions.astype(str).unique()) != {manifest.contrast.tested, manifest.contrast.reference}:
        raise ValueError("Subset condition values differ from the declared ctrl/stim contrast")
    donors_by_condition = {
        level: set(donors.loc[conditions.astype(str) == level].astype(str).unique())
        for level in (manifest.contrast.tested, manifest.contrast.reference)
    }
    if donors_by_condition[manifest.contrast.tested] != donors_by_condition[manifest.contrast.reference]:
        raise ValueError("Kang conditions do not contain the same paired donors")
    if len(donors_by_condition[manifest.contrast.tested]) != EXPECTED_DONORS:
        raise ValueError(f"Expected {EXPECTED_DONORS} paired donors, observed {len(donors_by_condition[manifest.contrast.tested])}")

    counts = adata.layers[manifest.counts_layer]
    values = counts.data if sparse.issparse(counts) else np.asarray(counts).reshape(-1)
    if not np.issubdtype(values.dtype, np.integer):
        if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
            raise ValueError("Manifest counts layer must contain finite raw integer counts")
    if (values < 0).any():
        raise ValueError("Manifest counts layer contains negative values")

    genes = adata.var["gene_name"].astype(str).tolist() if "gene_name" in adata.var else adata.var_names.astype(str).tolist()
    if len(genes) != len(set(genes)):
        raise ValueError("Subset gene names are not unique")
    required_columns = {"gene", "log2_fold_change", "ci_low", "ci_high", "p_value", "fdr_q_value"}
    missing_columns = sorted(required_columns.difference(full_table.columns))
    if missing_columns:
        raise ValueError(f"Full DE table is missing result columns: {missing_columns}")
    if full_table["gene"].isna().any() or full_table["gene"].astype(str).duplicated().any():
        raise ValueError("Full DE table gene identifiers are missing or duplicated")
    if set(full_table["gene"].astype(str)) != set(genes):
        raise ValueError("Full DE table gene universe does not exactly match the H5AD gene universe")
    if len(full_table) != EXPECTED_GENES:
        raise ValueError(f"Expected one full DE row for each of {EXPECTED_GENES} genes")
    if "condition_a" in full_table and set(full_table["condition_a"].dropna().astype(str)) != {manifest.contrast.tested}:
        raise ValueError("Full DE table tested level differs from the manifest")
    if "condition_b" in full_table and set(full_table["condition_b"].dropna().astype(str)) != {manifest.contrast.reference}:
        raise ValueError("Full DE table reference level differs from the manifest")

    cell_type = declared_subset.get("cell_type")
    if not cell_type:
        raise ValueError("Source summary lacks the declared Kang cell type")
    cell_type_col = str(declared_subset.get("cell_type_column", "cell_type"))
    if cell_type_col in adata.obs and set(adata.obs[cell_type_col].dropna().astype(str).unique()) != {str(cell_type)}:
        raise ValueError("H5AD cell type values do not match the source summary")
    return {
        "study_id": manifest.study_id,
        "cell_type": str(cell_type),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "donors": sorted(donors_by_condition[manifest.contrast.tested]),
        "condition_counts": {str(key): int(value) for key, value in conditions.value_counts().sort_index().items()},
        "counts_layer": manifest.counts_layer,
        "raw_integer_counts": True,
        "data_origin": "real_public_dataset",
    }


def _numeric_series(frame: Any, column: str, index: Any):
    import pandas as pd

    if column not in frame:
        return pd.Series(float("nan"), index=index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").reindex(index)


def _summarize_primary(primary: Any, expected_donors: list[str], alpha: float = PROTOCOL_ALPHA):
    """Independently derive gene/fold counts and impact values from primary rows."""
    import numpy as np
    import pandas as pd

    required = {"gene", "left_out_donor", "log2_fold_change", "fdr_q_value"}
    missing = sorted(required.difference(primary.columns))
    if missing:
        raise ValueError(f"LOO primary output is missing columns: {missing}")
    if primary[["gene", "left_out_donor"]].isna().any().any():
        raise ValueError("LOO primary output contains missing gene or donor identifiers")
    if primary.duplicated(["gene", "left_out_donor"]).any():
        raise ValueError("LOO primary output contains duplicate gene/donor rows")
    donor_values = set(primary["left_out_donor"].astype(str).unique())
    unexpected = donor_values.difference({"__full_model__", *expected_donors})
    if unexpected:
        raise ValueError(f"LOO primary output contains unrequested donors: {sorted(unexpected)}")
    full = primary.loc[primary["left_out_donor"].astype(str) == "__full_model__"].copy()
    folds = primary.loc[primary["left_out_donor"].astype(str) != "__full_model__"].copy()
    if full.empty or full["gene"].astype(str).duplicated().any():
        raise ValueError("LOO primary output must contain one complete full-model row per gene")
    full["gene"] = full["gene"].astype(str)
    folds["gene"] = folds["gene"].astype(str)
    folds["left_out_donor"] = folds["left_out_donor"].astype(str)
    full_by_gene = full.set_index("gene", drop=False)
    fold_groups = {gene: group.set_index("left_out_donor", drop=False) for gene, group in folds.groupby("gene", sort=False)}
    rows: list[dict[str, Any]] = []
    for gene, baseline in full_by_gene.iterrows():
        group = fold_groups.get(gene, pd.DataFrame()).reindex(expected_donors)
        lfc = pd.to_numeric(group.get("log2_fold_change"), errors="coerce")
        qvals = pd.to_numeric(group.get("fdr_q_value"), errors="coerce")
        if not isinstance(lfc, pd.Series):
            lfc = pd.Series(np.nan, index=expected_donors, dtype=float)
        if not isinstance(qvals, pd.Series):
            qvals = pd.Series(np.nan, index=expected_donors, dtype=float)
        lfc.index = expected_donors
        qvals.index = expected_donors
        finite_lfc = pd.Series(np.isfinite(lfc.to_numpy(dtype=float)), index=expected_donors)
        finite_q = pd.Series(np.isfinite(qvals.to_numpy(dtype=float)), index=expected_donors)
        baseline_lfc = pd.to_numeric(pd.Series([baseline.get("log2_fold_change")]), errors="coerce").iloc[0]
        baseline_q = pd.to_numeric(pd.Series([baseline.get("fdr_q_value")]), errors="coerce").iloc[0]
        baseline_finite_lfc = pd.notna(baseline_lfc) and math.isfinite(float(baseline_lfc))
        baseline_finite_q = pd.notna(baseline_q) and math.isfinite(float(baseline_q))
        sign_defined = baseline_finite_lfc and float(baseline_lfc) != 0.0
        fold_signs = np.sign(lfc.fillna(0.0).to_numpy(dtype=float))
        baseline_sign = np.sign(float(baseline_lfc)) if baseline_finite_lfc else 0.0
        same_sign = finite_lfc & (fold_signs == baseline_sign)
        if sign_defined:
            same_sign = pd.Series(same_sign, index=expected_donors) & finite_lfc
        else:
            same_sign = pd.Series(False, index=expected_donors)
        direction_flip = finite_lfc & (fold_signs == -baseline_sign) if sign_defined else pd.Series(False, index=expected_donors)
        zero_effect = finite_lfc & (fold_signs == 0) if sign_defined else pd.Series(False, index=expected_donors)
        q_significant = finite_q & qvals.lt(alpha)
        baseline_significant = bool(float(baseline_q) < alpha) if baseline_finite_q else None
        q_status_changed = (
            finite_q & q_significant.ne(baseline_significant)
            if baseline_significant is not None
            else pd.Series(False, index=expected_donors)
        )
        if baseline_finite_lfc:
            deltas = (lfc - float(baseline_lfc)).abs().where(finite_lfc)
        else:
            deltas = pd.Series(float("nan"), index=expected_donors)
        finite_deltas = deltas.dropna()
        if finite_deltas.empty:
            max_delta = None
            influential_donor = None
        else:
            max_delta = float(finite_deltas.max())
            tied = sorted(str(donor) for donor in finite_deltas.index[finite_deltas == finite_deltas.max()])
            influential_donor = tied[0]
        missing_effect = [str(donor) for donor in expected_donors if not bool(finite_lfc.loc[donor])]
        missing_q = [str(donor) for donor in expected_donors if not bool(finite_q.loc[donor])]
        flips = [str(donor) for donor in expected_donors if bool(direction_flip.loc[donor])]
        zero_effects = [str(donor) for donor in expected_donors if bool(zero_effect.loc[donor])]
        changed_q = [str(donor) for donor in expected_donors if bool(q_status_changed.loc[donor])]
        rows.append({
            "gene": str(gene),
            "full_log2_fold_change": float(baseline_lfc) if baseline_finite_lfc else None,
            "full_ci_low": _finite_or_none(baseline.get("ci_low")),
            "full_ci_high": _finite_or_none(baseline.get("ci_high")),
            "full_fdr_q_value": float(baseline_q) if baseline_finite_q else None,
            "alpha": float(alpha),
            "expected_donor_folds": len(expected_donors),
            "n_finite_loo_effects": int(finite_lfc.sum()),
            "finite_folds_of_8": f"{int(finite_lfc.sum())}/{len(expected_donors)}",
            "finite_loo_effect_donors": [str(donor) for donor in expected_donors if bool(finite_lfc.loc[donor])],
            "n_same_sign_loo_effects": int(same_sign.sum()) if sign_defined else None,
            "same_sign_folds_of_8": f"{int(same_sign.sum())}/{len(expected_donors)}" if sign_defined else None,
            "same_sign_loo_donors": [str(donor) for donor in expected_donors if bool(same_sign.loc[donor])],
            "n_loo_q_below_alpha": int(q_significant.sum()),
            "q_below_alpha_folds_of_8": f"{int(q_significant.sum())}/{len(expected_donors)}",
            "n_finite_loo_q_values": int(finite_q.sum()),
            "baseline_significant_at_alpha": baseline_significant,
            "n_q_status_changes": len(changed_q),
            "q_status_changed_donors": changed_q,
            "direction_flip_donors": flips,
            "zero_effect_donors": zero_effects,
            "missing_effect_donors": missing_effect,
            "missing_q_donors": missing_q,
            "max_abs_delta_lfc": max_delta,
            "most_influential_donor": influential_donor,
            "loo_coverage_status": "complete" if len(missing_effect) == 0 else "incomplete",
            "direction_status": (
                "undefined_full_effect" if not sign_defined
                else "direction_stable_8_of_8" if len(missing_effect) == 0 and int(same_sign.sum()) == len(expected_donors)
                else "direction_flip_and_incomplete_coverage" if flips and missing_effect
                else "direction_flip" if flips
                else "incomplete_coverage" if missing_effect
                else "zero_effect_nonstable" if zero_effects
                else "direction_stable_8_of_8"
            ),
        })
    return pd.DataFrame(rows)


def _finite_or_none(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _candidate_cards(summary_table: Any):
    import pandas as pd

    q = pd.to_numeric(summary_table["full_fdr_q_value"], errors="coerce")
    lfc = pd.to_numeric(summary_table["full_log2_fold_change"], errors="coerce")
    known_symbols = {gene.upper() for gene in KNOWN_POSITIVE_MARKERS}
    exploratory = ~summary_table["gene"].astype(str).str.upper().isin(known_symbols)
    selected = summary_table.loc[q.lt(PROTOCOL_ALPHA) & lfc.abs().ge(PROTOCOL_MIN_ABS_LFC) & exploratory].copy()
    selected["classification"] = "post_selection_exploratory_sensitivity"
    selected["full_model_candidate_rule"] = "full-model q < 0.05 and abs(log2 fold change) >= 1"
    selected["significance_interpretation"] = "Descriptive threshold changes only; not a false-positive verdict."
    selected["independent_validation"] = False
    selected["novelty_claimed"] = False
    selected["priority_reason"] = selected.apply(_priority_reason, axis=1)
    rank = {
        "direction_flip_and_incomplete_coverage": 0,
        "direction_flip": 0,
        "incomplete_coverage": 1,
        "zero_effect_nonstable": 2,
        "q_status_change": 3,
        "direction_stable_without_q_status_change": 4,
    }
    selected["_priority"] = selected["priority_reason"].map(rank).fillna(4)
    selected = selected.sort_values(["_priority", "max_abs_delta_lfc", "full_fdr_q_value", "gene"], ascending=[True, False, True, True], na_position="last")
    return selected.drop(columns=["_priority"]).reset_index(drop=True)


def _priority_reason(row: Any) -> str:
    direction = str(row.get("direction_status", ""))
    if "direction_flip" in direction:
        return direction
    if direction == "incomplete_coverage":
        return direction
    if direction == "zero_effect_nonstable":
        return direction
    if int(row.get("n_q_status_changes", 0) or 0) > 0:
        return "q_status_change"
    return "direction_stable_without_q_status_change"


def _known_marker_cards(full_table: Any, sensitivity_summary: Any, source_summary: dict[str, Any]):
    import pandas as pd

    by_gene = {str(row["gene"]).upper(): row for _, row in full_table.iterrows()}
    marker_report = source_summary.get("marker_direction_check") or {}
    report_by_gene = {
        str(row.get("gene", "")).upper(): row
        for row in marker_report.get("markers", [])
        if isinstance(row, dict)
    }
    sensitivity_by_gene = {
        str(row["gene"]).upper(): row
        for _, row in sensitivity_summary.iterrows()
    }
    rows = []
    for marker in KNOWN_POSITIVE_MARKERS:
        full = by_gene.get(marker.upper())
        report = report_by_gene.get(marker.upper(), {})
        sensitivity = sensitivity_by_gene.get(marker.upper())
        rows.append({
            "gene": marker,
            "classification": "known_positive_benchmark",
            "expected_direction": "positive log2 fold change for stim versus ctrl",
            "full_log2_fold_change": _finite_or_none(full.get("log2_fold_change")) if full is not None else None,
            "full_ci_low": _finite_or_none(full.get("ci_low")) if full is not None else None,
            "full_ci_high": _finite_or_none(full.get("ci_high")) if full is not None else None,
            "full_fdr_q_value": _finite_or_none(full.get("fdr_q_value")) if full is not None else None,
            "finite_folds_of_8": sensitivity.get("finite_folds_of_8") if sensitivity is not None else "0/8",
            "same_sign_folds_of_8": sensitivity.get("same_sign_folds_of_8") if sensitivity is not None else None,
            "q_below_alpha_folds_of_8": sensitivity.get("q_below_alpha_folds_of_8") if sensitivity is not None else "0/8",
            "direction_status": sensitivity.get("direction_status") if sensitivity is not None else "missing_from_loo",
            "max_abs_delta_lfc": sensitivity.get("max_abs_delta_lfc") if sensitivity is not None else None,
            "most_influential_donor": sensitivity.get("most_influential_donor") if sensitivity is not None else None,
            "direction_flip_donors": sensitivity.get("direction_flip_donors") if sensitivity is not None else [],
            "missing_effect_donors": sensitivity.get("missing_effect_donors") if sensitivity is not None else [],
            "n_q_status_changes": sensitivity.get("n_q_status_changes") if sensitivity is not None else None,
            "observed_direction_status": report.get("status", "missing_from_source_marker_report"),
            "source_marker_check_overall_status": marker_report.get("status"),
            "reference": marker_report.get("reference", MARKER_REFERENCE_URL),
            "interpretation": "Predeclared known-positive benchmark only; direction is not a novelty or discovery test.",
        })
    return pd.DataFrame(rows)


def _compare_full_baseline(full_table: Any, primary: Any, *, rtol: float, atol: float) -> dict[str, Any]:
    import numpy as np
    import pandas as pd

    required = ("log2_fold_change", "ci_low", "ci_high", "p_value", "fdr_q_value")
    full_rows = primary.loc[primary["left_out_donor"].astype(str) == "__full_model__"].copy()
    if full_rows["gene"].astype(str).duplicated().any():
        return {"status": "failed", "reason": "duplicate full-model gene rows"}
    left = full_table.set_index(full_table["gene"].astype(str)).sort_index()
    right = full_rows.set_index(full_rows["gene"].astype(str)).sort_index()
    if not left.index.equals(right.index):
        return {"status": "failed", "reason": "full-model and source full-DE gene identifiers differ"}
    metrics = {}
    failures = []
    for column in required:
        a = pd.to_numeric(left[column], errors="coerce").to_numpy(dtype=float)
        b = pd.to_numeric(right[column], errors="coerce").to_numpy(dtype=float)
        finite_mask_equal = np.array_equal(np.isfinite(a), np.isfinite(b))
        nan_mask_equal = np.array_equal(np.isnan(a), np.isnan(b))
        both = np.isfinite(a) & np.isfinite(b)
        delta = np.abs(a[both] - b[both])
        max_abs = float(delta.max()) if delta.size else 0.0
        passed = finite_mask_equal and nan_mask_equal and (not delta.size or bool(np.all(np.isclose(a[both], b[both], rtol=rtol, atol=atol))))
        metrics[column] = {
            "passed": passed,
            "finite_mask_equal": finite_mask_equal,
            "nan_mask_equal": nan_mask_equal,
            "n_finite_compared": int(both.sum()),
            "max_absolute_difference": max_abs,
            "rtol": float(rtol),
            "atol": float(atol),
        }
        if not passed:
            failures.append(column)
    return {"status": "passed" if not failures else "failed", "failed_columns": failures, "metrics": metrics}


def _check_fold_row_coverage(primary: Any, expected_donors: list[str], expected_genes: int) -> dict[str, Any]:
    """Require a complete gene table per fit while allowing valid statistic NA values."""
    import pandas as pd

    required = {"gene", "left_out_donor", "fit_status"}
    if not required.issubset(primary.columns):
        return {"status": "failed", "reason": f"missing columns: {sorted(required.difference(primary.columns))}"}
    full = primary.loc[primary["left_out_donor"].astype(str) == "__full_model__"]
    full_genes = set(full["gene"].astype(str))
    expected_set = set(expected_donors)
    folds = primary.loc[primary["left_out_donor"].astype(str) != "__full_model__"]
    actual_donors = set(folds["left_out_donor"].astype(str).unique())
    unexpected = sorted(actual_donors.difference(expected_set))
    donor_metrics: dict[str, Any] = {}
    for donor in expected_donors:
        rows = folds.loc[folds["left_out_donor"].astype(str) == donor]
        observed_genes = set(rows["gene"].astype(str))
        duplicate_rows = int(rows["gene"].astype(str).duplicated().sum())
        donor_metrics[donor] = {
            "n_rows": int(len(rows)),
            "n_genes": int(len(observed_genes)),
            "missing_gene_count": int(len(full_genes.difference(observed_genes))),
            "extra_gene_count": int(len(observed_genes.difference(full_genes))),
            "duplicate_gene_rows": duplicate_rows,
            "all_fit_status_estimated": bool(rows["fit_status"].astype(str).eq("estimated").all()) if not rows.empty else False,
            "complete_gene_table": bool(
                len(full_genes) == expected_genes
                and len(rows) == expected_genes
                and observed_genes == full_genes
                and duplicate_rows == 0
                and rows["fit_status"].astype(str).eq("estimated").all()
            ),
        }
    full_status_ok = bool(
        len(full) == expected_genes
        and len(full_genes) == expected_genes
        and full["fit_status"].astype(str).eq("full_model").all()
    )
    passed = full_status_ok and not unexpected and all(row["complete_gene_table"] for row in donor_metrics.values())
    return {
        "status": "passed" if passed else "failed",
        "expected_gene_rows_per_fit": int(expected_genes),
        "full_model_complete_gene_table": full_status_ok,
        "unexpected_donors": unexpected,
        "donors": donor_metrics,
        "non_finite_statistics_are_allowed_when_rows_are_present": True,
    }


def _build_report(summary: dict[str, Any]) -> str:
    context = summary["input_context"]
    fit = summary.get("capability", {})
    audit = summary.get("audit", {})
    cards = summary.get("candidate_counts", {})
    examples = summary.get("candidate_examples", {})
    integrity = summary.get("checks", {}).get("source_unchanged", {})
    lines = [
        "# Public Kang donor-sensitivity analysis",
        "",
        "## Question and scope",
        "",
        f"For **{context.get('cell_type', 'the declared cell type')}**, which stim-versus-ctrl gene effects change materially when one paired donor is removed?",
        "",
        "This run re-fits all genes for each of the eight paired donor removals. It is a sensitivity analysis, not independent validation, and does not establish a novel biological discovery. FDR is interpreted within this cell type's genome-wide stim-versus-ctrl family; no cross-cell-type specificity claim is made.",
        "",
        "## Reproducibility and execution",
        "",
        f"- Study: `{context.get('study_id')}`; cells: {context.get('n_cells')}; genes: {context.get('n_genes')}; paired donors: {len(context.get('donors', []))}.",
        f"- Capability: `{fit.get('method_used')}`; requested folds: {fit.get('metrics', {}).get('n_requested_fits')}; successful folds: {fit.get('metrics', {}).get('n_successful_fits')}; failed folds: {fit.get('metrics', {}).get('n_failed_fits')}.",
        f"- Independent ScientificAuditor: `{audit.get('status', 'not_run')}`.",
        f"- Source H5AD unchanged: `{integrity.get('unchanged')}`; machine elapsed seconds: {summary.get('machine_elapsed_seconds')}.",
        "- Human active minutes: `null`; human-time savings claim supported: `false`.",
        "",
        "## Predeclared result interpretation",
        "",
        f"Full-model candidates are selected using q < {PROTOCOL_ALPHA} and |log2 fold change| >= {PROTOCOL_MIN_ABS_LFC}; these are post-selection exploratory sensitivity results. Direction stability requires a finite LOO effect with the same nonzero sign in all 8 of 8 folds. q-threshold changes are reported separately and are not treated as false-positive verdicts.",
        f"- Candidate cards: {cards.get('n_candidates', 0)}.",
        f"- Candidates with any direction flip: {cards.get('n_direction_flip', 0)}.",
        f"- Candidates with incomplete effect coverage: {cards.get('n_incomplete_coverage', 0)}.",
        f"- Candidates with any q-threshold change: {cards.get('n_q_status_change', 0)}.",
        "- The seven Seurat tutorial markers are retained in a separate known-positive benchmark table and excluded from exploratory candidate cards.",
        "",
        "## Exploratory sensitivity examples",
        "",
    ]
    for key, title in (
        ("direction_flip", "Direction flips"),
        ("incomplete_coverage", "Incomplete effect coverage"),
        ("q_status_change", "q-threshold changes"),
    ):
        rows = examples.get(key, [])
        if rows:
            rendered = "; ".join(
                f"`{row['gene']}` (full LFC {row.get('full_log2_fold_change')}, "
                f"folds {row.get('finite_folds_of_8')}, same sign {row.get('same_sign_folds_of_8')}, "
                f"q<alpha {row.get('q_below_alpha_folds_of_8')}, max |ΔLFC| {row.get('max_abs_delta_lfc')}, "
                f"most influential donor {row.get('most_influential_donor')})"
                for row in rows[:5]
            )
            lines.append(f"- {title}: {rendered}.")
    if not any(examples.get(key) for key in ("direction_flip", "incomplete_coverage", "q_status_change")):
        lines.append("No direction-flip, incomplete-coverage, or q-threshold-change exploratory candidates were observed under the predeclared criteria.")
    lines.extend([
        "## Files",
        "",
        "`full_baseline.tsv` contains the validated input DE table. `loo_primary.tsv` and `loo_summary.tsv` preserve the capability artifacts. `gene_sensitivity_summary.tsv` independently summarizes every gene from primary fold rows; candidate and known-marker cards are separate. `audit.json`, `contract_environment.json`, and `run_summary.json` retain audit, design/environment, and execution details.",
        "",
        "A missing estimate or failed donor fit is recorded as incomplete; no missing fold is imputed or silently removed. See the TSV and JSON artifacts for donor-level details.",
        "",
    ])
    return "\n".join(lines)


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path, help="Accepted Kang cell-type subset H5AD")
    parser.add_argument("--manifest", required=True, type=Path, help="The subset's acceptance manifest JSON")
    parser.add_argument("--full-de", required=True, type=Path, help="Full-model PyDESeq2 TSV from the subset run")
    parser.add_argument("--source-summary", required=True, type=Path, help="Kang subset run summary JSON")
    parser.add_argument("--work-dir", required=True, type=Path, help="New absolute external directory for durable artifacts")
    return parser.parse_args(argv)


def _preflight_output_dir(work_dir: Path) -> None:
    work_dir.mkdir(parents=True, exist_ok=True)
    collisions = [str(work_dir / name) for name in OUTPUT_NAMES if (work_dir / name).exists()]
    registry_dir = work_dir / "artifact_registry"
    if registry_dir.exists() and any(registry_dir.iterdir()):
        collisions.append(str(registry_dir))
    if collisions:
        raise FileExistsError("Refusing to overwrite existing research-value outputs: " + ", ".join(collisions))


def run(args: argparse.Namespace) -> int:
    started = time.perf_counter()
    work_dir = _validate_external_work_dir(args.work_dir)
    _preflight_output_dir(work_dir)
    manifest_path = args.manifest.expanduser().resolve()
    data_path = args.data.expanduser().resolve()
    full_de_path = args.full_de.expanduser().resolve()
    source_summary_path = args.source_summary.expanduser().resolve()
    summary: dict[str, Any] = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "study_id": None,
        "scientific_discovery_claimed": False,
        "independent_validation": False,
        "human_active_minutes": None,
        "human_savings_claim_supported": False,
        "output_directory": str(work_dir),
        "input_paths": {
            "data": str(data_path), "manifest": str(manifest_path),
            "full_de": str(full_de_path), "source_summary": str(source_summary_path),
        },
        "checks": {},
    }
    try:
        for path in (manifest_path, data_path, full_de_path, source_summary_path):
            if not path.is_file():
                raise FileNotFoundError(f"Required input does not exist: {path}")
        data_sha_before = _sha256_file(data_path)
        manifest_table, manifest_dir = __import__("eacbp.acceptance.manifest", fromlist=["load_manifest"]).load_manifest(manifest_path)
        source_summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
        table_sha = _sha256_file(full_de_path)
        provenance = _sha256_assertions(
            data_sha256=data_sha_before,
            manifest_data_sha256=manifest_table.sha256,
            manifest_file_sha256=_sha256_file(manifest_path),
            source_summary=source_summary,
            full_de_sha256=table_sha,
        )
        import anndata as ad
        import pandas as pd

        adata = ad.read_h5ad(data_path)
        full_table = pd.read_csv(full_de_path, sep="\t", na_values=["NA"], keep_default_na=True)
        input_context = _validate_inputs(manifest=manifest_table, source_summary=source_summary, adata=adata, full_table=full_table)
        summary.update({"study_id": manifest_table.study_id, "input_context": input_context})
        summary["checks"]["input_provenance"] = provenance

        from eacbp.acceptance.runner import _check_reference_environment, _parameters
        environment = _check_reference_environment(manifest_table, manifest_dir)
        summary["checks"]["reference_environment"] = environment
        if environment.get("status") != "passed":
            raise RuntimeError("Fixed scientific runtime/lock verification failed")
        lock_summary = source_summary.get("lockfile", {})
        if lock_summary.get("sha256") != manifest_table.reference_environment.lockfile_sha256:
            raise RuntimeError("Source summary and manifest fixed-lock SHA256 differ")

        expected_donors = input_context["donors"]
        contract_parameters = {
            **_parameters(manifest_table),
            "study_id": manifest_table.study_id,
            "min_loo_donors": EXPECTED_DONORS,
            "scientific_question": (
                f"For {input_context['cell_type']}, which stim-versus-ctrl gene effects are sensitive to removing one paired donor?"
            ),
            "target_cell_type": input_context["cell_type"],
            "target_population": f"{input_context['cell_type']} cells from the eight paired Kang donors",
            "fdr_family": f"kang_2018:{input_context['cell_type']}:stim_vs_ctrl_all_genes",
            "sensitivity_plan": ["leave_one_donor_out", "remove_each_complete_paired_donor"],
        }
        from eacbp.artifact.registry import ArtifactRegistry
        from eacbp.capabilities.advanced_statistics import PyDESeq2LeaveOneDonorOutCapability
        from eacbp.schemas.artifact import ArtifactType
        from eacbp.schemas.task import TaskContract, TaskStatus

        registry_path = work_dir / "artifact_registry"
        registry = ArtifactRegistry(str(registry_path))
        input_uri = f"adata://{manifest_table.study_id}/kang_subset/v1"
        primary_uri = f"table://{manifest_table.study_id}/donor_sensitivity/v1"
        summary_uri = f"table://{manifest_table.study_id}/donor_sensitivity_summary/v1"
        registry.register(
            input_uri,
            adata,
            ArtifactType.ANNDATA,
            manifest_table.study_id,
            "verify_public_research_value",
            "accepted_kang_raw_counts_with_counts_layer",
            summary_metrics={
                "dataset_kind": "real",
                "source_sha256": data_sha_before,
                "manifest_sha256": _sha256_file(manifest_path),
                "source_summary_sha256": _sha256_file(source_summary_path),
                "full_de_sha256": table_sha,
                "counts_layer": manifest_table.counts_layer,
            },
        )
        contract = TaskContract(
            task_id=f"{manifest_table.study_id}_donor_sensitivity",
            capability="donor_sensitivity",
            method="pydeseq2_leave_one_donor_out_v1",
            input_artifacts=[input_uri],
            expected_outputs=[primary_uri, summary_uri],
            validation_requirements=[
                "advanced_statistics_integrity",
                "multiple_testing_correction",
                "pseudoreplication_audit",
            ],
            parameters=contract_parameters,
        )
        capability_result = PyDESeq2LeaveOneDonorOutCapability().execute(contract, registry)
        result_record = _json_safe(capability_result.model_dump(mode="json"))
        summary["capability"] = result_record
        primary = registry.load_payload(primary_uri) if primary_uri in capability_result.output_artifacts else None
        capability_summary = registry.load_payload(summary_uri) if summary_uri in capability_result.output_artifacts else None
        capability_ok = capability_result.status == TaskStatus.SUCCESS and not bool(capability_result.metrics.get("skipped"))
        summary["checks"]["capability_completed"] = {
            "status": "passed" if capability_ok else "failed",
            "task_status": capability_result.status.value,
            "skip_reason": capability_result.metrics.get("skip_reason"),
        }
        if primary is None or capability_summary is None:
            raise RuntimeError("Donor LOO capability did not publish both required durable output artifacts")
        row_coverage = _check_fold_row_coverage(primary, expected_donors, EXPECTED_GENES)
        primary.to_csv(work_dir / "loo_primary.tsv", sep="\t", index=False, na_rep="NA")
        capability_summary.to_csv(work_dir / "loo_summary.tsv", sep="\t", index=False, na_rep="NA")
        full_table.to_csv(work_dir / "full_baseline.tsv", sep="\t", index=False, na_rep="NA")
        gene_summary = _summarize_primary(primary, expected_donors)
        if set(gene_summary["gene"].astype(str)) != set(full_table["gene"].astype(str)):
            raise RuntimeError("LOO full-model rows do not cover the same all-gene universe as the supplied baseline")
        gene_summary.to_csv(work_dir / "gene_sensitivity_summary.tsv", sep="\t", index=False, na_rep="NA")
        candidate_cards = _candidate_cards(gene_summary)
        candidate_cards.to_csv(work_dir / "candidate_result_cards.tsv", sep="\t", index=False, na_rep="NA")
        known_cards = _known_marker_cards(full_table, gene_summary, source_summary)
        known_cards.to_csv(work_dir / "known_marker_result_cards.tsv", sep="\t", index=False, na_rep="NA")

        baseline_check = _compare_full_baseline(
            full_table, primary, rtol=manifest_table.statistics.rtol, atol=manifest_table.statistics.atol
        )
        summary["checks"]["full_model_matches_supplied_baseline"] = baseline_check
        from eacbp.auditor import ScientificAuditor

        audit_report = ScientificAuditor().audit_task(contract, capability_result, registry)
        audit_payload = _json_safe(audit_report.model_dump(mode="json"))
        audit_payload["status"] = "passed" if audit_report.overall_passed else "failed"
        (work_dir / "audit.json").write_text(
            json.dumps(audit_payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        fit_metrics = capability_result.metrics
        fold_coverage_passed = (
            capability_ok
            and int(fit_metrics.get("n_requested_fits", 0)) == EXPECTED_DONORS
            and int(fit_metrics.get("n_successful_fits", 0)) == EXPECTED_DONORS
            and int(fit_metrics.get("n_failed_fits", 0)) == 0
            and row_coverage.get("status") == "passed"
        )
        summary["checks"]["all_eight_donor_folds_complete"] = {
            "status": "passed" if fold_coverage_passed else "failed",
            "expected_donors": expected_donors,
            "n_requested_fits": fit_metrics.get("n_requested_fits"),
            "n_successful_fits": fit_metrics.get("n_successful_fits"),
            "n_failed_fits": fit_metrics.get("n_failed_fits"),
            "failed_fits": fit_metrics.get("failed_fits", {}),
            "fold_row_coverage": row_coverage,
            "genes_with_complete_finite_lfc_coverage": int((gene_summary["n_finite_loo_effects"] == EXPECTED_DONORS).sum()),
            "genes_with_incomplete_finite_lfc_coverage": int((gene_summary["n_finite_loo_effects"] != EXPECTED_DONORS).sum()),
            "non_estimable_gene_values_do_not_fail_fit_coverage": True,
        }
        summary["audit"] = {"status": audit_payload["status"], "overall_passed": bool(audit_report.overall_passed)}
        summary["candidate_counts"] = {
            "n_candidates": int(len(candidate_cards)),
            "n_direction_flip": int(candidate_cards["direction_status"].astype(str).str.contains("direction_flip").sum()) if not candidate_cards.empty else 0,
            "n_incomplete_coverage": int(candidate_cards["loo_coverage_status"].ne("complete").sum()) if not candidate_cards.empty else 0,
            "n_q_status_change": int(candidate_cards["n_q_status_changes"].fillna(0).astype(int).gt(0).sum()) if not candidate_cards.empty else 0,
        }
        example_columns = [
            "gene", "full_log2_fold_change", "full_fdr_q_value", "finite_folds_of_8",
            "same_sign_folds_of_8", "q_below_alpha_folds_of_8", "direction_status",
            "n_q_status_changes", "max_abs_delta_lfc", "most_influential_donor",
        ]
        summary["candidate_examples"] = {
            "direction_flip": candidate_cards.loc[
                candidate_cards["direction_status"].astype(str).str.contains("direction_flip")
            ].head(5)[example_columns].to_dict(orient="records") if not candidate_cards.empty else [],
            "incomplete_coverage": candidate_cards.loc[
                candidate_cards["loo_coverage_status"].ne("complete")
            ].head(5)[example_columns].to_dict(orient="records") if not candidate_cards.empty else [],
            "q_status_change": candidate_cards.loc[
                candidate_cards["n_q_status_changes"].fillna(0).astype(int).gt(0)
            ].head(5)[example_columns].to_dict(orient="records") if not candidate_cards.empty else [],
        }
        contract_environment = {
            "manifest": _json_safe(manifest_table.model_dump(mode="json")),
            "environment_check": environment,
            "task_contract": _json_safe(contract.model_dump(mode="json")),
            "resolved_inference_contract": _json_safe(capability_result.inference_contract.model_dump(mode="json")) if capability_result.inference_contract else None,
            "runtime_versions": environment.get("actual", {}),
            "protocol": {
                "alpha": PROTOCOL_ALPHA,
                "candidate_minimum_abs_log2_fold_change": PROTOCOL_MIN_ABS_LFC,
                "candidate_rule": "full-model q < 0.05 and abs(log2 fold change) >= 1",
                "known_positive_markers_excluded_from_exploratory_cards": list(KNOWN_POSITIVE_MARKERS),
                "direction_stable_rule": "finite LOO effect with same nonzero sign in all 8 of 8 paired-donor folds",
                "q_threshold_changes_are_descriptive_only": True,
                "fdr_family_scope": "one genome-wide DE contrast within this cell type",
                "independent_validation": False,
                "novelty_claimed": False,
            },
        }
        _write_json(contract_environment, work_dir / "contract_environment.json")

        data_sha_after = _sha256_file(data_path)
        source_unchanged = data_sha_before == data_sha_after
        summary["checks"]["source_unchanged"] = {
            "path": str(data_path),
            "sha256_before": data_sha_before,
            "sha256_after": data_sha_after,
            "unchanged": source_unchanged,
            "status": "passed" if source_unchanged else "failed",
        }
        summary["baseline_comparison"] = baseline_check
        summary["checks"]["independent_audit"] = summary["audit"]
        elapsed = time.perf_counter() - started
        summary["machine_elapsed_seconds"] = elapsed
        summary["human_active_minutes"] = None
        summary["human_savings_claim_supported"] = False
        summary["overall_status"] = "passed" if all((
            environment.get("status") == "passed",
            capability_ok,
            fold_coverage_passed,
            baseline_check.get("status") == "passed",
            audit_report.overall_passed,
            source_unchanged,
        )) else "failed"
        summary["status"] = summary["overall_status"]
        report = _build_report(summary)
        (work_dir / "report.md").write_text(report, encoding="utf-8")
        _write_json(summary, work_dir / "run_summary.json")
        print(json.dumps(_json_safe(summary), ensure_ascii=False, indent=2, allow_nan=False), flush=True)
        return 0 if summary["overall_status"] == "passed" else 1
    except Exception as exc:
        summary["status"] = "failed"
        summary["overall_status"] = "failed"
        summary["error"] = {"type": type(exc).__name__, "message": str(exc)}
        summary["machine_elapsed_seconds"] = time.perf_counter() - started
        summary["human_active_minutes"] = None
        summary["human_savings_claim_supported"] = False
        if not (work_dir / "run_summary.json").exists():
            _write_json(summary, work_dir / "run_summary.json")
        if not (work_dir / "report.md").exists():
            (work_dir / "report.md").write_text(
                "# Public Kang donor-sensitivity analysis\n\nRun failed before a complete analysis was available.\n\n"
                f"Error: `{type(exc).__name__}: {exc}`\n",
                encoding="utf-8",
            )
        print(f"verify_public_research_value failed: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        return 2


def main(argv: list[str] | None = None) -> int:
    return run(_parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
