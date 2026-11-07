"""Bounded, pre-execution QC perturbations and descriptive paired comparisons.

These helpers never choose a winning branch or change the primary analysis.
Callers execute branches through the normal audited study lifecycle and supply
only admitted QC, cluster and statistical results. Reanalysing the same cells
does not provide independent biological replication or validate cell damage.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import re
from typing import Any, Mapping

import numpy as np
import pandas as pd


QC_SENSITIVITY_POLICY_VERSION = "bounded_scaled_mad_v1"
_PARAMETERS = ("min_genes", "max_mito_pct")
_REASONS = ("retained", "low_gene_count", "high_mito_percentage",
            "low_gene_count_and_high_mito_percentage")


def _fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def _config_scope(config: Mapping[str, Any]) -> dict[str, Any]:
    clean = deepcopy(dict(config))
    settings = clean.get("autonomous_analysis", {})
    for name in ("qc_sensitivity_plan", "qc_sensitivity_branch"):
        settings.pop(name, None)
    return clean


def qc_config_sha256(config: Mapping[str, Any]) -> str:
    """Bind native analysis configuration while excluding self-referential receipts."""
    return _fingerprint(_config_scope(config))


def _normalized_groups(groups: pd.Series, ids: list[str]) -> list[str | None]:
    if not groups.index.is_unique or set(map(str, groups.index)) != set(ids):
        raise ValueError("Observed group identities must cover the original cell universe uniquely")
    aligned = groups.copy()
    aligned.index = aligned.index.astype(str)
    values = aligned.reindex(ids).astype("string").str.strip()
    return [None if pd.isna(value) or not value else str(value) for value in values]


def _metric_identity(ids, detected, totals, mito, available) -> str:
    return _fingerprint({"ids": list(ids), "n_genes": np.asarray(detected, dtype=int).tolist(),
                         "total_counts": np.asarray(totals, dtype=float).tolist(),
                         "mito_pct": np.where(available, mito, 0.).astype(float).tolist(),
                         "mito_pct_available": np.asarray(available, dtype=bool).tolist()})


def _thresholds(detected, mito, multiplier: int) -> dict[str, Any]:
    logs = np.log1p(detected)
    middle = np.median(logs)
    lower = middle - multiplier * 1.4826 * np.median(np.abs(logs - middle))
    maximum = 100.0
    if mito is not None:
        middle = np.median(mito)
        maximum = min(100., max(5., float(middle + multiplier * 1.4826 * np.median(np.abs(mito - middle)))))
    return {"min_genes": max(1, int(np.floor(np.expm1(lower)))), "max_mito_pct": maximum}


def _retained_identity(ids: list[str], retained) -> str:
    return _fingerprint([cell for cell, keep in zip(ids, retained) if keep])


def freeze_qc_sensitivity_plan(*, payload, observations: Mapping[str, Any],
                               source_sha256: str, baseline_config: Mapping[str, Any]) -> dict[str, Any]:
    """Freeze the primary rule and exactly two perturbations from raw counts.

    The source digest must be independently measured by the importer. No
    reference, expected annotation, desired retention or branch result is an
    accepted input. Observational labels cannot influence the thresholds.
    """
    from eacbp.application.autonomous import observe_input

    if not isinstance(source_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", source_sha256):
        raise ValueError("source_sha256 must be a lowercase SHA-256 digest")
    question = baseline_config.get("autonomous_analysis", {}).get("question")
    observed = observe_input(payload, question=question)
    if observed != observations:
        raise ValueError("Input observations changed before QC sensitivity planning")
    settings = baseline_config.get("autonomous_analysis", {})
    if settings.get("source_sha256") != source_sha256:
        raise ValueError("Baseline config and raw source identity disagree")
    source = observations["counts_source"]
    matrix = payload.layers["counts"] if source == "layers.counts" else payload.raw["X"] if source == "raw.X" else payload.X
    var = payload.raw["var"] if source == "raw.X" else payload.var
    gene_source = observations["gene_name_source"]
    genes = var.index.astype(str).str.strip().to_numpy() if gene_source == "var.index" else var[gene_source[4:]].astype(str).str.strip().to_numpy()
    mt = var["mt"].to_numpy(dtype=bool) if "mt" in var else np.char.startswith(np.char.lower(genes.astype(str)), "mt-")
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    mito = None
    if mt.any():
        mito = np.divide(np.asarray(matrix[:, mt].sum(axis=1)).ravel() * 100., totals,
                         out=np.zeros(len(totals)), where=totals > 0)
    ids = payload.obs.index.astype(str).tolist()
    if len(set(ids)) != len(ids):
        raise ValueError("String-normalized cell IDs must remain unique")
    main = {key: observations["qc"][key] for key in _PARAMETERS}
    configured = baseline_config.get("capability_parameters", {}).get("qc", {})
    if {key: configured.get(key) for key in _PARAMETERS} != main or main != _thresholds(detected, mito, 3):
        raise ValueError("Primary QC must remain the observed unmodified 3 scaled MAD policy")

    def spec(branch_id, multiplier):
        parameters = _thresholds(detected, mito, multiplier)
        branch_config = _config_scope(baseline_config)
        branch_config["capability_parameters"]["qc"].update(parameters)
        retained = detected >= parameters["min_genes"]
        if mito is not None:
            retained &= mito <= parameters["max_mito_pct"]
        n_retained = int(retained.sum())
        status = "planned"
        reason = "Only the two frozen QC thresholds change; downstream methods and scientific scope remain fixed"
        if n_retained < 3:
            status, reason = "not_estimable", "Fewer than three predicted retained cells; no threshold rescue or replacement estimand is allowed"
        elif branch_id != "primary" and parameters == main:
            status, reason = "duplicate_primary", "The observed distribution produces exactly the primary thresholds"
        return {"branch_id": branch_id, "scaled_mad_multiplier": multiplier,
                "parameters": parameters, "planned_status": status, "reason": reason,
                "config_sha256": _fingerprint(branch_config),
                "predicted_retained_cells": n_retained,
                "predicted_retained_ids_sha256": _retained_identity(ids, retained)}

    group_column = observations.get("obs_column_map", {}).get("cell_type") if observations.get("use_observed_cell_types") else None
    groups = _normalized_groups(payload.obs[group_column], ids) if group_column else None
    available = (totals > 0) if mito is not None else np.zeros(len(ids), dtype=bool)
    plan = {"schema_version": 1, "policy_version": QC_SENSITIVITY_POLICY_VERSION,
            "source_sha256": source_sha256, "counts_source": source,
            "baseline_config_sha256": qc_config_sha256(baseline_config),
            "cell_ids": ids, "n_features": int(matrix.shape[1]),
            "raw_qc_metrics_sha256": _metric_identity(ids, detected, totals, mito if mito is not None else np.zeros(len(ids)), available),
            "mitochondrial_metric_source": "var.mt" if "mt" in var and mt.any() else "gene_symbol_prefix" if mt.any() else "unavailable",
            "n_mitochondrial_genes": int(mt.sum()),
            "primary": spec("primary", 3), "branches": [spec("stricter_2mad", 2), spec("wider_4mad", 4)],
            "group_retention": {"authorized": group_column is not None, "source_column": group_column,
                                "labels_sha256": _fingerprint(groups) if groups is not None else None},
            "allowed_parameter_changes": ["capability_parameters.qc.min_genes", "capability_parameters.qc.max_mito_pct"],
            "interpretation": "Same-input sensitivity analysis, not independent validation. Low RNA alone does not establish cell damage. No branch replaces the primary or is selected using reference agreement."}
    plan["plan_sha256"] = _fingerprint(plan)
    plan["plan_id"] = "qc_sensitivity_" + plan["plan_sha256"][:16]
    return plan


def _validate_plan(plan: Mapping[str, Any]) -> None:
    clean = {key: value for key, value in plan.items() if key not in ("plan_id", "plan_sha256")}
    digest = _fingerprint(clean)
    if plan.get("policy_version") != QC_SENSITIVITY_POLICY_VERSION or plan.get("plan_sha256") != digest or plan.get("plan_id") != "qc_sensitivity_" + digest[:16]:
        raise ValueError("QC sensitivity plan identity or policy changed")


def qc_branch_config(plan: Mapping[str, Any], branch_id: str,
                     baseline_config: Mapping[str, Any]) -> dict[str, Any]:
    """Derive a branch from the frozen main configuration without retuning."""
    _validate_plan(plan)
    if qc_config_sha256(baseline_config) != plan["baseline_config_sha256"]:
        raise ValueError("Baseline analysis scope/config changed after QC sensitivity was frozen")
    branch = next((item for item in plan["branches"] if item["branch_id"] == branch_id), None)
    if branch is None:
        raise ValueError("Unknown frozen QC sensitivity branch")
    if branch["planned_status"] != "planned":
        raise ValueError("QC sensitivity branch is duplicate or not estimable; it must not be executed")
    result = deepcopy(dict(baseline_config))
    result["capability_parameters"]["qc"].update(branch["parameters"])
    result["autonomous_analysis"]["qc_sensitivity_branch"] = {
        "plan_id": plan["plan_id"], "plan_sha256": plan["plan_sha256"],
        "branch_id": branch_id, "source_sha256": plan["source_sha256"],
        "baseline_config_sha256": plan["baseline_config_sha256"],
        "parameters": deepcopy(branch["parameters"]),
        "interpretation": "Prespecified same-source QC perturbation; primary policy and estimand are unchanged"}
    return result


def _disposition(plan, spec, receipt):
    if receipt.get("source_sha256") != plan["source_sha256"]:
        raise ValueError("QC comparison uses a different raw source")
    if receipt.get("parameters") != spec["parameters"]:
        raise ValueError("Executed QC thresholds differ from the frozen branch")
    if receipt.get("config_sha256") != spec["config_sha256"]:
        raise ValueError("Executed analysis scope/config differs from the frozen QC branch")
    if receipt.get("qc_admitted") is not True:
        return None
    frame = pd.DataFrame(receipt.get("cell_disposition", {})).copy()
    required = {"ids", "n_genes", "total_counts", "mito_pct", "mito_pct_available", "retained", "reason", "low_gene_flag", "high_mito_flag"}
    if not required.issubset(frame):
        raise ValueError("Admitted QC comparison requires every input cell's complete disposition")
    frame["ids"] = frame["ids"].astype(str)
    if frame["ids"].duplicated().any() or set(frame["ids"]) != set(plan["cell_ids"]):
        raise ValueError("QC disposition changed the raw cell universe")
    frame = frame.set_index("ids").loc[plan["cell_ids"]]
    for name in ("retained", "mito_pct_available", "low_gene_flag", "high_mito_flag"):
        if frame[name].isna().any() or not frame[name].map(lambda value: isinstance(value, (bool, np.bool_))).all():
            raise ValueError("QC disposition flags must be complete booleans")
    detected = pd.to_numeric(frame["n_genes"], errors="raise").to_numpy(dtype=float)
    totals = pd.to_numeric(frame["total_counts"], errors="raise").to_numpy(dtype=float)
    mito = pd.to_numeric(frame["mito_pct"], errors="raise").to_numpy(dtype=float)
    available = frame["mito_pct_available"].to_numpy(dtype=bool)
    if not np.isfinite(detected).all() or np.any(detected < 0) or np.any(detected != np.floor(detected)) or not np.isfinite(totals).all() or np.any(totals < 0):
        raise ValueError("QC raw count metrics must be finite and nonnegative")
    expected_available = totals > 0 if plan["n_mitochondrial_genes"] else np.zeros(len(totals), dtype=bool)
    if not np.array_equal(available, expected_available) or not np.isfinite(mito[available]).all() or np.any((mito[available] < 0) | (mito[available] > 100)):
        raise ValueError("QC mitochondrial metric availability differs from the frozen raw rule")
    if not np.isnan(mito[~available]).all():
        raise ValueError("Unavailable QC mitochondrial metrics must remain explicitly missing")
    if _metric_identity(plan["cell_ids"], detected, totals, mito, available) != plan["raw_qc_metrics_sha256"]:
        raise ValueError("QC raw count/mitochondrial metrics changed between paired runs")
    low = detected < spec["parameters"]["min_genes"]
    high = np.where(available, mito, 0.) > spec["parameters"]["max_mito_pct"]
    retained = ~(low | high)
    reasons = np.full(len(frame), "retained", dtype=object)
    reasons[low], reasons[high], reasons[low & high] = _REASONS[1], _REASONS[2], _REASONS[3]
    if not np.array_equal(frame["low_gene_flag"], low) or not np.array_equal(frame["high_mito_flag"], high) or not np.array_equal(frame["retained"], retained) or not np.array_equal(frame["reason"].astype(str), reasons):
        raise ValueError("QC dispositions do not implement the frozen parameters")
    if _retained_identity(plan["cell_ids"], retained) != spec["predicted_retained_ids_sha256"]:
        raise ValueError("QC retained cell identity disagrees with the pre-execution plan")
    return frame


def _cluster_comparison(main, branch, main_frame, branch_frame):
    if main is None or branch is None:
        return {"status": "not_available", "reason": "Both admitted clustering outputs are required"}
    for labels, frame in ((main, main_frame), (branch, branch_frame)):
        if not isinstance(labels, pd.Series) or not labels.index.is_unique or not labels.index.astype(str).is_unique or labels.isna().any() or set(map(str, labels.index)) != set(frame.index[frame["retained"]]):
            raise ValueError("Cluster labels must cover exactly the admitted retained cells")
    ids = main_frame.index[main_frame["retained"] & branch_frame["retained"]]
    if len(ids) < 2:
        return {"status": "not_estimable", "n_shared_cells": len(ids), "reason": "Fewer than two shared retained cells"}
    first, second = main.copy(), branch.copy()
    first.index, second.index = first.index.astype(str), second.index.astype(str)
    a, b = first.loc[ids].astype(str).to_numpy(), second.loc[ids].astype(str).to_numpy()
    counts = pd.crosstab(a, b).to_numpy()
    pairs = lambda values: float(np.sum(values * (values - 1) / 2))
    together = pairs(counts)
    row, column = pairs(counts.sum(axis=1)), pairs(counts.sum(axis=0))
    expected = row * column / (len(ids) * (len(ids) - 1) / 2)
    maximum = (row + column) / 2
    ari = 1.0 if maximum == expected else (together - expected) / (maximum - expected)
    return {"status": "computed", "n_shared_cells": len(ids), "primary_shared_clusters": len(set(a)),
            "branch_shared_clusters": len(set(b)), "shared_cell_ari": float(ari),
            "interpretation": "Agreement of reclustering the same shared cells; not independent recovery of true cell identities"}


def _deg_rows(records):
    result = {}
    for record in records or []:
        if record.get("admitted") is not True:
            continue
        frame = record["table"]
        required = {"gene", "contrast_tested_level", "contrast_reference_level", "log2_fold_change", "fdr_q_value", "status", "alpha"}
        if not isinstance(frame, pd.DataFrame) or not required.issubset(frame):
            raise ValueError("DEG sensitivity requires native gene, oriented contrast, status, alpha and effect/FDR columns")
        for row in frame.to_dict("records"):
            if any(pd.isna(row[key]) or not str(row[key]).strip() for key in ("gene", "contrast_tested_level", "contrast_reference_level")):
                raise ValueError("DEG sensitivity requires an explicit observed categorical contrast")
            key = (str(record.get("target_cell_type") or "all_observations"), str(row["contrast_tested_level"]), str(row["contrast_reference_level"]), str(row["gene"]))
            if key in result:
                raise ValueError("Duplicate native DEG estimand/gene in QC sensitivity comparison")
            alpha = float(row["alpha"])
            if not math.isfinite(alpha) or not 0 < alpha < 1:
                raise ValueError("DEG alpha must be finite and between zero and one")
            if row["status"] not in ("estimated", "not_estimated", "estimated_fdr_missing"):
                raise ValueError("Unknown native DEG estimation status")
            for field in ("log2_fold_change", "fdr_q_value"):
                if pd.notna(row[field]) and not math.isfinite(float(row[field])):
                    raise ValueError("DEG sensitivity cannot hide an infinite statistical result")
            if pd.notna(row["fdr_q_value"]) and not 0 <= float(row["fdr_q_value"]) <= 1:
                raise ValueError("DEG FDR must be between zero and one when estimated")
            result[key] = row
    return result


def compare_deg_tables(baseline_tables, branch_tables) -> dict[str, Any]:
    """Compare the full gene universe, preserving missing/unestimated results."""
    primary, branch = _deg_rows(baseline_tables), _deg_rows(branch_tables)
    rows = []
    direction_counts = {"same_nonzero_direction": 0, "opposite_direction": 0,
                        "zero_in_either": 0, "not_estimable": 0, "missing_in_primary": 0, "missing_in_branch": 0}
    fdr_changes, paired_estimated = 0, 0
    number = lambda row, key: float(row[key]) if row is not None and pd.notna(row.get(key)) and math.isfinite(float(row[key])) else None
    for key in sorted(primary.keys() | branch.keys()):
        a, b = primary.get(key), branch.get(key)
        af, bf = number(a, "log2_fold_change"), number(b, "log2_fold_change")
        aq, bq = number(a, "fdr_q_value"), number(b, "fdr_q_value")
        status = "missing_in_primary" if a is None else "missing_in_branch" if b is None else "not_estimable"
        if a is not None and b is not None and af is not None and bf is not None and a["status"] != "not_estimated" and b["status"] != "not_estimated":
            paired_estimated += 1
            status = "zero_in_either" if af == 0 or bf == 0 else "same_nonzero_direction" if np.sign(af) == np.sign(bf) else "opposite_direction"
        direction_counts[status] += 1
        if a is not None and b is not None and a["alpha"] != b["alpha"]:
            raise ValueError("DEG significance alpha changed across QC branches")
        significant_a = aq < float(a["alpha"]) if aq is not None else None
        significant_b = bq < float(b["alpha"]) if bq is not None else None
        if significant_a is not None and significant_b is not None and significant_a != significant_b:
            fdr_changes += 1
        rows.append({"target_cell_type": key[0], "tested": key[1], "reference": key[2], "gene": key[3],
                     "direction_status": status, "primary_row_status": a["status"] if a is not None else "missing",
                     "branch_row_status": b["status"] if b is not None else "missing", "primary_log2_fold_change": af,
                     "branch_log2_fold_change": bf, "primary_fdr": aq, "branch_fdr": bq,
                     "primary_significant": significant_a, "branch_significant": significant_b,
                     "primary_design": {field: a[field] for field in ("design", "paired", "donor_col", "n_pseudobulk_samples", "n_donors_condition_a", "n_donors_condition_b") if a is not None and field in a},
                     "branch_design": {field: b[field] for field in ("design", "paired", "donor_col", "n_pseudobulk_samples", "n_donors_condition_a", "n_donors_condition_b") if b is not None and field in b},
                     "effect_change": bf - af if af is not None and bf is not None else None})
    return {"status": "computed" if rows else "not_available", "n_primary_rows": len(primary),
            "n_branch_rows": len(branch), "n_paired_estimated": paired_estimated,
            "direction_counts": direction_counts, "n_fdr_status_changes": fdr_changes, "rows": rows,
            "interpretation": "Paired same-input sensitivity with full missing-result accounting; no independent-validation or universal-stability claim"}


def compare_qc_sensitivity_runs(plan: Mapping[str, Any], *, baseline: Mapping[str, Any],
                                branches: Mapping[str, Mapping[str, Any]],
                                observed_groups: pd.Series | None = None) -> dict[str, Any]:
    """Build a detached derived report; never modify a native run or its status."""
    _validate_plan(plan)
    allowed = {item["branch_id"] for item in plan["branches"]}
    if not set(branches).issubset(allowed):
        raise ValueError("Unplanned QC branch cannot enter the sensitivity report")
    group_scope = plan["group_retention"]
    if observed_groups is not None and not group_scope["authorized"]:
        raise ValueError("Observed group retention was not authorized by the frozen question")
    groups = _normalized_groups(observed_groups, plan["cell_ids"]) if observed_groups is not None else None
    if groups is not None and _fingerprint(groups) != group_scope["labels_sha256"]:
        raise ValueError("Observed group labels changed after QC planning")
    main_frame = _disposition(plan, plan["primary"], baseline)
    comparisons = []
    for spec in plan["branches"]:
        branch_id = spec["branch_id"]
        receipt = branches.get(branch_id)
        item = {"branch_id": branch_id, "parameters": deepcopy(spec["parameters"]),
                "planned_status": spec["planned_status"], "run_status": receipt.get("run_status") if receipt else "not_run"}
        if spec["planned_status"] != "planned":
            if receipt is not None:
                raise ValueError("Duplicate/not-estimable branch must not have been executed")
            item.update(status="not_run", reason=spec["reason"])
        elif receipt is None:
            item.update(status="missing", reason="A planned branch has no native run receipt")
        else:
            frame = _disposition(plan, spec, receipt)
            if frame is None or main_frame is None:
                item.update(status="not_estimable", reason="Primary and branch QC must both be independently admitted")
            else:
                first, second = main_frame["retained"], frame["retained"]
                transitions = pd.crosstab(main_frame["reason"], frame["reason"])
                changed = first != second
                details = []
                for cell, a, b in zip(plan["cell_ids"], main_frame.to_dict("records"), frame.to_dict("records")):
                    details.append({"cell_id": cell, "n_genes": int(a["n_genes"]), "total_counts": float(a["total_counts"]),
                                    "mito_pct": float(a["mito_pct"]) if a["mito_pct_available"] else None,
                                    "mito_pct_available": bool(a["mito_pct_available"]), "primary_retained": bool(a["retained"]),
                                    "branch_retained": bool(b["retained"]), "primary_reason": str(a["reason"]), "branch_reason": str(b["reason"])})
                retention = []
                if groups is not None:
                    labels = pd.Series(groups, index=main_frame.index, dtype="string")
                    for group in sorted(set(value for value in groups if value is not None)):
                        mask = labels.eq(group).fillna(False).to_numpy(dtype=bool)
                        original = int(mask.sum())
                        retained_a, retained_b = int(first.to_numpy()[mask].sum()), int(second.to_numpy()[mask].sum())
                        retention.append({"observed_group": group, "original_cells": original,
                                          "primary_retained": retained_a, "branch_retained": retained_b,
                                          "primary_retention_fraction": retained_a / original, "branch_retention_fraction": retained_b / original})
                item.update(status="compared", primary_retained=int(first.sum()), branch_retained=int(second.sum()),
                            n_added=int((~first & second).sum()), n_removed=int((first & ~second).sum()), n_changed=int(changed.sum()),
                            disposition_transitions=[{"primary_reason": a, "branch_reason": b, "n_cells": int(transitions.loc[a, b])}
                                                     for a in transitions.index for b in transitions.columns if transitions.loc[a, b]],
                            cell_dispositions=details, observed_group_retention=retention,
                            group_retention_status="computed" if groups is not None else "not_supplied" if group_scope["authorized"] else "not_authorized",
                            clustering=_cluster_comparison(baseline.get("clusters"), receipt.get("clusters"), main_frame, frame),
                            deg=compare_deg_tables(baseline.get("deg_tables"), receipt.get("deg_tables")),
                            downstream_summaries={"primary": deepcopy(baseline.get("summaries", {})), "branch": deepcopy(receipt.get("summaries", {}))})
        comparisons.append(item)
    result = {"schema_version": 1, "plan_id": plan["plan_id"], "plan_sha256": plan["plan_sha256"],
              "source_sha256": plan["source_sha256"], "primary_parameters": deepcopy(plan["primary"]["parameters"]),
              "primary_run_status": baseline.get("run_status"),
              "status": "completed" if main_frame is not None and baseline.get("run_status") in ("success", "completed") and all(
                  item["status"] == "not_run" or item["status"] == "compared" and item["run_status"] in ("success", "completed") for item in comparisons) else "partial",
              "completion_scope": "Prespecified QC perturbation execution/comparison; this status does not establish scientific robustness or full question coverage",
              "branches": comparisons, "interpretation": plan["interpretation"]}
    # Native unestimated statistics become explicit nulls above. The report
    # itself must never acquire NaN/Inf or unsupported hidden objects.
    _fingerprint(result)
    return result
