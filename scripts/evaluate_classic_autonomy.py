"""External, post-run acceptance for the classic autonomous exercises.

Expected labels, genes and direction edges live only in the evaluation directory.
This program never calls planning/execution or writes run configuration. A passed
native run and a passed external scientific benchmark are distinct outcomes.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import sys

import numpy as np
import pandas as pd
from scipy import sparse

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def check(name, passed, **evidence):
    return {"name": name, "passed": bool(passed), **evidence}


def _summary_value(value):
    """Compare JSON summaries while permitting equivalent timezone renderings."""
    if isinstance(value, dict):
        return {key: _summary_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_summary_value(item) for item in value]
    if isinstance(value, str) and re.match(r"^\d{4}-\d{2}-\d{2}T", value):
        try:
            date = datetime.fromisoformat(value.replace("Z", "+00:00"))
            if date.tzinfo is not None:
                return date.astimezone(timezone.utc).isoformat()
        except ValueError:
            pass
    return value


def _completed_native_snapshot(summary, snapshot):
    if _summary_value(summary) != _summary_value(snapshot.summary):
        raise ValueError("Native summary does not match its validated snapshot")
    if snapshot.summary.get("status") != "success" or snapshot.summary.get("analysis_request_status") != "complete":
        raise ValueError("Native snapshot must complete before answer-bearing references are opened")
    contexts = list(snapshot.audit_contexts.values())
    admitted = {context.task_id for context in contexts if context.receipt_present and context.status == "passed"
                and context.overall_passed and not context.stop_rule_triggered}
    if not snapshot.task_results or any(task.status.value != "success" or task.task_id not in admitted
                                        for task in snapshot.task_results):
        raise ValueError("Every native task needs successful execution and persisted audit admission before references")
    if not contexts or any(not context.receipt_present or context.status != "passed" or not context.overall_passed
                           or context.stop_rule_triggered for context in contexts):
        raise ValueError("Native audit admission must complete before references")


def _container_path(value):
    if not isinstance(value, str):
        raise ValueError("Unsupported container path in execution receipt")
    path = PurePosixPath(value)
    if not path.is_absolute() or ".." in path.parts or "\\" in value:
        raise ValueError("Unsupported container path in execution receipt")
    return path


def _host_path(value):
    # Resolve on the actual evaluation host, including Windows junctions and
    # symlinks. Foreign-host path syntax is deliberately unsupported.
    path = Path(value)
    if not path.is_absolute():
        raise ValueError("Mount source must be an absolute path on the evaluation host")
    return path.resolve(strict=True)


def _mount(value, *, volume=False):
    if volume:
        match = re.fullmatch(r"(.+):(/[^:]+)(?::(ro|rw))?", value)
        if not match:
            raise ValueError("Unsupported --volume shape in execution receipt")
        source, target, mode = match.groups()
        return {"source": _host_path(source), "target": _container_path(target), "readonly": mode == "ro"}
    aliases = {"src": "source", "dst": "target", "destination": "target", "ro": "readonly"}
    fields = {}
    for item in value.split(","):
        key, separator, field = item.partition("=")
        key = aliases.get(key, key)
        if key not in {"type", "source", "target", "readonly"} or key in fields:
            raise ValueError("Unknown or repeated --mount field in execution receipt")
        if not separator:
            if key != "readonly":
                raise ValueError("Incomplete --mount field in execution receipt")
            field = "true"
        fields[key] = field
    if fields.get("type") != "bind" or not {"source", "target"}.issubset(fields):
        raise ValueError("Only explicit bind mounts are supported by this evaluator")
    if fields.get("readonly", "false") not in {"true", "false", "1", "0"}:
        raise ValueError("Unknown mount readonly mode")
    return {"source": _host_path(fields["source"]), "target": _container_path(fields["target"]),
            "readonly": fields.get("readonly", "false") in {"true", "1"}}


def _execution_boundary(project, run_root, reference_root, snapshot, expected_image):
    """Validate local runner records, not a signed remote execution attestation."""
    if not isinstance(expected_image, str) or not re.fullmatch(r"sha256:[0-9a-f]{64}", expected_image):
        raise ValueError("An independently checked expected-image sha256 digest is required")
    receipt_path = run_root.parent / f"{run_root.name}-execution-command.json"
    if not receipt_path.is_file():
        raise ValueError("Missing execution receipt for this exact native run")
    receipt_hash = sha256(receipt_path)
    execution = json.loads(receipt_path.read_text(encoding="utf-8"))
    command = execution.get("command")
    if execution.get("project") != project or execution.get("answer_reference_mounted") is not False:
        raise ValueError("Execution receipt project/reference declaration does not match")
    if not isinstance(command, list) or any(not isinstance(item, str) for item in command) or len(command) < 4:
        raise ValueError("Unsupported execution command receipt")
    if Path(command[0]).name not in {"docker", "docker.exe"} or command[1] != "run":
        raise ValueError("Only the isolated Docker runner is supported")
    options, mounts = {}, []
    value_options = {"--network", "--user", "--entrypoint", "--env", "-e", "--mount", "--volume", "-v"}
    i = 2
    while i < len(command) and command[i].startswith("-"):
        token = command[i]
        name, equal, value = token.partition("=")
        if name in {"--rm", "--read-only"} and not equal:
            if name in options:
                raise ValueError("Repeated Docker execution option")
            options[name] = True
            i += 1
            continue
        if name not in value_options:
            raise ValueError(f"Unsupported Docker option: {name}")
        if not equal:
            i += 1
            if i >= len(command):
                raise ValueError("Missing Docker option value")
            value = command[i]
        if name in {"--mount", "--volume", "-v"}:
            mounts.append(_mount(value, volume=name != "--mount"))
        else:
            name = "--env" if name == "-e" else name
            if name in options:
                raise ValueError("Repeated Docker execution option")
            options[name] = value
        i += 1
    if i >= len(command) or command[i] != expected_image:
        raise ValueError("Execution receipt image differs from the independently checked frozen image")
    image = command[i]
    if options != {"--rm": True, "--read-only": True, "--network": "none", "--user": "10001:10001",
                   "--entrypoint": "/opt/venv/bin/eacbp", "--env": "PYTHONPATH=/linux-deps:/deps"}:
        raise ValueError("Execution receipt lacks the supported isolated runner settings")
    arguments = command[i + 1:]
    if not arguments or arguments[0] != "analyze" or len(arguments[1:]) % 2:
        raise ValueError("Unsupported native analyze invocation")
    flags = {}
    allowed = {"--question", "--data", "--species", "--tissue", "--run-dir", "--threads", "--max-runtime-hours"}
    for key, value in zip(arguments[1::2], arguments[2::2]):
        if key not in allowed or key in flags:
            raise ValueError("Unknown or repeated public analyze option")
        flags[key] = value
    required = allowed - {"--max-runtime-hours"}
    if not required.issubset(flags):
        raise ValueError("Incomplete question-only analyze invocation")
    autonomy = snapshot.config.get("autonomous_analysis", {})
    design = snapshot.manifest.biological_design
    source = json.loads((run_root / "source.json").read_text(encoding="utf-8"))
    if (flags["--question"] != autonomy.get("question") or flags["--species"] != design.species
            or flags["--tissue"] != design.tissue or flags["--data"] != source.get("path")
            or flags["--threads"] != str(snapshot.config.get("threads"))):
        raise ValueError("Execution receipt question/source/context differs from the native run")
    if "--max-runtime-hours" in flags and float(flags["--max-runtime-hours"]) != snapshot.manifest.constraints.max_runtime_hours:
        raise ValueError("Execution receipt runtime constraint differs from the native run")
    reference = reference_root.resolve(strict=True)
    by_target = {}
    for mount in mounts:
        host = mount["source"]
        if host == reference or reference in host.parents or host in reference.parents:
            raise ValueError("Execution mount exposes reference answers through a path or alias")
        if mount["target"] in by_target:
            raise ValueError("Repeated mount target")
        by_target[mount["target"]] = mount
    root = run_root.parent.parent.resolve(strict=True)
    data_target = _container_path(flags["--data"])
    expected_targets = {data_target, PurePosixPath("/results"), PurePosixPath("/scratch"),
                        PurePosixPath("/deps"), PurePosixPath("/linux-deps")}
    if set(by_target) != expected_targets:
        raise ValueError("Unknown, overlapping or missing runner mount targets")
    input_mount = by_target[data_target]
    if not input_mount["readonly"] or not input_mount["source"].is_file() or input_mount["source"].parent != (root / "datasets").resolve(strict=True):
        raise ValueError("The native input must be the single read-only curriculum dataset file")
    for target, relative, readonly in (("/results", "runs", False), ("/scratch", "cache", False),
                                       ("/deps", "dependencies/site-packages", True),
                                       ("/linux-deps", "dependencies/linux-site-packages", True)):
        mount = by_target[PurePosixPath(target)]
        if mount["source"] != (root / relative).resolve(strict=True) or not mount["source"].is_dir() or mount["readonly"] != readonly:
            raise ValueError(f"Unknown runner mount source or mode for {target}")
    container_run = _container_path(flags["--run-dir"])
    if container_run.parent != PurePosixPath("/results") or container_run.name != run_root.name:
        raise ValueError("Execution receipt belongs to a different native run directory")
    digest = sha256(input_mount["source"])
    if any(str(source.get(key, "")).removeprefix("sha256:") != digest
           for key in ("sha256", "sha256_before_import", "sha256_after_import")) or digest != str(autonomy.get("source_sha256", "")).removeprefix("sha256:"):
        raise ValueError("Execution receipt source file/hash differs from the native import")
    return {"path": receipt_path, "sha256": receipt_hash, "expected_image": expected_image, "receipt_image": image,
            "public_flags": sorted(flags), "source_path": str(input_mount["source"]), "source_sha256": digest,
            "attestation_scope": "Local runner receipt consistency; not a signed remote execution proof"}


def counts_equal(actual, expected):
    if actual.shape != expected.shape:
        return False
    if sparse.issparse(actual) or sparse.issparse(expected):
        return bool((sparse.csr_matrix(actual) != sparse.csr_matrix(expected)).nnz == 0)
    return bool(np.array_equal(actual, expected))


def label_series(data, key):
    obs = data.obs
    if key not in obs or not obs.index.is_unique:
        raise ValueError(f"Unique cell IDs and obs[{key!r}] are required")
    values = obs[key]
    if values.isna().any() or values.astype(str).str.strip().eq("").any():
        raise ValueError(f"Missing {key} labels cannot be scored as a successful population")
    return pd.Series(values.astype(str).to_numpy(), index=obs.index.astype(str))


def population_retention(truth, retained_ids, thresholds):
    if not truth.index.is_unique or not pd.Index(retained_ids).is_unique:
        raise ValueError("Repeated cell IDs are not independent reference coverage")
    retained = truth.index.isin(pd.Index(retained_ids).astype(str))
    per_class = {}
    for label in sorted(truth.unique()):
        mask = truth.to_numpy() == label
        n, kept = int(mask.sum()), int((mask & retained).sum())
        per_class[label] = {"source_cells": n, "retained_cells": kept,
                            "excluded_cells": n - kept, "retention": kept / n}
    fraction = float(retained.mean()) if len(truth) else 0.0
    return {"overall_retention": fraction, "per_class": per_class,
            "checks": [check("overall_reference_retention", fraction >= thresholds["overall_reference_retention_min"]),
                       check("every_reference_class_retained", bool(per_class) and all(
                           row["retention"] >= thresholds["per_class_retention_min"] for row in per_class.values()))]}


def clustering_comparison(truth, prediction, thresholds):
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    shared = truth.index.intersection(prediction.index, sort=False)
    retention = population_retention(truth, prediction.index, thresholds)
    if not len(shared):
        raise ValueError("No exact reference barcodes match the newly computed clusters")
    ari = float(adjusted_rand_score(truth.loc[shared], prediction.loc[shared]))
    nmi = float(normalized_mutual_info_score(truth.loc[shared], prediction.loc[shared]))
    table = pd.crosstab(truth.loc[shared], prediction.loc[shared])
    return {**retention, "n_shared_cells": len(shared), "adjusted_rand_index": ari,
            "normalized_mutual_information": nmi,
            "crosstab": table.to_dict(), "checks": retention["checks"] + [
                check("clustering_ari", ari >= thresholds["adjusted_rand_index_min"]),
                check("clustering_nmi", nmi >= thresholds["normalized_mutual_information_min"])]}


def marker_recovery(truth, prediction, markers, expectations, thresholds):
    required = {"cluster", "gene", "rank", "mean_log1p_difference"}
    if not required.issubset(markers):
        raise ValueError(f"Native descriptive marker table is missing {sorted(required - set(markers))}")
    shared = truth.index.intersection(prediction.index, sort=False)
    contingency = pd.crosstab(truth.loc[shared], prediction.loc[shared])
    classes = {}
    for label, genes in expectations.items():
        cluster = str(contingency.loc[label].idxmax()) if label in contingency.index else None
        purity = float(contingency.loc[label, cluster] / contingency[cluster].sum()) if cluster else 0.0
        rows = markers.loc[markers["cluster"].astype(str).eq(cluster)]
        hits = sorted(set(rows.loc[(rows["rank"] <= thresholds["marker_rank_max"])
                                   & (rows["mean_log1p_difference"] > 0), "gene"].astype(str)).intersection(genes))
        classes[label] = {"evaluation_only_dominant_cluster": cluster, "cluster_reference_purity": purity,
                          "recovered_expected_markers": hits,
                          "supported": bool(hits and purity >= thresholds["marker_reference_cluster_purity_min"])}
    supported = sum(row["supported"] for row in classes.values())
    return {"classes": classes, "n_supported_classes": supported,
            "checks": [check("independent_tutorial_marker_recovery", supported >= thresholds["minimum_supported_marker_classes"])]}


def velocity_direction(graph, labels, edges, thresholds):
    graph = sparse.csr_matrix(graph, dtype=np.float64, copy=True)
    if graph.shape != (len(labels), len(labels)) or not np.isfinite(graph.data).all() or (graph.data < 0).any():
        raise ValueError("Positive velocity graph must be finite, nonnegative and cell-aligned")
    graph.sum_duplicates()
    graph.eliminate_zeros()
    totals = np.asarray(graph.sum(axis=1)).ravel()
    normalized = sparse.diags(np.divide(1., totals, out=np.zeros_like(totals), where=totals > 0)) @ graph
    output = []
    values = labels.to_numpy()
    for source, target in edges:
        left, right = values == source, values == target
        if not left.any() or not right.any():
            output.append({"source": source, "target": target, "passed": False, "reason": "A reference state is absent"})
            continue
        forward = float(np.asarray(normalized[left][:, right].sum(axis=1)).mean())
        reverse = float(np.asarray(normalized[right][:, left].sum(axis=1)).mean())
        n_edges = int(graph[left][:, right].nnz + graph[right][:, left].nnz)
        index = (forward - reverse) / (forward + reverse) if forward + reverse > 0 else 0.0
        output.append({"source": source, "target": target, "forward_flow": forward, "reverse_flow": reverse,
                       "direction_index": index, "cross_state_edges": n_edges,
                       "passed": bool(index >= thresholds["direction_index_min"]
                                      and n_edges >= thresholds["minimum_cross_state_edges"])})
    return {"edges": output, "checks": [check("published_core_velocity_direction", bool(output) and all(row["passed"] for row in output))]}


def _source_subset(source, data, matrix_key="counts"):
    ids = pd.Index(data.obs.index.astype(str))
    genes = pd.Index(data.var.index.astype(str))
    source_ids = pd.Index(source.obs_names.astype(str))
    source_genes = pd.Index(source.var_names.astype(str))
    rows, columns = source_ids.get_indexer(ids), source_genes.get_indexer(genes)
    valid = ids.is_unique and genes.is_unique and (rows >= 0).all() and (columns >= 0).all()
    matrix = data.layers.get(matrix_key)
    expected = source.X[rows][:, columns] if valid else None
    return bool(valid and matrix is not None and counts_equal(matrix, expected))


def _task_artifacts(snapshot, registry, capability, *, kind=None, primary=False, load=True):
    output = []
    for task in snapshot.task_results:
        if task.capability != capability:
            continue
        for uri in task.output_artifacts:
            meta = snapshot.artifacts[uri]
            if primary and uri != task.output_artifacts[0]:
                continue
            if kind is None or meta.type.value == kind:
                output.append((task, meta, registry.get(uri)[1] if load else None))
    return output


def _one_scdata(snapshot, registry, capability):
    from eacbp.capabilities.sc_data import SCData
    artifacts = _task_artifacts(snapshot, registry, capability, kind="AnnData")
    if len(artifacts) != 1:
        raise ValueError(f"Expected exactly one shared {capability} AnnData; found {len(artifacts)}")
    task, meta, payload = artifacts[0]
    return task, meta, SCData.from_payload(payload)


def _target(meta):
    return meta.parameters.get("target_cell_type") or meta.parameters.get("inference_contract", {}).get("target_cell_type")


def _kang_targets(source, spec):
    obs = source.obs
    labels, donor, condition = (spec[key] for key in ("label_key", "donor_key", "condition_key"))
    requested = [spec["tested_condition"], spec["reference_condition"]]
    targets = {}
    for label in sorted(obs[labels].dropna().astype(str).unique()):
        group = obs.loc[obs[labels].astype(str).eq(label)]
        complete = not group[[donor, condition]].isna().any().any()
        sets = {level: set(group.loc[group[condition].astype(str).eq(level), donor].astype(str)) for level in requested}
        first, second = map(set, sets.values())
        overlap = bool(first & second)
        eligible = complete and min(len(first), len(second)) >= spec["eligibility"]["minimum_donors_per_condition"] and (not overlap or first == second)
        targets[label] = {"eligible": eligible, "paired": overlap,
                          "donors_by_condition": {key: sorted(value) for key, value in sets.items()}}
    return targets


def _independent_donor_effects(source, spec, label, genes, retained_ids):
    obs = source.obs
    source_ids = pd.Index(source.obs_names.astype(str))
    positions = pd.Index(source.var_names.astype(str)).get_indexer(genes)
    if (positions < 0).any():
        raise ValueError("The pinned Kang source lacks a prespecified response gene")
    donors = sorted(obs.loc[obs[spec["label_key"]].astype(str).eq(label), spec["donor_key"]].astype(str).unique())
    estimates = []
    for donor in donors:
        normalized = []
        for condition in (spec["tested_condition"], spec["reference_condition"]):
            mask = (obs[spec["label_key"]].astype(str).eq(label)
                    & obs[spec["donor_key"]].astype(str).eq(donor)
                    & obs[spec["condition_key"]].astype(str).eq(condition)).to_numpy()
            mask &= source_ids.isin(retained_ids)
            if not mask.any():
                break
            total = float(source.X[mask].sum())
            if total <= 0:
                break
            normalized.append(np.asarray(source.X[mask][:, positions].sum(axis=0)).ravel() * 10000. / total)
        if len(normalized) == 2:
            estimates.append(np.log2((normalized[0] + .5) / (normalized[1] + .5)))
    if len(estimates) != len(donors):
        return {"complete_donor_coverage": False, "n_complete_donors": len(estimates), "expected_donors": len(donors)}
    values = np.stack(estimates)
    return {"complete_donor_coverage": True, "n_complete_donors": len(donors),
            "normalization": "Independent retained-source donor sums; CP10K with fixed 0.5 offset, not a p value or replacement DE model",
            "genes": {gene: {"mean_log_ratio": float(values[:, i].mean()),
                              "positive_donor_fraction": float((values[:, i] > 0).mean())} for i, gene in enumerate(genes)}}


def leaveout_integrity(folds, summary, donors, genes, *, paired):
    required = {"left_out_donor", "gene", "removed_sample_count", "n_pseudobulk_samples"}
    summary_required = {"gene", "n_requested_fits", "n_successful_fits", "all_requested_fits_completed", "complete_fit_coverage"}
    if not required.issubset(folds) or not summary_required.issubset(summary):
        return {"checks": [check("complete_donor_leaveout_coverage", False), check("leaveout_summary_complete", False)]}
    expected_folds = set(donors) | {"__full_model__"}
    groups = folds.groupby("left_out_donor", sort=True)
    all_genes = set(genes)
    complete = (set(groups.groups) == expected_folds and not folds[["left_out_donor", "gene"]].duplicated().any()
                and all(set(group["gene"].astype(str)) == all_genes for _, group in groups))
    baseline = folds.loc[folds["left_out_donor"].eq("__full_model__"), "n_pseudobulk_samples"].dropna().unique()
    removed = 2 if paired else 1
    removed_coverage = len(baseline) == 1 and all(
        group["removed_sample_count"].eq(0 if name == "__full_model__" else removed).all()
        and group["n_pseudobulk_samples"].eq(baseline[0] if name == "__full_model__" else baseline[0] - removed).all()
        for name, group in groups)
    summary_complete = (not summary["gene"].duplicated().any() and set(summary["gene"].astype(str)) == all_genes
                        and summary["n_requested_fits"].eq(len(donors)).all()
                        and summary["n_successful_fits"].eq(len(donors)).all()
                        and summary["all_requested_fits_completed"].eq(True).all()
                        and summary["complete_fit_coverage"].eq(1.).all())
    return {"expected_folds": sorted(expected_folds), "observed_folds": sorted(map(str, groups.groups)),
            "checks": [check("complete_donor_leaveout_coverage", complete and removed_coverage),
                       check("leaveout_summary_complete", summary_complete)]}


def evaluate_kang(snapshot, registry, source, clustered, spec):
    threshold = spec["thresholds"]
    truth = label_series(source, spec["label_key"])
    retention = population_retention(truth, clustered.obs.index.astype(str), threshold)
    targets = _kang_targets(source, spec)
    expected = {name for name, row in targets.items() if row["eligible"]}
    de = _task_artifacts(snapshot, registry, "deg", kind="Table")
    # Read each large detailed LOO table only while evaluating its target.
    sensitivity = _task_artifacts(snapshot, registry, "donor_sensitivity", kind="Table", primary=True, load=False)
    de_map = {_target(meta): (task, meta, table) for task, meta, table in de}
    sensitivity_map = {_target(meta): (task, meta, table) for task, meta, table in sensitivity}
    checks = retention["checks"] + [check("all_eligible_observed_targets_analyzed", bool(expected)
               and set(de_map) == expected and len(de_map) == len(de)),
               check("all_eligible_donor_sensitivity_branches", set(sensitivity_map) == expected and len(sensitivity_map) == len(sensitivity))]
    populations = {}
    for label in sorted(expected):
        row = {"eligibility": targets[label], "checks": []}
        if label not in de_map:
            row["checks"].append(check("native_deg_present", False))
            populations[label] = row
            continue
        task, meta, table = de_map[label]
        required = {"gene", "log2_fold_change", "fdr_q_value", "condition_a", "condition_b", "statistical_unit", "paired", "status"}
        if not required.issubset(table) or table["gene"].duplicated().any():
            raise ValueError(f"Invalid full-gene donor DEG table for {label}")
        branch = registry.load_scdata(task.input_artifacts[0])
        genes = set((branch.var["gene_name"] if "gene_name" in branch.var else branch.var.index).astype(str))
        source_genes = set(source.var_names.astype(str))
        row["feature_disposition"] = {"source_feature_count": len(source_genes), "retained_native_feature_count": len(genes),
            "excluded_before_deg_count": len(source_genes - genes),
            "excluded_before_deg_features": sorted(source_genes - genes),
            "interpretation": "Native QC exclusions precede the declared DEG feature family; all retained native features must remain in the full table, including nonestimable rows"}
        row["checks"].append(check("full_declared_feature_coverage", set(table["gene"].astype(str)) == genes))
        branch_labels = source.obs[spec["label_key"]].reindex(branch.obs.index.astype(str))
        row["checks"].append(check("target_branch_raw_counts_and_membership", _source_subset(source, branch)
                                  and branch_labels.notna().all() and branch_labels.astype(str).eq(label).all()))
        row["checks"].append(check("observed_donor_estimand", table["condition_a"].eq(spec["tested_condition"]).all()
                                  and table["condition_b"].eq(spec["reference_condition"]).all()
                                  and table["statistical_unit"].eq("donor_pseudobulk").all()
                                  and table["paired"].eq(targets[label]["paired"]).all()
                                  and task.method_used == "pydeseq2_pseudobulk_v1"))
        q = pd.to_numeric(table["fdr_q_value"], errors="coerce")
        row["checks"].append(check("fdr_range_and_nonestimability", q.dropna().between(0, 1).all()
                                  and (table.loc[q.isna(), "status"].astype(str) != "estimated").all()))
        expected_genes = spec["expected_ifn_genes"]
        effects = table.set_index("gene").reindex(expected_genes)
        lfc = pd.to_numeric(effects["log2_fold_change"], errors="coerce")
        significant = pd.to_numeric(effects["fdr_q_value"], errors="coerce").le(threshold["fdr_max"]) & lfc.gt(0)
        n_positive = int((np.isfinite(lfc) & lfc.gt(0)).sum())
        independent = _independent_donor_effects(source, spec, label, expected_genes, clustered.obs.index.astype(str))
        donor_genes = independent.get("genes", {})
        independent_positive = sum(value["mean_log_ratio"] > 0 for value in donor_genes.values())
        concordant = sum(value["positive_donor_fraction"] >= threshold["donor_positive_fraction_min"] for value in donor_genes.values())
        row.update(n_positive_ifn_genes=n_positive, n_positive_fdr_ifn_genes=int(significant.sum()), independent_donor_effects=independent,
                   n_independent_positive_genes=independent_positive, n_donor_concordant_genes=concordant,
                   ifn_effects={gene: {"log2_fold_change": float(lfc.loc[gene]) if np.isfinite(lfc.loc[gene]) else None,
                                       "fdr_q_value": float(effects.loc[gene, "fdr_q_value"]) if pd.notna(effects.loc[gene, "fdr_q_value"]) else None} for gene in expected_genes})
        row["checks"].extend([check("ifn_effect_direction", n_positive >= threshold["ifn_positive_gene_min"]),
                               check("independent_donor_direction", independent["complete_donor_coverage"]
                                     and independent_positive >= threshold["ifn_independent_pseudobulk_positive_gene_min"])])
        if label in spec["tutorial_anchor_populations"]:
            row["checks"].extend([check("tutorial_anchor_fdr", int(significant.sum()) >= threshold["anchor_fdr_positive_gene_min"]),
                                   check("tutorial_anchor_donor_concordance", concordant >= threshold["anchor_donor_concordant_gene_min"])])
        if label in sensitivity_map:
            sensitivity_task, sensitivity_meta, _ = sensitivity_map[label]
            folds = registry.get(sensitivity_meta.uri)[1]
            donors = set(targets[label]["donors_by_condition"][spec["tested_condition"]])
            summaries = [uri for uri in sensitivity_task.output_artifacts[1:] if snapshot.artifacts[uri].type.value == "Table"]
            if len(summaries) == 1:
                summary_table = registry.get(summaries[0])[1]
                loo = leaveout_integrity(folds, summary_table, donors, genes, paired=targets[label]["paired"])
                row["leaveout_integrity"] = loo
                row["checks"] += loo["checks"]
            else:
                row["checks"].append(check("leaveout_summary_complete", False))
            del folds
        populations[label] = row
    checks.extend(item for row in populations.values() for item in row["checks"])
    return {**retention, "observed_target_eligibility": targets, "populations": populations, "checks": checks}


def evaluate(project, run_root, reference_root, *, expected_image=None):
    run_root, reference_root = Path(run_root), Path(reference_root)
    # Do not open answer-bearing references until native scientific execution
    # has published its final autonomous summary.
    summary = json.loads((run_root / "summary.json").read_text(encoding="utf-8"))
    if summary.get("status") != "success" or summary.get("analysis_request_status") != "complete":
        raise ValueError("Native autonomous analysis must complete before answer-bearing references are opened")
    import anndata as ad
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.evidence.snapshot import load_study_snapshot, _verify_registry_scope

    protected = [run_root / name for name in ("snapshot.json", "summary.json", "run_config.json", "manifest.json",
                 "config.json", "source.json", "autonomy.json", "input_preparation.json", "artifacts/.artifact_registry.json")]
    native_hashes_before = {str(path.relative_to(run_root)): sha256(path) for path in protected if path.is_file()}
    if not (run_root / "artifacts" / ".artifact_registry.json").is_file():
        raise ValueError("A completed native artifact registry is required")
    snapshot = load_study_snapshot(run_root / "snapshot.json")
    registry = ArtifactRegistry(str(run_root / "artifacts"))
    _verify_registry_scope(snapshot, registry)
    _completed_native_snapshot(summary, snapshot)
    execution_boundary = _execution_boundary(project, run_root, reference_root, snapshot, expected_image)
    # The native snapshot, actual artifact/audit hashes, summary and isolated
    # invocation have now been validated without opening answer-bearing files.
    protocol_path = reference_root / "blind_acceptance_protocol.json"
    reference_receipt_path = reference_root / "source_receipt.json"
    reference_receipt_hash = sha256(reference_receipt_path)
    receipt = json.loads(reference_receipt_path.read_text(encoding="utf-8-sig"))
    protocol = json.loads(protocol_path.read_text(encoding="utf-8-sig"))
    if sha256(protocol_path) != receipt.get("blind_protocol_sha256"):
        raise ValueError("The frozen external protocol hash does not match its source receipt")
    spec = protocol["projects"][project]
    source_map = {entry["key"]: entry for entry in receipt["sources"]}
    selected_keys = {spec["source_key"], spec.get("reference_key", spec["source_key"])}
    sources = {}
    for key in selected_keys:
        entry = source_map[key]
        path = reference_root / entry["file"]
        if sha256(path) != entry["sha256"]:
            raise ValueError(f"Pinned source changed: {key}")
        sources[key] = ad.read_h5ad(path)
        if list(sources[key].shape) != entry["expected_shape"]:
            raise ValueError(f"Pinned source shape mismatch: {key}")
    autonomy = snapshot.config.get("autonomous_analysis", {})
    source = sources[spec["source_key"]]
    raw_candidates = [meta for meta in snapshot.artifacts.values() if meta.created_by_task == "task_000_ingest"]
    if len(raw_candidates) != 1:
        raise ValueError("One native ingested raw artifact is required")
    raw = registry.load_scdata(raw_candidates[0].uri)
    _, _, clustered = _one_scdata(snapshot, registry, "clustering")
    original_axes = (list(raw.obs.index.astype(str)) == list(source.obs_names.astype(str))
                     and list(raw.var.index.astype(str)) == list(source.var_names.astype(str)))
    forbidden = set(protocol["common"]["forbidden_answer_parameters"])
    supplied = snapshot.config.get("capability_parameters", {})
    found_forbidden = sorted({key for parameters in supplied.values() if isinstance(parameters, dict) for key in parameters if key in forbidden})
    unwanted_fields = [f"obs.{key}" for key in raw.obs if any(word in str(key).casefold() for word in ("cluster", "leiden", "louvain", "pseudotime", "velocity"))]
    if project != "kang2018":
        unwanted_fields += [f"obs.{key}" for key in raw.obs if key in {"cell_type", "celltype", "S_score", "G2M_score"}]
    unwanted_fields += [f"var.{key}" for key in raw.var if str(key).casefold().startswith(("highly_variable", "velocity_", "fit_"))]
    unwanted_fields += [f"obsm.{key}" for key in raw.obsm] + [f"obsp.{key}" for key in raw.obsp]
    unwanted_fields += [f"uns.{key}" for key in raw.uns if key in {"neighbors", "pca", "iroot", "paga"} or str(key).endswith("_colors")]
    common_checks = [check("completed_native_run", summary.get("status") == "success" and summary.get("analysis_request_status") == "complete"),
                     check("deterministic_rule_driven_autonomy", autonomy.get("autonomy_kind") == protocol["common"]["required_autonomy_kind"]),
                     check("fixed_public_input_hash", str(autonomy.get("source_sha256", "")).removeprefix("sha256:") == source_map[spec["source_key"]]["sha256"]),
                     check("frozen_scientific_question_and_context", autonomy.get("question") == spec["native_question"]
                           and snapshot.manifest.biological_design.species == spec["species"]
                           and snapshot.manifest.biological_design.tissue == spec["tissue"]),
                     check("no_answer_configuration", not found_forbidden, forbidden_fields=found_forbidden),
                     check("prior_analysis_fields_quarantined", not unwanted_fields, remaining_fields=unwanted_fields),
                     check("exact_original_raw_axes_and_counts", original_axes and _source_subset(source, raw)),
                     check("retained_cluster_counts_preserved", _source_subset(source, clustered)),
                     check("native_audit_admission", bool(snapshot.audit_contexts) and all(
                         context.receipt_present and context.status == "passed" and context.overall_passed
                         and not context.stop_rule_triggered for context in snapshot.audit_contexts.values())),
                     check("no_reference_simulation", not any(node.is_simulated for node in snapshot.evidence_graph.evidence_nodes.values()))]
    common_checks.append(check("question_only_execution_boundary", True,
                               public_flags=execution_boundary["public_flags"], reference_mount_absent=True,
                               expected_image=execution_boundary["expected_image"], receipt_image=execution_boundary["receipt_image"],
                               attestation_scope=execution_boundary["attestation_scope"]))
    if project == "kang2018":
        source_obs = source.obs
        observations = autonomy.get("observations", {})
        columns = observations.get("obs_column_map", {})
        exact_design = all(columns.get(role) in source_obs and canonical in raw.obs
                           and raw.obs[canonical].astype(str).tolist() == source_obs[columns[role]].astype(str).tolist()
                           for role, canonical in (("donor", "donor"), ("condition", "condition"), ("cell_type", "cell_type")))
        common_checks.append(check("permitted_observed_design_metadata_preserved", exact_design))
    marker_artifacts = _task_artifacts(snapshot, registry, "cluster_markers", kind="Table")
    common_checks.append(check("native_descriptive_marker_table", len(marker_artifacts) == 1))
    if project == "kang2018":
        science = evaluate_kang(snapshot, registry, source, clustered, spec)
    else:
        reference = sources[spec["reference_key"]]
        truth, predicted = label_series(reference, spec["label_key"]), label_series(clustered, "cluster")
        science = clustering_comparison(truth, predicted, spec["thresholds"])
        science["checks"].append(check("all_expected_reference_classes", set(truth.unique()) == set(spec["expected_reference_labels"])))
        if project == "pbmc3k":
            if marker_artifacts:
                markers = marker_recovery(truth, predicted, marker_artifacts[0][2], spec["expected_markers"], spec["thresholds"])
                science["marker_recovery"] = markers
                science["checks"] += markers["checks"]
        else:
            task, _, velocity = _one_scdata(snapshot, registry, "rna_velocity")
            annotated = velocity.uns.get("rna_velocity", {})
            kinetic_counts = all(name in velocity.layers and counts_equal(velocity.layers[name],
                source.layers[original][source.obs_names.get_indexer(velocity.obs.index)][:, source.var_names.get_indexer(velocity.var.index)])
                for name, original in (("spliced_counts", "spliced"), ("unspliced_counts", "unspliced")))
            finite_scores = all(key in velocity.obs and np.isfinite(velocity.obs[key].to_numpy(dtype=float)).all()
                                for key in ("velocity_pseudotime", "root_cells", "end_points", "velocity_confidence"))
            confidence = float(np.median(velocity.obs["velocity_confidence"])) if finite_scores else None
            n_velocity_genes = int(velocity.var["velocity_genes"].sum()) if "velocity_genes" in velocity.var else 0
            science["checks"] += [check("true_stochastic_kinetics", task.method_used == "scvelo_stochastic_v1" and annotated.get("mode") == "stochastic"
                                         and annotated.get("direction_source") == "stochastic_kinetic_model" and annotated.get("caller_root_used") is False),
                                   check("kinetic_raw_counts_preserved", kinetic_counts and _source_subset(source, velocity)),
                                   check("velocity_finite_model_scores", finite_scores),
                                   check("enough_usable_velocity_genes", n_velocity_genes >= spec["thresholds"]["velocity_genes_min"]),
                                   check("local_velocity_confidence", confidence is not None and confidence >= spec["thresholds"]["median_velocity_confidence_min"])]
            if "velocity_graph" not in velocity.uns:
                raise ValueError("No kinetic velocity graph was persisted")
            aligned = truth.reindex(velocity.obs.index.astype(str))
            if aligned.isna().any():
                raise ValueError("Velocity cells do not align to the fixed source")
            direction = velocity_direction(velocity.uns["velocity_graph"], aligned, spec["published_direction_edges"], spec["thresholds"])
            science["direction"] = direction
            science["checks"] += direction["checks"]
            science.update(n_velocity_genes=n_velocity_genes, median_velocity_confidence=confidence,
                           time_interpretation="Unitless model inference; not measured developmental time or experimentally observed fate",
                           state_score_summaries={label: {key: float(velocity.obs.loc[aligned.to_numpy() == label, key].median())
                               for key in ("root_cells", "end_points", "velocity_pseudotime", "velocity_confidence")}
                               for label in sorted(aligned.unique())})
    native_hashes_after = {str(path.relative_to(run_root)): sha256(path) for path in protected if path.is_file()}
    common_checks.append(check("native_run_metadata_unchanged_by_evaluation", native_hashes_before == native_hashes_after))
    common_checks.append(check("execution_receipt_unchanged_after_evaluation",
                               sha256(execution_boundary["path"]) == execution_boundary["sha256"]))
    common_checks.append(check("pinned_references_unchanged_after_evaluation",
        all(sha256(reference_root / source_map[key]["file"]) == source_map[key]["sha256"] for key in selected_keys)
        and sha256(protocol_path) == receipt["blind_protocol_sha256"]
        and sha256(reference_receipt_path) == reference_receipt_hash))
    all_checks = common_checks + science["checks"]
    passed = all(row["passed"] for row in all_checks)
    engineering_passed = all(row["passed"] for row in common_checks)
    return {"schema_version": 1, "project": project, "run_root": str(run_root),
            "evaluated_at_utc": datetime.now(timezone.utc).isoformat(), "protocol_id": protocol["protocol_id"],
            "protocol_sha256": sha256(protocol_path), "source_receipt_sha256": reference_receipt_hash,
            "snapshot_sha256": sha256(run_root / "snapshot.json"),
            "evaluator_sha256": sha256(Path(__file__)), "execution_receipt_sha256": execution_boundary["sha256"],
            "expected_image": execution_boundary["expected_image"], "receipt_image": execution_boundary["receipt_image"],
            "execution_attestation_scope": execution_boundary["attestation_scope"],
            "status": "passed" if passed else "failed", "engineering_status": "passed" if engineering_passed else "failed",
            "scientific_benchmark_status": ("passed" if all(row["passed"] for row in science["checks"]) else "failed") if engineering_passed else "unverified",
            "scope": spec["scientific_scope"], "autonomy_limits": protocol["autonomy_limits"],
            "reference_used_only_after_execution": True,
            "native_config_modified_by_evaluator": native_hashes_before != native_hashes_after,
            "protected_native_file_hashes_before": native_hashes_before, "protected_native_file_hashes_after": native_hashes_after,
            "common_checks": common_checks, "scientific_evaluation": science,
            "failed_checks": [row["name"] for row in all_checks if not row["passed"]]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", choices=("pbmc3k", "kang2018", "pancreas_velocity"), required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--expected-image", required=True, help="Independently inspected frozen image sha256 digest; an engineering provenance check")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.output.exists():
        parser.error("output already exists; choose a new external evaluation receipt")
    if os.name == "nt" and args.output.resolve().drive.upper() == "C:":
        parser.error("public-data evaluation outputs must not be stored on C:")
    try:
        result = evaluate(args.project, args.run_root, args.reference_root, expected_image=args.expected_image)
        exit_code = 0 if result["status"] == "passed" else 1
    except Exception as exc:
        result = {"schema_version": 1, "project": args.project, "status": "invalid",
                  "error": f"{type(exc).__name__}: {exc}", "scientific_benchmark_status": "unverified",
                  "native_config_modified_by_evaluator": False}
        exit_code = 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"project": args.project, "status": result["status"], "output": str(args.output)}, ensure_ascii=False))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
