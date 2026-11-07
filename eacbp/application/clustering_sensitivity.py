"""Prespecified, reference-free clustering perturbations and paired descriptions.

This policy never chooses a winning resolution, replaces the primary result, or
equates reproducibility on shared cells with biological validation.
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
from scipy import sparse
from scipy.optimize import linear_sum_assignment


CLUSTERING_SENSITIVITY_POLICY_VERSION = "bounded_reference_free_clustering_v1"
_RECEIPTS = ("clustering_sensitivity_plan", "clustering_sensitivity_branch",
             "qc_sensitivity_plan", "qc_sensitivity_branch")
_SOURCE_FACTS = ("source_sha256", "observations", "requirement_contract", "input_preparation")
_CLUSTERING_DEFAULTS = {"resolution": 1.0, "random_seed": 42, "n_neighbors": 15}
_INTERPRETATION = ("Prespecified dependent sensitivity descriptions on shared cells; "
                   "not independent biological validation, cell-type truth, viability, "
                   "or a rule for selecting a winning branch. The primary is unchanged.")


def _fingerprint(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(",", ":"), allow_nan=False).encode("utf-8")).hexdigest()


def _digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def _config_scope(config: Mapping[str, Any]) -> dict[str, Any]:
    clean = deepcopy(dict(config))
    settings = clean.get("autonomous_analysis", {})
    for name in _RECEIPTS:
        settings.pop(name, None)
    return clean


def clustering_config_sha256(config: Mapping[str, Any]) -> str:
    """Hash the actual configuration excluding only four self-referential receipts."""
    return _fingerprint(_config_scope(config))


def _scientific_config_sha256(config: Mapping[str, Any]) -> str:
    clean = _config_scope(config)
    settings = clean.get("autonomous_analysis", {})
    for name in _SOURCE_FACTS:
        settings.pop(name, None)
    parameters = clean.setdefault("capability_parameters", {}).setdefault("clustering", {})
    for key, value in _CLUSTERING_DEFAULTS.items():
        parameters.setdefault(key, value)
    return _fingerprint(clean)


def clustering_scientific_config_sha256(config: Mapping[str, Any]) -> str:
    """Expose the immutable scientific scope used in native execution receipts."""
    return _scientific_config_sha256(config)


def raw_count_fingerprint(matrix) -> str:
    """Bind integer count values, shape and dtype independently of sparse layout."""
    canonical = sparse.csr_matrix(matrix, copy=True)
    canonical.sum_duplicates()
    canonical.sort_indices()
    canonical.eliminate_zeros()
    values = canonical.data
    if (not np.isfinite(values).all() or np.any(values < 0)
            or not np.equal(values, np.rint(values)).all()):
        raise ValueError("Clustering sensitivity requires finite nonnegative integer counts")
    digest = hashlib.sha256()
    digest.update(json.dumps({"shape": list(canonical.shape), "dtype": canonical.dtype.str},
                             sort_keys=True).encode("utf-8"))
    for array in (canonical.indptr.astype("<i8"), canonical.indices.astype("<i8"),
                  canonical.data.astype(canonical.dtype.newbyteorder("<"), copy=False)):
        digest.update(np.ascontiguousarray(array).tobytes())
    return digest.hexdigest()


def feature_axis_fingerprint(var: pd.DataFrame, gene_name_source: str) -> str:
    """Bind the complete ordered feature identities and gene names, not an HVG mask."""
    ids = var.index.astype(str).tolist()
    if len(set(ids)) != len(ids):
        raise ValueError("String-normalized feature IDs must be unique")
    genes = (var.index.astype(str).str.strip().tolist() if gene_name_source == "var.index"
             else var[gene_name_source[4:]].astype(str).str.strip().tolist())
    return _fingerprint({"feature_ids": ids, "genes": genes})


def _raw_source(payload, observations):
    source = observations["counts_source"]
    if source == "layers.counts":
        return payload.layers["counts"], payload.var
    if source == "raw.X":
        return payload.raw["X"], payload.raw["var"]
    return payload.X, payload.var


def _retained_ids(payload, observations, qc_parameters) -> list[str]:
    matrix, var = _raw_source(payload, observations)
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
    keep = detected >= qc_parameters["min_genes"]
    genes = (var.index.astype(str).str.strip().to_numpy() if observations["gene_name_source"] == "var.index"
             else var[observations["gene_name_source"][4:]].astype(str).str.strip().to_numpy())
    mt = (var["mt"].to_numpy(dtype=bool) if "mt" in var else
          np.char.startswith(np.char.lower(genes.astype(str)), "mt-"))
    if mt.any():
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        mito = np.divide(np.asarray(matrix[:, mt].sum(axis=1)).ravel() * 100., totals,
                         out=np.zeros(len(totals)), where=totals > 0)
        keep &= mito <= qc_parameters["max_mito_pct"]
    return [cell for cell, kept in zip(payload.obs.index.astype(str), keep) if kept]


def freeze_clustering_sensitivity_plan(*, payload, observations: Mapping[str, Any],
                                      source_sha256: str, baseline_config: Mapping[str, Any]) -> dict[str, Any]:
    """Freeze at most eight perturbations using counts and observed donor metadata."""
    from eacbp.application.autonomous import observe_input, prepare_input
    from eacbp.capabilities.biological_units import complete_identity

    _digest(source_sha256, "source_sha256")
    settings = baseline_config.get("autonomous_analysis", {})
    if settings.get("intent") != "exploration" or baseline_config.get("method_profile") != "standard":
        raise ValueError("Clustering sensitivity v1 supports standard native exploration only")
    if baseline_config.get("advanced_analysis") or any(
            name not in {"cluster_markers", "qc_review"} and value is not False
            for name, value in baseline_config.get("analysis_extensions", {}).items()):
        raise ValueError("Clustering sensitivity v1 does not support advanced or unrelated extensions")
    if any(baseline_config.get("analysis_extensions", {}).get(name) in (None, False)
           for name in ("cluster_markers", "qc_review")):
        raise ValueError("Clustering sensitivity v1 requires native cluster markers and QC review")
    if baseline_config.get("method_overrides", {}).get("clustering", "scanpy_leiden_umap_v1") != "scanpy_leiden_umap_v1":
        raise ValueError("Clustering sensitivity requires native Scanpy Leiden")
    question = settings.get("question")
    if observe_input(payload, question=question) != observations or settings.get("observations") != observations:
        raise ValueError("Raw input observations differ from the frozen baseline")
    if settings.get("source_sha256") != source_sha256:
        raise ValueError("Baseline and raw source identities differ")
    parameters = baseline_config.get("capability_parameters", {})
    configured = {key: parameters.get("clustering", {}).get(key, value)
                  for key, value in _CLUSTERING_DEFAULTS.items()}
    if parameters.get("clustering", {}).get("marker_reference"):
        raise ValueError("Reference-free clustering sensitivity cannot use a marker-reference panel")
    if configured != _CLUSTERING_DEFAULTS or parameters.get("integration", {}).get("random_seed", 42) != 42:
        raise ValueError("Primary clustering must use resolution 1, 15 neighbors and seed 42; PCA/integration seed remains 42")
    qc = {key: parameters.get("qc", {}).get(key) for key in ("min_genes", "max_mito_pct")}
    if qc != {key: observations["qc"][key] for key in qc}:
        raise ValueError("Primary QC must retain the original complete-input 3 scaled MAD thresholds")
    ids = payload.obs.index.astype(str).tolist()
    if len(set(ids)) != len(ids):
        raise ValueError("String-normalized cell IDs must be unique")
    matrix, var = _raw_source(payload, observations)
    axis_sha256 = feature_axis_fingerprint(var, observations["gene_name_source"])
    prepared, _ = prepare_input(payload, settings)
    canonical_observations = observe_input(prepared, question=question)

    def spec(branch_id, positions, clustering, kind, reason, omitted_donor=None):
        selected = [ids[int(index)] for index in positions]
        full = len(positions) == len(ids)
        subset = payload if full else prepared.subset_obs(np.isin(np.arange(len(ids)), positions))
        branch_observations = (dict(observations) if full else
                               observe_input(subset, question=question) if len(selected) >= 3 else None)
        count_observations = branch_observations or canonical_observations
        branch_matrix, _ = _raw_source(subset, count_observations)
        retained = _retained_ids(prepared if full else subset,
                                 canonical_observations if full else count_observations, qc)
        config = _config_scope(baseline_config)
        config.setdefault("capability_parameters", {}).setdefault("clustering", {}).update(clustering)
        return {"branch_id": branch_id, "kind": kind, "source_kind": "original" if full else "derived_subset",
                "selected_positions": [int(index) for index in positions],
                "selected_ids": selected, "selected_ids_sha256": _fingerprint(selected),
                "n_source_cells": len(selected), "n_features": int(matrix.shape[1]),
                "raw_counts_sha256": raw_count_fingerprint(branch_matrix),
                "feature_axis_sha256": axis_sha256, "expected_observations": deepcopy(branch_observations),
                "clustering_parameters": dict(clustering), "qc_parameters": dict(qc),
                "scientific_config_sha256": _scientific_config_sha256(config),
                "predicted_retained_ids": retained, "predicted_retained_ids_sha256": _fingerprint(retained),
                "planned_status": "planned" if len(retained) >= 3 else "not_estimable",
                "reason": reason if len(retained) >= 3 else "Fewer than three cells survive fixed primary QC; no rescue is permitted",
                "omitted_donor": omitted_donor}

    all_positions = np.arange(len(ids))
    primary = spec("primary", all_positions, _CLUSTERING_DEFAULTS, "primary", "Unmodified frozen primary analysis")
    branches = [
        spec("resolution_0_5", all_positions, {**_CLUSTERING_DEFAULTS, "resolution": .5}, "resolution", "Only Leiden resolution changes"),
        spec("resolution_2_0", all_positions, {**_CLUSTERING_DEFAULTS, "resolution": 2.}, "resolution", "Only Leiden resolution changes"),
        spec("clustering_seed_17", all_positions, {**_CLUSTERING_DEFAULTS, "random_seed": 17}, "clustering_stage_seed", "Neighbor graph, Leiden and UMAP seed changes; PCA/integration remains seed 42"),
        spec("uniform_sample_80pct", np.sort(np.random.default_rng(20261003).choice(
            len(ids), size=int(math.floor(.8 * len(ids))), replace=False)), _CLUSTERING_DEFAULTS,
             "source_subsample", "Uniform raw-source 80% cell subsample with seed 20261003, before QC; original ordering and all features preserved"),
    ]
    donor_column = observations.get("obs_column_map", {}).get("donor")
    donor_record = {"source_column": donor_column, "donors": [], "status": "not_eligible",
                    "reason": "Observed biological donor identity is unavailable; sample or library labels are not promoted",
                    "donor_independence_established": False}
    if donor_column:
        values = payload.obs[donor_column]
        if not complete_identity(values):
            donor_record["reason"] = "Observed donor metadata are incomplete"
        else:
            donors = sorted(values.astype(str).str.strip().unique().tolist())
            donor_record["donors"] = donors
            if not 2 <= len(donors) <= 4:
                donor_record["reason"] = "The bounded complete donor leave-one-out policy requires two to four observed donors"
            else:
                donor_record.update(status="planned", reason="Every observed donor is omitted once; labels do not establish sampling independence")
                labels = values.astype(str).str.strip().to_numpy()
                slugs = {donor: donor if re.fullmatch(r"[A-Za-z][A-Za-z0-9_-]{0,63}", donor)
                         else "donor_" + hashlib.sha256(donor.encode("utf-8")).hexdigest()[:16]
                         for donor in donors}
                folded = [slug.casefold() for slug in slugs.values()]
                for donor, slug in list(slugs.items()):
                    if folded.count(slug.casefold()) > 1:
                        slugs[donor] = "donor_" + hashlib.sha256(donor.encode("utf-8")).hexdigest()[:16]
                if len({slug.casefold() for slug in slugs.values()}) != len(donors):
                    raise ValueError("Observed donor identities collide in their bounded filesystem namespace")
                for donor in donors:
                    branches.append(spec(f"leave_out_{slugs[donor]}", np.flatnonzero(labels != donor),
                                         _CLUSTERING_DEFAULTS, "donor_leave_one_out",
                                         "One complete observed donor is removed from raw source before fixed QC", donor))
    plan = {"schema_version": 1, "policy_version": CLUSTERING_SENSITIVITY_POLICY_VERSION,
            "source_sha256": source_sha256, "baseline_config_sha256": clustering_config_sha256(baseline_config),
            "counts_source": observations["counts_source"], "cell_ids": ids,
            "n_features": int(matrix.shape[1]), "feature_axis_sha256": axis_sha256,
            "primary": primary, "branches": branches, "donor_perturbation": donor_record,
            "subsample_seed": 20261003, "small_cluster_minimum": 20,
            "winner_selection_permitted": False, "biological_stability_claim_supported": False,
            "interpretation": _INTERPRETATION}
    plan["plan_sha256"] = _fingerprint(plan)
    plan["plan_id"] = "clustering_sensitivity_" + plan["plan_sha256"][:16]
    return plan


def validate_clustering_sensitivity_plan(plan: Mapping[str, Any]) -> None:
    clean = {key: value for key, value in plan.items() if key not in {"plan_id", "plan_sha256"}}
    digest = _fingerprint(clean)
    if (plan.get("policy_version") != CLUSTERING_SENSITIVITY_POLICY_VERSION
            or plan.get("plan_sha256") != digest
            or plan.get("plan_id") != "clustering_sensitivity_" + digest[:16]):
        raise ValueError("Clustering sensitivity plan identity or policy changed")


def _branch_spec(plan, branch_id):
    validate_clustering_sensitivity_plan(plan)
    values = [plan["primary"], *plan["branches"]]
    spec = next((item for item in values if item["branch_id"] == branch_id), None)
    if spec is None:
        raise ValueError("Unknown frozen clustering sensitivity branch")
    return spec


def clustering_branch_config(plan: Mapping[str, Any], branch_id: str,
                             baseline_config: Mapping[str, Any], *,
                             source_sha256: str | None = None) -> dict[str, Any]:
    """Bind an actual derived-file digest after preparation, before native execution."""
    from eacbp.application.analysis_requirements import freeze_question_requirements

    spec = _branch_spec(plan, branch_id)
    if clustering_config_sha256(baseline_config) != plan["baseline_config_sha256"]:
        raise ValueError("Baseline configuration changed after the clustering plan was frozen")
    if spec["planned_status"] != "planned":
        raise ValueError("A not-estimable frozen branch cannot be executed")
    if spec["source_kind"] == "original":
        source_sha256 = plan["source_sha256"] if source_sha256 is None else source_sha256
        if source_sha256 != plan["source_sha256"]:
            raise ValueError("A full-source branch must use the original source")
    else:
        _digest(source_sha256, "derived source_sha256")
        if source_sha256 == plan["source_sha256"]:
            raise ValueError("A subset cannot claim the original complete source digest")
    result = _config_scope(baseline_config)
    result.setdefault("capability_parameters", {}).setdefault("clustering", {}).update(spec["clustering_parameters"])
    settings = result["autonomous_analysis"]
    observation = deepcopy(spec["expected_observations"])
    settings.update(source_sha256=source_sha256, observations=observation,
                    requirement_contract=freeze_question_requirements(
                        question=settings["question"], observations=observation, source_sha256=source_sha256),
                    input_preparation={key: observation[key] for key in ("counts_source", "gene_name_source", "obs_column_map")})
    settings["clustering_sensitivity_branch"] = {
        "plan_id": plan["plan_id"], "plan_sha256": plan["plan_sha256"], "branch_id": branch_id,
        "parent_source_sha256": plan["source_sha256"], "source_sha256": source_sha256,
        "selected_ids_sha256": spec["selected_ids_sha256"], "raw_counts_sha256": spec["raw_counts_sha256"],
        "feature_axis_sha256": spec["feature_axis_sha256"], "qc_parameters": dict(spec["qc_parameters"]),
        "clustering_parameters": dict(spec["clustering_parameters"]), "interpretation": _INTERPRETATION}
    validate_clustering_branch_config(plan, branch_id, result)
    return result


def validate_clustering_branch_config(plan, branch_id, config) -> None:
    """Reject changed methods, QC, normalisation, integration or source facts."""
    from eacbp.application.analysis_requirements import freeze_question_requirements

    spec = _branch_spec(plan, branch_id)
    if _scientific_config_sha256(config) != spec["scientific_config_sha256"]:
        raise ValueError("Branch scientific configuration differs from the frozen allowed delta")
    settings = config.get("autonomous_analysis", {})
    source_sha256 = _digest(settings.get("source_sha256"), "branch source_sha256")
    if spec["source_kind"] == "original" and source_sha256 != plan["source_sha256"]:
        raise ValueError("Full-source branch source identity changed")
    observation = spec["expected_observations"]
    preparation = {key: observation[key] for key in ("counts_source", "gene_name_source", "obs_column_map")}
    if settings.get("observations") != observation or settings.get("input_preparation") != preparation:
        raise ValueError("Derived source observations/preparation differ from the frozen subset")
    expected = freeze_question_requirements(question=settings.get("question"), observations=observation,
                                            source_sha256=source_sha256)
    if settings.get("requirement_contract") != expected:
        raise ValueError("Branch question requirement contract is not bound to its actual source")


def _verified_labels(plan, spec, receipt):
    if receipt.get("admitted") is not True:
        return None
    if spec["planned_status"] != "planned":
        raise ValueError("A not-estimable frozen run cannot claim admitted clustering")
    if receipt.get("run_status") != "success":
        raise ValueError("An admitted clustering run must have successful native status")
    for key in ("selected_ids_sha256", "raw_counts_sha256", "feature_axis_sha256"):
        if receipt.get(key) != spec[key]:
            raise ValueError(f"Admitted branch {key} differs from the frozen source subset")
    if receipt.get("raw_selected_ids") != spec["selected_ids"]:
        raise ValueError("Admitted raw source membership/order differs from the frozen subset")
    if receipt.get("qc_parameters") != spec["qc_parameters"]:
        raise ValueError("Clustering perturbation changed the frozen primary QC thresholds")
    if receipt.get("clustering_parameters") != spec["clustering_parameters"]:
        raise ValueError("Executed clustering parameters differ from the frozen branch")
    if receipt.get("scientific_config_sha256") != spec["scientific_config_sha256"]:
        raise ValueError("Admitted scientific configuration changed beyond the frozen delta")
    actual_config = _digest(receipt.get("config_sha256"), "actual branch config_sha256")
    if actual_config != receipt.get("expected_config_sha256"):
        raise ValueError("Executed native configuration differs from its pre-execution binding")
    source = _digest(receipt.get("source_sha256"), "actual branch source_sha256")
    if source != receipt.get("expected_source_sha256") or (
            spec["source_kind"] == "original" and source != plan["source_sha256"]):
        raise ValueError("Executed raw source differs from its pre-execution binding")
    if receipt.get("retained_ids") != spec["predicted_retained_ids"]:
        raise ValueError("Admitted QC retained membership/order differs from the fixed raw-count rule")
    labels = receipt.get("clusters")
    if (not isinstance(labels, pd.Series) or not labels.index.is_unique
            or not labels.index.astype(str).is_unique or labels.isna().any()
            or set(labels.index.astype(str)) != set(spec["predicted_retained_ids"])):
        raise ValueError("Clusters must uniquely cover every admitted retained cell")
    normalized = labels.copy()
    normalized.index = normalized.index.astype(str)
    return normalized.astype(str).reindex(spec["predicted_retained_ids"])


def _paired_descriptions(primary, branch, minimum):
    shared = primary.index[primary.index.isin(branch.index)]
    main_clusters = sorted(primary.unique())
    branch_clusters = sorted(branch.unique())
    counts = pd.crosstab(primary.loc[shared], branch.loc[shared]).reindex(
        index=main_clusters, columns=branch_clusters, fill_value=0).to_numpy(dtype=np.int64)
    shared_main, shared_branch = counts.sum(axis=1), counts.sum(axis=0)
    unions = shared_main[:, None] + shared_branch[None, :] - counts
    jaccard = np.divide(counts, unions, out=np.zeros(counts.shape, dtype=float), where=unions > 0)
    assignment = {}
    if len(shared):
        active_rows, active_cols = np.flatnonzero(shared_main), np.flatnonzero(shared_branch)
        rows, cols = linear_sum_assignment(-jaccard[np.ix_(active_rows, active_cols)])
        assignment = {int(active_rows[row]): int(active_cols[col]) for row, col in zip(rows, cols)}
    descriptions = []
    for row, cluster in enumerate(main_clusters):
        total = int(primary.eq(cluster).sum())
        common = int(shared_main[row])
        col = assignment.get(row)
        intersection = int(counts[row, col]) if col is not None else 0
        no_overlap = common == 0
        descriptions.append({"primary_cluster": cluster, "primary_total_cells": total,
                             "shared_cells": common, "shared_coverage": common / total,
                             "small_support": common < minimum,
                             "small_primary_support": total < minimum,
                             "support_status": "not_evaluable" if no_overlap else "small_support" if common < minimum else "described",
                             "matched_branch_cluster": branch_clusters[col] if col is not None else None,
                             "matched_branch_shared_cells": (None if no_overlap else int(shared_branch[col]) if col is not None else 0),
                             "matching_status": "no_shared_cells" if no_overlap else "matched" if col is not None else "unmatched",
                             "intersection_cells": intersection,
                             "union_cells": (None if no_overlap else int(unions[row, col]) if col is not None else common),
                             "shared_jaccard": None if no_overlap else float(jaccard[row, col]) if col is not None else 0.,
                             "shared_recall": None if no_overlap else intersection / common,
                             "shared_precision": (None if no_overlap else intersection / int(shared_branch[col]) if col is not None else 0.)})
    branch_descriptions = [{"branch_cluster": cluster, "shared_cells": int(shared_branch[column]),
                            "matched_primary_cluster": main_clusters[next(row for row, col in assignment.items() if col == column)]
                            if column in assignment.values() else None,
                            "matching_status": "matched" if column in assignment.values() else "unmatched"}
                           for column, cluster in enumerate(branch_clusters)]
    result = {"n_shared_cells": len(shared), "shared_cell_ids_sha256": _fingerprint(shared.tolist()),
              "primary_total_cells": len(primary), "branch_total_cells": len(branch),
              "primary_shared_fraction": len(shared) / len(primary) if len(primary) else None,
              "branch_shared_fraction": len(shared) / len(branch) if len(branch) else None,
              "has_small_primary_clusters": any(row["small_primary_support"] for row in descriptions),
              "has_small_shared_clusters": any(row["small_support"] for row in descriptions),
              "primary_total_clusters": len(main_clusters),
              "branch_total_clusters": len(branch_clusters), "primary_shared_clusters": int(np.count_nonzero(shared_main)),
              "branch_shared_clusters": int(np.count_nonzero(shared_branch)), "cluster_matches": descriptions,
              "branch_cluster_matches": branch_descriptions,
              "matching_rule": "Sort observed common labels; one-to-one Hungarian assignment maximizes total common-cell Jaccard; assigned zero-overlap pairs retain their actual denominators",
              "interpretation": _INTERPRETATION}
    if len(shared) < 2:
        return {**result, "status": "not_estimable", "shared_cell_ari": None, "shared_cell_nmi": None,
                "degenerate_partition": True, "reason": "Fewer than two shared retained cells"}
    nonempty = counts[np.ix_(shared_main > 0, shared_branch > 0)]
    n = len(shared)
    pairs = lambda array: float(np.sum(array.astype(float) * (array - 1) / 2))
    expected = pairs(shared_main) * pairs(shared_branch) / (n * (n - 1) / 2)
    maximum = (pairs(shared_main) + pairs(shared_branch)) / 2
    ari = 1. if maximum == expected else (pairs(counts) - expected) / (maximum - expected)
    probability = nonempty / n
    row_p, col_p = probability.sum(axis=1), probability.sum(axis=0)
    positive = probability > 0
    mutual = float(np.sum(probability[positive] * np.log((probability / (row_p[:, None] * col_p[None, :]))[positive])))
    entropy = lambda values: -float(np.sum(values * np.log(values)))
    denominator = (entropy(row_p) + entropy(col_p)) / 2
    nmi = 1. if denominator == 0 else mutual / denominator
    return {**result, "status": "computed", "shared_cell_ari": float(ari), "shared_cell_nmi": float(nmi),
            "degenerate_partition": len(row_p) < 2 or len(col_p) < 2,
            "reason": "A single-cluster partition does not provide population-resolution evidence" if len(row_p) < 2 or len(col_p) < 2 else None}


def compare_clustering_sensitivity_runs(plan: Mapping[str, Any], *, baseline, branches) -> dict[str, Any]:
    """Describe only admitted frozen runs; keep missing/failed branches explicit."""
    validate_clustering_sensitivity_plan(plan)
    if set(branches) - {item["branch_id"] for item in plan["branches"]}:
        raise ValueError("Unexpected unfrozen clustering branch")
    primary = _verified_labels(plan, plan["primary"], baseline)
    results = []
    for spec in plan["branches"]:
        receipt = branches.get(spec["branch_id"])
        item = {"branch_id": spec["branch_id"], "kind": spec["kind"], "planned_status": spec["planned_status"],
                "clustering_parameters": spec["clustering_parameters"], "omitted_donor": spec["omitted_donor"]}
        if spec["planned_status"] != "planned":
            item.update(status="not_estimable", reason=spec["reason"])
        elif receipt is None:
            item.update(status="not_available", reason="The frozen native branch was not supplied")
        else:
            labels = _verified_labels(plan, spec, receipt)
            if labels is None:
                item.update(status="not_available", reason="The native branch was not independently admitted")
            elif primary is None:
                item.update(status="not_available", reason="The primary native result was not independently admitted")
            else:
                item.update(_paired_descriptions(primary, labels, plan["small_cluster_minimum"]))
        if receipt is not None:
            item["run_status"] = receipt.get("run_status")
            item["provenance"] = deepcopy(receipt.get("provenance", {}))
        results.append(item)
    complete = primary is not None and all(item["status"] == "computed" for item in results if item["planned_status"] == "planned")
    return {"schema_version": 1, "policy_version": plan["policy_version"], "plan_id": plan["plan_id"],
            "plan_sha256": plan["plan_sha256"], "status": "completed" if complete else "partial",
            "completion_scope": "Execution and paired descriptions of all estimable frozen branches only",
            "branches": results, "donor_perturbation": deepcopy(plan["donor_perturbation"]),
            "winner_selected": False, "primary_replaced": False,
            "biological_stability_claim_supported": False, "interpretation": _INTERPRETATION}
