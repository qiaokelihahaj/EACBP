"""QC perturbations remain prespecified, paired and independently auditable."""

from copy import deepcopy
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.application.autonomous import observe_input, prepare_input
from eacbp.application.qc_sensitivity import (
    compare_deg_tables, compare_qc_sensitivity_runs, freeze_qc_sensitivity_plan,
    qc_branch_config, qc_config_sha256,
)
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.computational import ComputationalValidator
from eacbp.capabilities.qc import QCCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


SOURCE = "a" * 64


def _input(*, groups=False, mitochondrial=True):
    detected = [0, 10, 40, 80, 100, 130, 150, 170, 190, 210, 230, 250, 270, 290, 310]
    counts = np.zeros((len(detected), 400), dtype=np.int32)
    for i, n in enumerate(detected):
        counts[i, :n] = 1
        if n and i >= 4:
            counts[i, 0] += i * 2
    counts[1, 1] += 1
    obs = pd.DataFrame(index=[f"c{i}" for i in range(len(detected))])
    if groups:
        obs["celltype"] = ["low_RNA"] * 5 + ["other"] * 10
    genes = ["MT-A" if mitochondrial else "G0", *[f"G{i}" for i in range(1, 400)]]
    return SCData(sparse.csr_matrix(counts), obs, pd.DataFrame(index=genes))


def _freeze(payload, *, groups=False):
    question = "Explore existing observed cell populations" if groups else "Discover cell populations"
    observation = observe_input(payload, question=question)
    config = {"mode": "real", "method_profile": "standard", "threads": 1,
              "autonomous_analysis": {"question": question, "source_sha256": SOURCE,
                                      "observations": observation},
              "analysis_extensions": {"cluster_markers": True},
              "capability_parameters": {"qc": {name: observation["qc"][name]
                                                 for name in ("min_genes", "max_mito_pct")}}}
    plan = freeze_qc_sensitivity_plan(payload=payload, observations=observation,
                                      source_sha256=SOURCE, baseline_config=config)
    return plan, config


def _native_receipt(tmp_path, payload, config):
    prepared, _ = prepare_input(payload, config["autonomous_analysis"])
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://qc_sensitivity/input/v1"
    registry.register(uri, prepared, ArtifactType.ANNDATA, "qc_sensitivity", "input", "test")
    parameters = config["capability_parameters"]["qc"]
    contract = TaskContract(task_id="qc", capability="qc", input_artifacts=[uri], parameters=parameters)
    result = QCCapability().execute(contract, registry)
    audit = ComputationalValidator().audit(contract, result, registry)
    assert audit.overall_passed
    output = registry.load_scdata(result.output_artifacts[0])
    return {"source_sha256": SOURCE, "config_sha256": qc_config_sha256(config),
            "parameters": deepcopy(parameters), "qc_admitted": True,
            "cell_disposition": output.uns["qc"]["cell_disposition"],
            "run_status": "success", "summaries": {"qc": result.metrics}}


def _pair(tmp_path, *, groups=False):
    payload = _input(groups=groups)
    plan, config = _freeze(payload, groups=groups)
    baseline = _native_receipt(tmp_path / "primary", payload, config)
    branches = {spec["branch_id"]: _native_receipt(
        tmp_path / spec["branch_id"], payload, qc_branch_config(plan, spec["branch_id"], config))
        for spec in plan["branches"] if spec["planned_status"] == "planned"}
    return payload, plan, config, baseline, branches


def test_plan_preserves_primary_and_freezes_two_raw_distribution_perturbations():
    payload = _input()
    original = payload.to_dict()
    plan, config = _freeze(payload)
    unchanged = deepcopy(config)
    assert [branch["scaled_mad_multiplier"] for branch in plan["branches"]] == [2, 4]
    assert plan["primary"]["parameters"] == config["capability_parameters"]["qc"]
    strict, wide = plan["branches"]
    assert strict["parameters"]["min_genes"] >= plan["primary"]["parameters"]["min_genes"] >= wide["parameters"]["min_genes"] >= 1
    assert strict["parameters"]["max_mito_pct"] <= plan["primary"]["parameters"]["max_mito_pct"] <= wide["parameters"]["max_mito_pct"]
    branch_config = qc_branch_config(plan, "wider_4mad", config)
    assert config == unchanged
    assert branch_config["capability_parameters"]["qc"] == wide["parameters"]
    assert branch_config["analysis_extensions"] == config["analysis_extensions"]
    np.testing.assert_array_equal(payload.X.toarray(), original["X"].toarray())
    assert json.loads(json.dumps(plan, allow_nan=False)) == plan


def test_annotation_reference_and_desired_labels_cannot_retune_thresholds():
    payload = _input(groups=True)
    main, _ = _freeze(payload)
    poisoned = payload.copy()
    poisoned.obs["celltype"] = "desired_reference_answer"
    poisoned.obs["percent_mito"] = 100.
    poisoned.uns["reference"] = {"desired_retention": 1.0, "winning_qc": 1}
    poisoned.obsm["X_umap"] = np.ones((payload.n_obs, 2))
    plan, _ = _freeze(poisoned)
    assert plan == main
    assert not plan["group_retention"]["authorized"]


@pytest.mark.parametrize("mutation", ["source", "primary_threshold", "scope", "observations"])
def test_freezing_rejects_changed_source_primary_or_observation(mutation):
    payload = _input()
    _, config = _freeze(payload)
    observation = deepcopy(config["autonomous_analysis"]["observations"])
    source = SOURCE
    if mutation == "source":
        source = "b" * 64
    elif mutation == "primary_threshold":
        config["capability_parameters"]["qc"]["min_genes"] += 1
    elif mutation == "scope":
        config["threads"] = float("nan")
    else:
        observation["n_cells"] -= 1
    with pytest.raises(ValueError):
        freeze_qc_sensitivity_plan(payload=payload, observations=observation,
                                    source_sha256=source, baseline_config=config)


def test_frozen_plan_and_scope_cannot_be_edited_or_replaced():
    plan, config = _freeze(_input())
    changed = deepcopy(config)
    changed["analysis_extensions"]["cluster_markers"] = False
    with pytest.raises(ValueError, match="scope/config changed"):
        qc_branch_config(plan, "wider_4mad", changed)
    changed_plan = deepcopy(plan)
    changed_plan["branches"][1]["parameters"]["min_genes"] = 1
    with pytest.raises(ValueError, match="identity"):
        qc_branch_config(changed_plan, "wider_4mad", config)
    with pytest.raises(ValueError, match="Unknown"):
        qc_branch_config(plan, "manual_reference_best", config)
    config["autonomous_analysis"]["qc_sensitivity_plan"] = plan
    branch = qc_branch_config(plan, "wider_4mad", config)
    assert qc_config_sha256(branch) == plan["branches"][1]["config_sha256"]


def test_identical_distributions_and_unestimable_branches_are_explicit():
    uniform = SCData(np.ones((6, 5)), pd.DataFrame(index=list("abcdef")), pd.DataFrame(index=list("ABCDE")))
    plan, config = _freeze(uniform)
    assert all(spec["planned_status"] == "duplicate_primary" for spec in plan["branches"])
    with pytest.raises(ValueError, match="duplicate"):
        qc_branch_config(plan, "stricter_2mad", config)
    empty = SCData(np.zeros((6, 5)), uniform.obs, uniform.var)
    plan, _ = _freeze(empty)
    assert plan["primary"]["planned_status"] == "not_estimable"
    assert plan["primary"]["predicted_retained_cells"] == 0
    assert all(spec["planned_status"] == "not_estimable" for spec in plan["branches"])
    assert all(spec["parameters"]["min_genes"] == 1 for spec in plan["branches"])


def test_native_admitted_runs_compare_full_dispositions_and_authorized_groups(tmp_path):
    payload, plan, _, baseline, branches = _pair(tmp_path, groups=True)
    before = deepcopy(baseline)
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches,
                                         observed_groups=payload.obs["celltype"])
    assert report["status"] == "completed"
    assert len(report["branches"]) == 2
    assert all(len(row["cell_dispositions"]) == payload.n_obs for row in report["branches"])
    wide = report["branches"][1]
    assert wide["n_added"] > 0 and wide["n_removed"] == 0
    assert sum(row["n_cells"] for row in wide["disposition_transitions"]) == payload.n_obs
    assert wide["cell_dispositions"][0]["primary_retained"] is False
    assert wide["cell_dispositions"][0]["branch_retained"] is False
    assert wide["cell_dispositions"][0]["mito_pct"] is None
    assert wide["observed_group_retention"][0]["observed_group"] == "low_RNA"
    assert wide["group_retention_status"] == "computed"
    assert "not independent" in report["interpretation"]
    assert baseline["parameters"] == before["parameters"]
    np.testing.assert_array_equal(baseline["cell_disposition"]["retained"], before["cell_disposition"]["retained"])
    json.dumps(report, allow_nan=False)


@pytest.mark.parametrize("mutation", ["source", "config", "parameters", "cell_ids", "raw_metrics", "availability", "retained", "reason", "flags"])
def test_comparison_rejects_mismatched_or_forged_branch_evidence(tmp_path, mutation):
    _, plan, _, baseline, branches = _pair(tmp_path)
    receipt = branches["wider_4mad"]
    table = receipt["cell_disposition"]
    if mutation == "source":
        receipt["source_sha256"] = "b" * 64
    elif mutation == "config":
        receipt["config_sha256"] = "b" * 64
    elif mutation == "parameters":
        receipt["parameters"]["min_genes"] = 0
    elif mutation == "cell_ids":
        table["ids"][0] = "unexpected_cell"
    elif mutation == "raw_metrics":
        table["total_counts"][1] += 1
    elif mutation == "availability":
        table["mito_pct_available"][0] = True
    elif mutation == "retained":
        table["retained"][0] = True
    elif mutation == "reason":
        table["reason"][0] = "retained"
    else:
        table["retained"] = table["retained"].astype(str)
    with pytest.raises(ValueError):
        compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)


def test_group_labels_are_optional_and_cannot_supply_posthoc_reference_answers(tmp_path):
    payload, plan, _, baseline, branches = _pair(tmp_path)
    with pytest.raises(ValueError, match="not authorized"):
        compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches,
                                     observed_groups=pd.Series("reference", index=payload.obs.index))
    payload, plan, _, baseline, branches = _pair(tmp_path / "authorized", groups=True)
    altered = payload.obs["celltype"].copy()
    altered.iloc[0] = "answer"
    with pytest.raises(ValueError, match="labels changed"):
        compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches, observed_groups=altered)
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    assert report["branches"][0]["group_retention_status"] == "not_supplied"


def test_missing_and_unadmitted_runs_remain_partial(tmp_path):
    _, plan, _, baseline, branches = _pair(tmp_path)
    branches.pop("stricter_2mad")
    branches["wider_4mad"]["qc_admitted"] = False
    branches["wider_4mad"]["run_status"] = "failed"
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    assert report["status"] == "partial"
    assert [row["status"] for row in report["branches"]] == ["missing", "not_estimable"]
    assert report["branches"][1]["run_status"] == "failed"


def test_failed_downstream_run_cannot_make_the_sensitivity_round_complete(tmp_path):
    _, plan, _, baseline, branches = _pair(tmp_path)
    branches["wider_4mad"]["run_status"] = "failed"
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    assert report["status"] == "partial"
    assert report["branches"][1]["status"] == "compared"
    assert report["branches"][1]["run_status"] == "failed"


def test_shared_cell_ari_handles_relabeling_and_excluded_cells(tmp_path):
    _, plan, _, baseline, branches = _pair(tmp_path)
    main = pd.DataFrame(baseline["cell_disposition"])
    primary_ids = main.loc[main.retained, "ids"]
    baseline["clusters"] = pd.Series([str(i % 2) for i in range(len(primary_ids))], index=primary_ids)
    for receipt in branches.values():
        frame = pd.DataFrame(receipt["cell_disposition"])
        ids = frame.loc[frame.retained, "ids"]
        shared = baseline["clusters"].reindex(ids)
        receipt["clusters"] = shared.map({"0": "renamed_A", "1": "renamed_B"}).fillna("extra")
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    for row in report["branches"]:
        assert row["clustering"]["shared_cell_ari"] == 1.
        assert row["clustering"]["n_shared_cells"] <= len(primary_ids)
    branches["wider_4mad"]["clusters"].loc["unknown_cell"] = "fake"
    with pytest.raises(ValueError, match="exactly"):
        compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)


def test_unknown_mito_and_authoritative_count_layer_are_preserved(tmp_path):
    payload = _input(mitochondrial=False)
    payload.layers["counts"] = payload.X.copy()
    payload.X = payload.X.astype(float) / 3
    plan, config = _freeze(payload)
    assert plan["counts_source"] == "layers.counts"
    assert plan["mitochondrial_metric_source"] == "unavailable"
    assert plan["primary"]["parameters"]["max_mito_pct"] == 100.
    baseline = _native_receipt(tmp_path / "main", payload, config)
    branches = {spec["branch_id"]: _native_receipt(tmp_path / spec["branch_id"], payload,
                                                    qc_branch_config(plan, spec["branch_id"], config))
                for spec in plan["branches"] if spec["planned_status"] == "planned"}
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    assert report["status"] == "completed"
    assert all(row["mito_pct"] is None for row in report["branches"][0]["cell_dispositions"])


def test_raw_feature_space_and_explicit_mito_annotation_use_the_same_rule(tmp_path):
    source = _input()
    raw_var = source.var.copy()
    raw_var["mt"] = [True] + [False] * (len(raw_var) - 1)
    raw_var.index = [f"ENSG{i}" for i in range(len(raw_var))]
    payload = SCData(np.full((source.n_obs, 3), .3), source.obs,
                     pd.DataFrame(index=["normalized_A", "normalized_B", "normalized_C"]),
                     raw={"X": source.X.copy(), "var": raw_var})
    plan, config = _freeze(payload)
    assert plan["counts_source"] == "raw.X"
    assert plan["n_features"] == 400
    assert plan["mitochondrial_metric_source"] == "var.mt"
    baseline = _native_receipt(tmp_path / "main", payload, config)
    branches = {spec["branch_id"]: _native_receipt(tmp_path / spec["branch_id"], payload,
                                                    qc_branch_config(plan, spec["branch_id"], config))
                for spec in plan["branches"] if spec["planned_status"] == "planned"}
    report = compare_qc_sensitivity_runs(plan, baseline=baseline, branches=branches)
    assert report["status"] == "completed"
    assert report["branches"][1]["n_added"] > 0


def _deg(genes, effects, *, fdr=None, target="T", tested="stim", reference="ctrl", statuses=None, alpha=.05):
    frame = pd.DataFrame({"gene": genes, "contrast_tested_level": tested,
                          "contrast_reference_level": reference, "log2_fold_change": effects,
                          "fdr_q_value": fdr if fdr is not None else [.01] * len(genes),
                          "status": statuses if statuses is not None else ["estimated"] * len(genes), "alpha": alpha})
    return {"target_cell_type": target, "table": frame, "admitted": True}


def test_deg_full_union_exposes_sign_reversals_missing_and_unestimated_genes():
    first = _deg(["same", "flip", "zero", "unestimated", "missing"], [1., -1., 0., np.nan, 2.],
                 fdr=[.01, .01, .01, np.nan, .01], statuses=["estimated"] * 3 + ["not_estimated", "estimated"])
    second = _deg(["same", "flip", "zero", "unestimated", "new"], [2., 2., 2., 1., 2.],
                  fdr=[.01, .4, .01, np.nan, .01])
    result = compare_deg_tables([first], [second])
    assert result["direction_counts"] == {"same_nonzero_direction": 1, "opposite_direction": 1,
                                            "zero_in_either": 1, "not_estimable": 1,
                                            "missing_in_primary": 1, "missing_in_branch": 1}
    assert result["n_fdr_status_changes"] == 1
    assert result["n_paired_estimated"] == 3
    missing = next(row for row in result["rows"] if row["gene"] == "missing")
    assert missing["branch_log2_fold_change"] is None and missing["branch_row_status"] == "missing"
    unestimated = next(row for row in result["rows"] if row["gene"] == "unestimated")
    assert unestimated["primary_log2_fold_change"] is None and unestimated["primary_significant"] is None
    json.dumps(result, allow_nan=False)


def test_deg_oriented_contrast_and_population_cannot_silently_pair():
    first = _deg(["G"], [1.])
    second = _deg(["G"], [-1.], tested="ctrl", reference="stim")
    result = compare_deg_tables([first], [second])
    assert result["n_paired_estimated"] == 0
    assert result["direction_counts"]["missing_in_primary"] == 1
    assert result["direction_counts"]["missing_in_branch"] == 1
    with pytest.raises(ValueError, match="Duplicate"):
        compare_deg_tables([first, first], [second])
    second = _deg(["G"], [1.], alpha=.1)
    with pytest.raises(ValueError, match="alpha changed"):
        compare_deg_tables([first], [second])
    first["admitted"] = False
    assert compare_deg_tables([first], [])["status"] == "not_available"


@pytest.mark.parametrize("field,value", [("log2_fold_change", np.inf), ("fdr_q_value", -.1), ("fdr_q_value", 1.1), ("status", "made_up")])
def test_deg_invalid_numbers_or_status_are_rejected(field, value):
    invalid = _deg(["G"], [1.])
    invalid["table"].loc[0, field] = value
    with pytest.raises(ValueError):
        compare_deg_tables([invalid], [])
