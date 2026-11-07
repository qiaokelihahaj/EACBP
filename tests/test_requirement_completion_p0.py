"""Scientific scope must use admitted immutable results in every lifecycle."""
from copy import deepcopy
from types import SimpleNamespace

import pytest

from eacbp.application.analysis_requirements import freeze_question_requirements
from eacbp.application.study_service import _autonomous_request_summary
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport
from eacbp.capabilities import create_default_capability_registry
from eacbp.schemas.audit import stable_fingerprint
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import ScientificResult, ScientificResultStatus, TaskResult, TaskStatus


def coverage_case(question="Compare existing observed cell populations"):
    source = "a" * 64
    observation = {"conditions": ["ctrl", "stim"], "condition_metadata_complete": True,
                   "contrast": {"tested": "stim", "reference": "ctrl"},
                   "use_observed_cell_types": True, "observed_cell_types": ["B", "T"],
                   "obs_column_map": {"donor": "donor"},
                   "targets": [{"cell_type": target, "eligible": True} for target in ("B", "T")]}
    frozen = freeze_question_requirements(question=question, observations=observation, source_sha256=source)
    manifest = StudyManifest(study_id="coverage", biological_design=BiologicalDesign(
        species="human", tissue="blood", target_cell_types=["B", "T"]))
    config = {"autonomous_analysis": {"question": question, "intent": "condition_comparison",
              "source_sha256": source, "scientific_gaps": [], "requirement_contract": frozen},
              "analysis_extensions": {"cluster_markers": True}}
    history, reports, records, commits = [], [], [], {}

    def add(capability, target=None, *, network_kind=None):
        task_id = "_".join(filter(None, (capability, target, network_kind)))
        parameters = {"cell_type": target} if capability == "subset_cells" else {"target_cell_type": target} if target else {}
        if capability in {"deg", "donor_sensitivity", "functional_activity"}:
            parameters.update(condition_a="stim", condition_b="ctrl")
        if network_kind:
            parameters["network_kind"] = network_kind
        inputs = [f"table://coverage/subset_cells_{target}/v1"] if capability in {"deg", "donor_sensitivity", "functional_activity"} else []
        result = TaskResult(task_id=task_id, capability=capability, method_used="unit_test_method",
                            status=TaskStatus.SUCCESS, input_artifacts=inputs,
                            output_artifacts=[f"table://coverage/{task_id}/v1"],
                            metrics={"all_requested_fits_completed": True} if capability == "donor_sensitivity" else {})
        if capability in {"deg", "donor_sensitivity", "functional_activity"}:
            result.scientific_result = ScientificResult(status=ScientificResultStatus.ESTIMATED_INCONCLUSIVE,
                summary="Estimated with no supported features", n_features_tested=3, n_features_supported=0)
        context = {"capability": capability, "method": result.method_used, "parameters": parameters}
        payload = ArtifactRegistry._result_payload(result)
        record = SimpleNamespace(task_id=task_id, signature=task_id, contract=context,
            contract_fingerprint=stable_fingerprint(context), result_fingerprint=ArtifactRegistry._result_identity(payload),
            input_artifacts=list(inputs), output_artifacts=list(result.output_artifacts), artifact_hashes={},
            overall_passed=True, stop_rule_triggered=False)
        history.append(result)
        reports.append(ValidationReport(auditor_name="independent", target_task_id=task_id, overall_passed=True))
        records.append(record)
        commits[task_id] = {"result": deepcopy(payload), "output_hashes": {}}
        return result, record

    for capability in ("dataset_audit", "qc", "normalization", "integration", "clustering", "cluster_markers"):
        add(capability)
    for target in ("B", "T"):
        for capability in ("subset_cells", "deg", "donor_sensitivity"):
            add(capability, target)
    registry = SimpleNamespace(list_audit_records=lambda **_: records, get_task_commit=lambda signature: commits.get(signature))
    orchestrator = SimpleNamespace(artifact_registry=registry, task_history=history, audit_reports=reports,
                                  capability_registry=create_default_capability_registry())
    return SimpleNamespace(config=config, manifest=manifest, orchestrator=orchestrator, history=history,
                           records=records, commits=commits, add=add)


def certify(case):
    return _autonomous_request_summary({"status": "success"}, manifest=case.manifest,
                                      config=case.config, orchestrator=case.orchestrator)


def test_requested_complex_dimensions_cannot_be_completed_by_deg_alone():
    case = coverage_case("Compare pathways, TF regulators and communication in existing observed cell populations")
    summary = certify(case)
    assert summary["status"] == "success" and summary["analysis_request_status"] == "partial"
    unmet = {check["dimension"] for check in summary["requirement_completion_checks"] if check["status"] != "passed"}
    assert unmet == {"pathway", "regulator", "communication"}
    assert summary["analysis_completion_checks"] == summary["requirement_completion_checks"]


@pytest.mark.parametrize("state", ["estimated_inconclusive", "not_estimable", "assumptions_failed"])
def test_admitted_scientific_states_remain_distinct(state):
    case = coverage_case()
    result = next(result for result in case.history if result.task_id == "deg_B")
    result.scientific_result = ScientificResult(status=state, summary="Audited scientific state",
        n_features_tested=3 if state == "estimated_inconclusive" else 0, n_features_supported=0)
    case.commits[result.task_id]["result"] = ArtifactRegistry._result_payload(result)
    summary = certify(case)
    check = next(check for check in summary["requirement_completion_checks"]
                 if check["dimension"] == "population_expression" and check["target_cell_type"] == "B")
    assert check["status"] == ("passed" if state == "estimated_inconclusive" else state)
    assert summary["analysis_request_status"] == ("complete" if state == "estimated_inconclusive" else "partial")


@pytest.mark.parametrize("tamper", ["input", "scientific_state", "metrics", "missing_commit", "output_hash"])
def test_old_admission_cannot_be_borrowed_by_modified_results(tamper):
    case = coverage_case()
    result = next(result for result in case.history if result.task_id == "deg_B")
    if tamper == "input":
        result.input_artifacts = ["table://coverage/subset_cells_T/v1"]
    elif tamper == "scientific_state":
        result.scientific_result.n_features_tested = 100
    elif tamper == "metrics":
        result.metrics["condition_a"] = "unobserved"
    elif tamper == "missing_commit":
        case.commits.pop(result.task_id)
    else:
        case.commits[result.task_id]["output_hashes"] = {result.output_artifacts[0]: "changed"}
    assert certify(case)["analysis_request_status"] == "partial"


def test_pathway_resource_cannot_answer_regulator_requirement():
    case = coverage_case("Compare pathways and TF regulators in existing observed cell populations")
    for target in ("B", "T"):
        case.add("functional_activity", target, network_kind="pathway")
    checks = certify(case)["requirement_completion_checks"]
    assert all(check["status"] == "passed" for check in checks if check["dimension"] == "pathway")
    assert all(check["status"] != "passed" for check in checks if check["dimension"] == "regulator")


@pytest.mark.parametrize("field", ["question", "source_sha256"])
def test_frozen_contract_must_bind_saved_request(field):
    case = coverage_case()
    case.config["autonomous_analysis"][field] = "changed"
    with pytest.raises(ValueError, match="saved question or source"):
        certify(case)


def test_legacy_scope_is_not_reparsed_or_extended():
    case = coverage_case("Compare pathways and TF regulators in existing observed cell populations")
    request = case.config["autonomous_analysis"]
    request.pop("requirement_contract")
    before = deepcopy(case.config)
    summary = certify(case)
    assert summary["analysis_request_status"] == "complete"
    assert "requirement_contract_id" not in summary
    assert case.config == before


def test_post_commit_dataset_planning_does_not_invalidate_scientific_admission():
    case = coverage_case()
    case.history[0].metrics["planning_decisions"] = [{"decision": "valid_post_commit_plan"}]
    assert certify(case)["analysis_request_status"] == "complete"


def test_generic_replicate_field_cannot_establish_biological_donors():
    import pandas as pd
    from eacbp.application.autonomous import observe_input
    from eacbp.capabilities.sc_data import SCData
    import numpy as np

    obs = pd.DataFrame({"replicate": ["r1", "r2", "r3", "r4"],
                        "condition": ["ctrl", "ctrl", "stim", "stim"], "cell_type": ["T"] * 4})
    observation = observe_input(SCData(np.ones((4, 3)), obs, pd.DataFrame(index=["A", "B", "C"])),
                                question="Compare existing observed cell populations")
    assert observation["obs_column_map"]["donor"] is None
    assert not observation["targets"]
