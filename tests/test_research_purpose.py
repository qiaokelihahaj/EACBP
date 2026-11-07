"""Research purpose constrains the plan and preserves unanswered requirements."""
import pytest

from eacbp.application.analysis_requirements import freeze_plan_requirements, evaluate_requirement_coverage
from eacbp.orchestrator.planning import preview_study_plan
from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec
from eacbp.schemas.task import TaskContract


def study():
    return StudyManifest(study_id="purpose", biological_design=BiologicalDesign(species="human", tissue="synthetic"),
                         data=DataSpec(raw_artifact_uri="adata://purpose/raw/v1"))


@pytest.mark.parametrize("purpose,expected", [("overview", set()), ("deg", {"deg"}),
    ("abundance", {"differential_abundance"}), ("trajectory", {"trajectory_inference"}), ("batch", set())])
def test_purpose_removes_unrequested_inference(purpose, expected):
    preview = preview_study_plan(study(), {"research_purpose": purpose, "method_profile": "baseline"})
    assert preview["valid"], preview
    caps = {t["capability"] for t in preview["tasks"]}
    assert caps == {"dataset_audit", "qc", "normalization", "integration", "clustering"} | expected
    ids = {t["task_id"] for t in preview["tasks"]}
    assert all(set(t["depends_on"]) <= ids for t in preview["tasks"])


def test_legacy_unspecified_scope_keeps_existing_plan():
    preview = preview_study_plan(study(), {"method_profile": "baseline"})
    assert {"deg", "differential_abundance", "trajectory_inference"} <= {t["capability"] for t in preview["tasks"]}


def test_explicit_consumer_keeps_its_required_producers():
    preview = preview_study_plan(study(), {"research_purpose": "overview", "method_profile": "baseline", "include_knowledge": True})
    assert preview["valid"], preview
    caps = {t["capability"] for t in preview["tasks"]}
    assert {"knowledge_retrieval", "deg"} <= caps
    assert "differential_abundance" not in caps


def test_missing_primary_capability_remains_an_unanswered_requirement():
    contract = freeze_plan_requirements(purpose="annotation", tasks=[], source_sha256="a" * 64)
    result = evaluate_requirement_coverage(contract, [], pipeline_status="success")
    assert result["analysis_request_status"] == "partial"
    assert result["requirement_completion_checks"][0]["missing_capabilities"] == ["cell_annotation"]


def test_admitted_different_task_cannot_satisfy_frozen_request():
    task = TaskContract(task_id="requested", capability="clustering", expected_outputs=["adata://purpose/cluster/v1"])
    contract = freeze_plan_requirements(purpose="overview", tasks=[task], source_sha256="a" * 64)
    result = {"task_id": "unrelated", "capability": "clustering", "status": "success", "method_used": "kmeans_marker_embedding_v1", "output_artifacts": task.expected_outputs}
    coverage = evaluate_requirement_coverage(contract, [(result, {"task_id": "unrelated"})])
    assert coverage["analysis_request_status"] == "partial"


def test_same_task_id_cannot_change_frozen_population_binding():
    task = TaskContract(task_id="subset", capability="subset_cells", parameters={"cell_type": "A"},
                        expected_outputs=["adata://purpose/subset/v1"])
    contract = freeze_plan_requirements(purpose="overview", tasks=[task], source_sha256="a" * 64)
    changed = task.model_copy(deep=True)
    changed.parameters["cell_type"] = "B"
    result = {"task_id": "subset", "capability": "subset_cells", "status": "success", "method_used": "subset_cells_v1",
              "output_artifacts": task.expected_outputs, "metrics": {"subset_cell_type": "B"}}
    coverage = evaluate_requirement_coverage(contract, [(result, changed.model_dump())])
    assert coverage["requirement_completion_checks"][0]["status"] != "passed"


def test_competing_question_scopes_are_rejected():
    preview = preview_study_plan(study(), {"research_purpose": "overview", "autonomous_analysis": {"intent": "exploration"}})
    assert not preview["valid"]


def test_annotation_requires_a_configured_annotation_task():
    preview = preview_study_plan(study(), {"research_purpose": "annotation"})
    assert not preview["valid"]
    assert "cell_annotation" in str(preview["errors"])


def test_explicit_target_fate_keeps_trajectory_producers():
    manifest = study()
    manifest.biological_design.target_cell_types = ["A", "B"]
    preview = preview_study_plan(manifest, {"research_purpose": "overview", "method_profile": "standard",
        "target_parameters": {target: {"trajectory_inference": {"root_cell_id": target + "1"},
                                       "fate_mapping": {"terminal_states": {"end": [target + "2"]}}} for target in ["A", "B"]}})
    assert preview["valid"], preview
    assert sum(task["capability"] == "fate_mapping" for task in preview["tasks"]) == 2
    assert sum(task["capability"] == "trajectory_inference" for task in preview["tasks"]) == 2


def test_frozen_task_completion_does_not_invent_a_question_dimension_requirement():
    task = TaskContract(task_id="fate", capability="fate_mapping", expected_outputs=["table://purpose/fate/v1"])
    clustering = TaskContract(task_id="cluster", capability="clustering", expected_outputs=["adata://purpose/clustering/v1"])
    contract = freeze_plan_requirements(purpose="overview", tasks=[clustering, task], source_sha256="a" * 64)
    admitted = [({"task_id": t.task_id, "capability": t.capability, "status": "success", "method_used": "tested",
                 "output_artifacts": t.expected_outputs, "metrics": {}}, t.model_dump()) for t in [clustering, task]]
    assert evaluate_requirement_coverage(contract, admitted)["analysis_request_status"] == "complete"


@pytest.mark.parametrize("options,capability", [({"run_spacell_adapter": True}, "spacell_microenvironment_analysis"),
    ({"adapters": ["spacell"]}, "spacell_microenvironment_analysis"),
    ({"run_compound_perturbation": True}, "compound_perturbation_simulation")])
def test_explicit_unavailable_branch_is_not_silently_dropped(options, capability):
    preview = preview_study_plan(study(), {"research_purpose": "overview", **options})
    assert not preview["valid"]
    assert capability in str(preview["errors"])


def test_explicit_manifest_hypothesis_keeps_knowledge_and_its_dependencies():
    manifest = study()
    manifest.hypotheses.user_provided = ["Inspect the supplied hypothesis"]
    preview = preview_study_plan(manifest, {"research_purpose": "overview"})
    assert preview["valid"], preview
    assert {"knowledge_retrieval", "deg"} <= {t["capability"] for t in preview["tasks"]}


def test_advanced_deg_respects_explicit_sensitivity_opt_out():
    preview = preview_study_plan(study(), {"research_purpose": "deg", "advanced_analysis": True,
                                          "analysis_extensions": {"donor_sensitivity": False}})
    assert preview["valid"], preview
    assert "donor_sensitivity" not in {t["capability"] for t in preview["tasks"]}
