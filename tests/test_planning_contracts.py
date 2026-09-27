"""Extension scope and static method contracts agree before computation."""

import pytest
from pydantic import ValidationError
from unittest.mock import patch

from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.planning import build_study_tasks, preview_study_plan
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract


def study(targets=("Microglia", "Neurons")):
    return StudyManifest(study_id="contracts", biological_design=BiologicalDesign(
        species="human", tissue="brain", target_cell_types=list(targets)))


class Extension(BaseCapability):
    def __init__(self, scope, plan_factory):
        super().__init__("extension", "extension_v1", descriptor=CapabilityDescriptor(
            "extension", "extension_v1", scope=scope, plan_factory=plan_factory))

    def execute(self, *args):
        raise AssertionError("planning must not execute capabilities")


def target_plan(*, manifest, parameters, tasks):
    subset = next(task for task in tasks if task.capability == "subset_cells")
    output = f"json://{manifest.study_id}/extension/v1"
    return [
        TaskContract(task_id="extension_first", capability="extension",
                     input_artifacts=subset.expected_outputs,
                     expected_outputs=[output], parameters=parameters),
        TaskContract(task_id="extension_second", capability="extension",
                     input_artifacts=[output], depends_on=["extension_first"],
                     expected_outputs=[f"json://{manifest.study_id}/extension_summary/v1"]),
    ]


def test_registered_extension_uses_same_target_expansion_and_parameter_overrides():
    registry = create_default_capability_registry()
    registry.register(Extension("per_target", target_plan))
    tasks = build_study_tasks(study(), {
        "analysis_extensions": {"extension": {"count": 2}},
        "target_parameters": {"Neurons": {"extension": {"count": 3}}},
    }, registry)
    selected = {task.task_id: task for task in tasks if task.capability == "extension"}
    assert set(selected) == {f"extension_{stage}_{target}"
                            for stage in ("first", "second") for target in ("microglia", "neurons")}
    for target, count in (("microglia", 2), ("neurons", 3)):
        first, second = selected[f"extension_first_{target}"], selected[f"extension_second_{target}"]
        assert first.input_artifacts == [f"adata://contracts/microglia_subset/{target}/v5"]
        assert first.parameters["count"] == count
        assert first.depends_on == [f"task_006_subset_{target}"]
        assert second.depends_on == [first.task_id]
        assert second.input_artifacts == first.expected_outputs
        assert first.parameters["target_branch"] == target
    assert registry.describe("extension").to_dict()["scope"] == "per_target"


def test_shared_scope_overrides_legacy_target_parameter_markers():
    registry = create_default_capability_registry()
    registry.register(Extension("shared", lambda **kwargs: TaskContract(
        task_id="shared_extension", capability="extension",
        parameters={"target_dependent": True, "target_cell_type": "label_only"},
        expected_outputs=["json://contracts/shared/v1"])))
    tasks = build_study_tasks(study(), {"analysis_extensions": {"extension": True}}, registry)
    selected = [task for task in tasks if task.capability == "extension"]
    assert len(selected) == 1
    assert selected[0].task_id == "shared_extension"
    assert selected[0].parameters["target_cell_type"] == "label_only"
    assert "target_branch" not in selected[0].parameters


def test_shared_extension_must_explicitly_choose_inputs_from_target_branches():
    registry = create_default_capability_registry()
    registry.register(Extension("shared", target_plan))
    with pytest.raises(ValueError, match="Shared task.*unexpanded per-target"):
        build_study_tasks(study(), {"analysis_extensions": {"extension": True}}, registry)


def test_per_target_extension_keeps_legacy_single_target_uri():
    registry = create_default_capability_registry()
    registry.register(Extension("per_target", target_plan))
    tasks = build_study_tasks(study(("Microglia",)), {"analysis_extensions": {"extension": True}}, registry)
    task = next(task for task in tasks if task.task_id == "extension_first")
    assert task.expected_outputs == ["json://contracts/extension/v1"]
    assert task.input_artifacts == ["adata://contracts/microglia_subset/v5"]


@pytest.mark.parametrize("targets", [(), ("Microglia",), ("Microglia", "Neurons")])
def test_unified_assembly_preserves_existing_builtin_contracts(targets):
    manifest = study(targets)
    state = {"method_profile": "baseline", "analysis_extensions": {"donor_sensitivity": True}}
    old = ComputationalDAGPlanner.order_tasks(ComputationalDAGPlanner.build_study_plan(manifest, state))
    new = build_study_tasks(manifest, state, create_default_capability_registry())
    assert [task.model_dump() for task in new] == [task.model_dump() for task in old]


@pytest.mark.parametrize("capability,parameters,message", [
    ("qc", {"min_gene": 10}, "min_gene"),
    ("qc", {"min_cells": 1}, "min_cells"),
    ("qc", {"min_genes": "10"}, "min_genes"),
    ("qc", {"max_mito_pct": 101}, "max_mito_pct"),
    ("normalization", {"target_sum": 0}, "target_sum"),
    ("normalization", {"n_top_genes": True}, "n_top_genes"),
    ("deg", {"allow_x_as_counts": "false"}, "allow_x_as_counts"),
    ("deg", {"condition_a": "A"}, "supplied together"),
    ("deg", {"condition_a": "A", "condition_b": "A"}, "must be different"),
    ("deg", {"pseudocount": float("nan")}, "pseudocount"),
])
def test_preview_rejects_invalid_method_parameters(capability, parameters, message):
    preview = preview_study_plan(study(("Microglia",)), {
        "method_profile": "baseline", "capability_parameters": {capability: parameters}})
    assert not preview["valid"]
    assert any(message in error["message"] for error in preview["errors"])


def test_parameter_normalization_preserves_omissions_context_and_legacy_aliases():
    registry = create_default_capability_registry()
    descriptor = registry.describe("normalization", "sc_normalize_log1p_v1")
    assert descriptor.validate_parameters({}) == {}
    context = {"study_id": "s", "target_branch": "neurons", "fdr_family": "target:neurons",
               "target_provenance": {"shared_preprocessing": ["qc"]}, "target_sum": 1000.0}
    assert descriptor.validate_parameters(context) == context
    assert descriptor.validate_types
    assert registry.describe("deg", "deg_pseudobulk_v1").scope == "per_target"
    assert registry.describe("deg").validate_parameters({"batch_col": "batch"}) == {"batch_col": "batch"}
    with pytest.raises(ValidationError):
        registry.describe("deg", "pydeseq2_pseudobulk_v1").validate_parameters({"paired": "false"})


def test_cellbender_static_options_rejected_before_resource_hashing(tmp_path, monkeypatch):
    import eacbp.orchestrator.advanced_plan as advanced_plan
    source, executable = tmp_path / "raw.h5", tmp_path / "cellbender.exe"
    source.write_bytes(b"raw")
    executable.write_bytes(b"executable")
    monkeypatch.setattr(advanced_plan, "pin_resource_files", lambda _: pytest.fail("must reject before hashing"))
    for extra in ({"timeout_sec": 0}, {"allow_overwrite": "false"}, {"timeot_sec": 2},
                  {"extra_args": ["--output=elsewhere.h5"]}):
        preview = preview_study_plan(study(()), {"analysis_extensions": {"background_removal": {
            "unfiltered_input_path": str(source), "executable": str(executable),
            "output_path": str(tmp_path / "output.h5"), **extra}}})
        assert not preview["valid"]
        assert preview["errors"]


def test_scope_declaration_rejects_typographical_errors():
    with pytest.raises(ValueError, match="scope"):
        CapabilityDescriptor("extension", "extension_v1", scope="per-target")


def test_registered_target_extension_executes_and_audits_each_branch(tmp_path):
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.auditor.base import ValidationReport, ValidationSeverity
    from eacbp.capabilities.registry import CapabilityRegistry
    from eacbp.orchestrator.loop import ScientificOrchestrator
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.evidence import EvidenceNode, EvidenceType
    from eacbp.schemas.task import TaskResult, TaskStatus

    def plan(**kwargs):
        return TaskContract(task_id="extension", capability="extension",
                            expected_outputs=["json://contracts/extension/v1"])

    def audit(contract, result, registry):
        payload = registry.load_payload(result.output_artifacts[0])
        report = ValidationReport(auditor_name="target_fixture", target_task_id=contract.task_id)
        report.add_check("target_matches", payload["target"] == contract.parameters["target_cell_type"],
                         ValidationSeverity.ERROR, "Persisted target matches its contract")
        return report

    def evidence(contract, result, report, registry):
        return [EvidenceNode(evidence_id=f"E_{contract.task_id}", type=EvidenceType.QC_METRICS,
                             summary="Independent target fixture", source_task_id=contract.task_id,
                             source_artifact_uris=result.output_artifacts, audit_passed=True)]

    class ExecutableExtension(Extension):
        def execute(self, contract, registry):
            output = contract.expected_outputs[0]
            registry.register(output, {"target": contract.parameters["target_cell_type"]},
                              ArtifactType.JSON, "contracts", contract.task_id, "target_fixture")
            return TaskResult(task_id=contract.task_id, capability="extension", method_used="extension_v1",
                              status=TaskStatus.SUCCESS, output_artifacts=[output])

    capability = ExecutableExtension("per_target", plan)
    capability.descriptor.validator = audit
    capability.descriptor.evidence_extractor = evidence
    capability.descriptor.required_audit_ids = ("target_matches",)
    capabilities = CapabilityRegistry()
    capabilities.register(capability)
    orchestrator = ScientificOrchestrator(ArtifactRegistry(str(tmp_path)), capabilities)
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[]):
        result = orchestrator.run_study(study(), {"analysis_extensions": {"extension": True}})
    assert result["status"] == "success", result
    assert len(orchestrator.task_history) == 2
    nodes = list(orchestrator.evidence_graph.evidence_nodes.values())
    assert len(nodes) == 2
    assert {node.biological_context["target_branch"] for node in nodes} == {"microglia", "neurons"}
    assert {node.source_artifact_uris[0] for node in nodes} == {
        "json://contracts/extension/microglia/v1", "json://contracts/extension/neurons/v1"}


@pytest.mark.parametrize("capability,method", [
    ("qc", "sc_qc_v1"), ("normalization", "library_size_log1p_v1"),
    ("deg", "donor_pseudobulk_welch_v1"), ("deg", "pydeseq2_pseudobulk_v1"),
])
def test_matrix_contracts_accept_spatial_payloads_but_reject_json(tmp_path, capability, method):
    import numpy as np
    import pandas as pd
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.capabilities.sc_data import SCData
    from eacbp.schemas.artifact import ArtifactType

    artifacts = ArtifactRegistry(str(tmp_path))
    spatial = "adata://contracts/spatial/v1"
    wrong_type = "json://contracts/non_matrix/v1"
    data = SCData(np.ones((3, 2)), pd.DataFrame(index=["a", "b", "c"]),
                  pd.DataFrame(index=["g1", "g2"]), obsm={"spatial": np.ones((3, 2))})
    artifacts.register(spatial, data, ArtifactType.SPATIAL_DATA, "contracts", "input", "fixture")
    artifacts.register(wrong_type, {}, ArtifactType.JSON, "contracts", "input", "fixture")
    capabilities = create_default_capability_registry()
    contract = TaskContract(task_id="matrix", capability=capability, method=method, input_artifacts=[spatial])
    capabilities.prepare_contract(contract, artifacts)
    contract.input_artifacts = [wrong_type]
    with pytest.raises(ValueError, match="incompatible inputs"):
        capabilities.prepare_contract(contract, artifacts)


@pytest.mark.parametrize("capability,parameters,output", [
    ("qc", {"min_genes": 1}, "adata://contracts/qc/v1"),
    ("normalization", {}, "adata://contracts/normalized/v2"),
])
def test_spatial_preprocessing_preserves_coordinates(tmp_path, capability, parameters, output):
    import numpy as np
    import pandas as pd
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.capabilities.sc_data import SCData
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskStatus

    coordinates = np.array([[0., 0.], [1., 0.], [0., 1.]])
    data = SCData(np.array([[1., 2.], [2., 3.], [3., 1.]]), pd.DataFrame(index=["a", "b", "c"]),
                  pd.DataFrame(index=["g1", "g2"]), obsm={"spatial": coordinates})
    artifacts = ArtifactRegistry(str(tmp_path))
    source = "adata://contracts/raw/v1"
    artifacts.register(source, data, ArtifactType.SPATIAL_DATA, "contracts", "input", "fixture")
    contract = TaskContract(task_id="preprocess", capability=capability, input_artifacts=[source],
                            expected_outputs=[output], parameters=parameters)
    result = create_default_capability_registry().execute_contract(contract, artifacts)
    assert result.status == TaskStatus.SUCCESS, result.error_message
    assert artifacts.get_metadata(source).type == ArtifactType.SPATIAL_DATA
    assert artifacts.get_metadata(output).type == ArtifactType.ANNDATA
    assert np.array_equal(artifacts.load_payload(output).obsm["spatial"], coordinates)
