"""Regression checks for provenance shared by built-in and plugin evidence."""

from types import SimpleNamespace
from unittest.mock import patch

import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport, ValidationSeverity
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.evidence.provenance import normalize_evidence_candidates, resolve_artifact_provenance
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _candidate(contract, result, report, registry):
    return [EvidenceNode(
        evidence_id=f"E_{contract.task_id}", type=EvidenceType.QC_METRICS,
        summary="Observed fixture quality", source_task_id=contract.task_id,
        source_artifact_uris=result.output_artifacts, audit_passed=True,
        # Deliberately incorrect: extractors do not own lineage decisions.
        data_origin_uris=["json://s/unrelated/v1"],
    )]


def _audit(contract, result, registry):
    report = ValidationReport(auditor_name="fixture", target_task_id=contract.task_id)
    report.add_check("fixture_positive", registry.get(result.output_artifacts[0])[1]["value"] > 0,
                     ValidationSeverity.ERROR, "Read persisted output")
    return report


class _Capability(BaseCapability):
    def __init__(self, extractor=_candidate):
        super().__init__("qc", "fixture_v1", descriptor=CapabilityDescriptor(
            "qc", "fixture_v1", input_types=[ArtifactType.JSON], output_types=[ArtifactType.JSON],
            validator=_audit, evidence_extractor=extractor, required_audit_ids=["fixture_positive"],
        ))

    def execute(self, contract, registry):
        registry.register(contract.expected_outputs[0], {"value": 1}, ArtifactType.JSON,
                          "s", contract.task_id, "fixture", parent_uris=contract.input_artifacts)
        return TaskResult(task_id=contract.task_id, status=TaskStatus.SUCCESS,
                          capability="qc", method_used="fixture_v1",
                          input_artifacts=contract.input_artifacts,
                          output_artifacts=contract.expected_outputs)


def _fixture(tmp_path, *, simulated=True, creator="t"):
    registry = ArtifactRegistry(str(tmp_path))
    roots = ["json://s/raw_a/v1", "json://s/raw_b/v1"]
    for uri in roots:
        registry.register(uri, {"value": 1}, ArtifactType.JSON, "s", "ingest", "ingest",
                          summary_metrics={"is_simulated": simulated and uri == roots[0]})
    output = "json://s/output/v1"
    registry.register(output, {"value": 1}, ArtifactType.JSON, "s", creator, "fixture",
                      parent_uris=roots)
    contract = TaskContract(task_id="t", capability="qc", method="fixture_v1",
                            input_artifacts=roots, expected_outputs=[output])
    result = TaskResult(task_id="t", capability="qc", method_used="fixture_v1",
                        status=TaskStatus.SUCCESS, input_artifacts=roots, output_artifacts=[output])
    report = _audit(contract, result, registry)
    return registry, roots, contract, result, report


@pytest.mark.parametrize("extractor", [None, _candidate], ids=["built_in", "descriptor"])
def test_all_extractors_receive_artifact_owned_provenance(tmp_path, extractor):
    registry, roots, contract, result, report = _fixture(tmp_path)
    capabilities = CapabilityRegistry()
    capabilities.register(_Capability(extractor))
    orchestrator = ScientificOrchestrator(registry, capabilities)

    nodes = orchestrator.extract_evidence_from_result(contract, result, report)

    assert nodes
    assert all(node.data_origin_uris == roots for node in nodes)
    assert all(node.is_simulated for node in nodes)
    assert all(node.source_task_id == contract.task_id and node.audit_passed for node in nodes)


def test_normalization_does_not_mutate_plugin_candidate(tmp_path):
    registry, roots, contract, result, report = _fixture(tmp_path)
    candidate = _candidate(contract, result, report, registry)[0]
    normalized = normalize_evidence_candidates([candidate], contract, result, report, registry)[0]
    assert normalized.data_origin_uris == roots
    assert normalized.is_simulated
    assert candidate.data_origin_uris == ["json://s/unrelated/v1"]
    assert not candidate.is_simulated


def test_output_list_cannot_assign_another_tasks_artifact(tmp_path):
    registry, _, contract, result, report = _fixture(tmp_path, creator="other")
    with pytest.raises(ValueError, match="created by a different task"):
        normalize_evidence_candidates(_candidate(contract, result, report, registry),
                                      contract, result, report, registry)


def test_source_verification_cannot_be_self_certified(tmp_path):
    registry, _, contract, result, report = _fixture(tmp_path)
    candidate = _candidate(contract, result, report, registry)[0]
    candidate.source_verified = True
    candidate.metrics["doi"] = "10.1234/a-plausible-identifier"
    with pytest.raises(ValueError, match="independent source verification"):
        normalize_evidence_candidates([candidate], contract, result, report, registry)


@pytest.mark.parametrize("marker", ["result", "candidate"])
def test_additional_simulation_markers_are_preserved(tmp_path, marker):
    registry, _, contract, result, report = _fixture(tmp_path, simulated=False)
    candidate = _candidate(contract, result, report, registry)[0]
    if marker == "result":
        result.metrics["is_simulated"] = True
    else:
        candidate.is_simulated = True
    assert normalize_evidence_candidates([candidate], contract, result, report, registry)[0].is_simulated


def test_failed_check_cannot_be_overridden_by_report_summary(tmp_path):
    registry, _, contract, result, report = _fixture(tmp_path)
    report.add_check("reject", False, ValidationSeverity.ERROR, "Rejected")
    report.overall_passed = True
    assert normalize_evidence_candidates(_candidate(contract, result, report, registry),
                                         contract, result, report, registry) == []


def test_shared_lineage_is_folded_across_plugin_tasks(tmp_path):
    registry, _, contract, result, report = _fixture(tmp_path, simulated=False)
    first = normalize_evidence_candidates(_candidate(contract, result, report, registry),
                                          contract, result, report, registry)[0]
    first.score = 0.2
    contract2 = contract.model_copy(update={"task_id": "t2", "expected_outputs": ["json://s/output2/v1"]})
    result2 = _Capability().execute(contract2, registry)
    second = normalize_evidence_candidates(_candidate(contract2, result2, _audit(contract2, result2, registry), registry),
                                           contract2, result2, _audit(contract2, result2, registry), registry)[0]
    second.score = 0.9
    assert ConfidenceCalculator.calculate([first], []) == ConfidenceCalculator.calculate([first, second], [])


def test_lineage_resolver_rejects_cycles_and_missing_ancestors():
    a, b = "json://s/a/v1", "json://s/b/v1"
    metadata = {
        a: SimpleNamespace(uri=a, parent_uris=[b], summary_metrics={}),
        b: SimpleNamespace(uri=b, parent_uris=[a], summary_metrics={}),
    }
    with pytest.raises(ValueError, match="cycle"):
        resolve_artifact_provenance([a], metadata.__getitem__)
    del metadata[b]
    with pytest.raises(KeyError):
        resolve_artifact_provenance([a], metadata.__getitem__)


def test_lineage_resolver_collects_shared_roots_and_intermediate_markers():
    a, b, c, d = [f"json://s/{name}/v1" for name in "abcd"]
    metadata = {
        a: SimpleNamespace(uri=a, parent_uris=[], summary_metrics={}),
        b: SimpleNamespace(uri=b, parent_uris=[a], summary_metrics={"is_simulated": True}),
        c: SimpleNamespace(uri=c, parent_uris=[a], summary_metrics={}),
        d: SimpleNamespace(uri=d, parent_uris=[b, c], summary_metrics={}),
    }
    provenance = resolve_artifact_provenance([d, b], metadata.__getitem__)
    assert provenance.root_uris == (a,)
    assert provenance.is_simulated


def test_plugin_run_propagates_simulation_into_claim_and_summary(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    root = "json://s/raw/v1"
    registry.register(root, {"value": 1}, ArtifactType.JSON, "s", "ingest", "ingest",
                      summary_metrics={"is_simulated": True})
    capabilities = CapabilityRegistry()
    capabilities.register(_Capability())
    orchestrator = ScientificOrchestrator(registry, capabilities)
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="test"))
    contract = TaskContract(task_id="t", capability="qc", method="fixture_v1",
                            input_artifacts=[root], expected_outputs=["json://s/output/v1"])
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[contract]):
        summary = orchestrator.run_study(manifest, {"method_overrides": {"qc": "fixture_v1"}})
    assert summary["status"] == "success", summary
    assert summary["is_simulated"]
    node = next(iter(orchestrator.evidence_graph.evidence_nodes.values()))
    claim = next(iter(orchestrator.evidence_graph.claim_nodes.values()))
    assert node.data_origin_uris == [root]
    assert claim.is_simulated
    assert claim.statement.startswith("[SIMULATED DATA]")
