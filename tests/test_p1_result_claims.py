"""Negative scientific results retain provenance without becoming positive claims."""
import pytest

from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.graph import EvidenceGraph
from eacbp.schemas.evidence import (
    EvidenceNode, EvidenceType, EvidencePolarity, LanguageTier, ClaimType,
)


def result_node(**overrides):
    fields = dict(evidence_id="E", type=EvidenceType.STATISTICAL_RESULT,
        polarity=EvidencePolarity.NEUTRAL, summary="No tested feature met the declared FDR threshold.",
        source_task_id="T", source_artifact_uris=["table://study/deg/v1"],
        inference_contract_id="inference-v1:fixture", audit_passed=True,
        metrics={"scientific_status": "estimated_inconclusive"})
    fields.update(overrides)
    return EvidenceNode(**fields)


def test_inconclusive_summary_is_exact_observation_with_zero_discovery_confidence():
    graph = EvidenceGraph()
    node = result_node()
    graph.add_evidence(node)
    claim = ClaimEngine(graph).create_result_summary("C", "E")
    assert claim.statement == node.summary
    assert claim.language_tier == LanguageTier.LEVEL_1_OBSERVATION
    assert claim.claim_type == ClaimType.RESULT_SUMMARY
    assert claim.inference_contract_ids == [node.inference_contract_id]
    assert set(claim.confidence.model_dump().values()) == {0.0}


@pytest.mark.parametrize("statement", ["The treatment has no effect.", "The groups are equivalent."])
def test_neutral_result_cannot_support_positive_or_equivalence_claim(statement):
    graph = EvidenceGraph()
    graph.add_evidence(result_node())
    with pytest.raises(ValueError, match="Neutral"):
        ClaimEngine(graph).create_claim(claim_id="C", statement=statement,
            language_tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE,
            support_evidence_ids=["E"])


@pytest.mark.parametrize("overrides", [
    {"audit_passed": False}, {"inference_contract_id": None},
    {"source_artifact_uris": []},
    {"metrics": {"scientific_status": "assumptions_failed"}},
])
def test_unadmitted_result_cannot_become_result_claim(overrides):
    graph = EvidenceGraph()
    graph.add_evidence(result_node(**overrides))
    with pytest.raises(ValueError):
        ClaimEngine(graph).create_result_summary("C", "E")


def test_significant_claim_carries_inference_identity():
    graph = EvidenceGraph()
    node = result_node(type=EvidenceType.PSEUDOBULK_DEG,
        polarity=EvidencePolarity.SUPPORTING, summary="Gene A differs between conditions.",
        metrics={"fdr_q_value": .01})
    graph.add_evidence(node)
    claim = ClaimEngine(graph).create_claim(claim_id="C", statement=node.summary,
        language_tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE, support_evidence_ids=["E"])
    assert claim.inference_contract_ids == [node.inference_contract_id]


def test_report_preserves_unknown_assumptions_and_failed_admission(tmp_path):
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.report.markdown_report import ScientificReportGenerator
    from eacbp.schemas.study import StudyManifest, BiologicalDesign
    from eacbp.schemas.task import TaskResult, TaskStatus, ScientificResult, AssumptionAssessment
    result = TaskResult(task_id="T", status=TaskStatus.SCIENTIFIC_FAILURE,
        capability="deg", method_used="pydeseq2", scientific_result=ScientificResult(
            status="assumptions_failed", summary="Design is confounded.", assumptions=[
                AssumptionAssessment(name="donor_independence", status="unknown",
                    reason="Metadata cannot establish biological independence.")]))
    report = ScientificReportGenerator(
        StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="kidney")),
        EvidenceGraph(), ArtifactRegistry(str(tmp_path)), [result])
    text = report.generate_markdown()
    assert "assumptions_failed" in text
    assert "not admitted" in text
    assert "donor_independence | unknown" in text
    assert "Non-significance does not establish" in text


@pytest.mark.parametrize("mismatch", ["candidate", "artifact"])
def test_normalization_rejects_foreign_inference_identity(tmp_path, mismatch):
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.auditor.base import ValidationReport
    from eacbp.evidence.provenance import normalize_evidence_candidates
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    registry.register(uri, {"value": 1}, ArtifactType.JSON, "study", "T", "fixture",
        summary_metrics={"inference_contract_id": "foreign" if mismatch == "artifact" else "expected"})
    task = TaskContract(task_id="T", capability="deg")
    result = TaskResult(task_id="T", capability="deg", method_used="fixture",
        status=TaskStatus.SUCCESS, output_artifacts=[uri], inference_contract_id="expected")
    node = result_node(source_artifact_uris=[uri],
        inference_contract_id="foreign" if mismatch == "candidate" else "expected")
    report = ValidationReport(auditor_name="fixture", target_task_id="T")
    with pytest.raises(ValueError, match="inference contract differs"):
        normalize_evidence_candidates([node], task, result, report, registry)
