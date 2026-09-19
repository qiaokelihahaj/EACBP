from unittest.mock import patch
import pandas as pd
import pytest
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport
from eacbp.auditor.statistical import StatisticalValidator
from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.graph import EvidenceGraph
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType, LanguageTier
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def test_environment_fingerprint_records_code_and_dependencies():
    from eacbp.orchestrator.checkpoint import execution_environment
    environment = execution_environment()
    assert len(environment["source_sha256"]) == 64
    assert environment["python"]
    assert any(key.lower() == "numpy" for key in environment["distributions"])


def test_unrelated_sentence_cannot_borrow_significant_evidence():
    graph = EvidenceGraph()
    node = EvidenceNode(evidence_id="E", type=EvidenceType.PSEUDOBULK_DEG,
        summary="GeneA has higher expression in treated versus control.", source_task_id="T",
        audit_passed=True, source_artifact_uris=["table://s/deg/v1"], metrics={"fdr_q_value": .01},
        biological_context={"gene": "GeneA"})
    graph.add_evidence(node)
    with pytest.raises(ValueError, match="restate"):
        ClaimEngine(graph).create_claim(claim_id="C", statement="GeneB has higher expression in disease.",
            language_tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE, support_evidence_ids=["E"])


def test_verified_but_unrelated_knowledge_cannot_upgrade_interpretation():
    graph = EvidenceGraph()
    graph.add_evidence(EvidenceNode(evidence_id="E", type=EvidenceType.PSEUDOBULK_DEG,
        summary="GeneA has higher expression.", source_task_id="T", audit_passed=True,
        source_artifact_uris=["table://s/deg/v1"], metrics={"fdr_q_value": .01}, biological_context={"gene": "GeneA"}))
    graph.add_evidence(EvidenceNode(evidence_id="K", type=EvidenceType.LITERATURE_SUPPORT,
        summary="GeneB reference", source_task_id="K", audit_passed=True, source_verified=True,
        source_artifact_uris=["json://s/reference/v1"], biological_context={"gene": "GeneB"}))
    with pytest.raises(ValueError, match="biological context"):
        ClaimEngine(graph).create_claim(claim_id="C", statement="GeneA has higher expression.",
            language_tier=LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION, support_evidence_ids=["E", "K"])


def manifest(sid="review"):
    return StudyManifest(study_id=sid, biological_design=BiologicalDesign(species="human", tissue="kidney"))


def deg(reg, q=.9):
    uri = "table://review/deg/v1"
    reg.register(uri, pd.DataFrame([{"gene": "KidneyGene", "log2_fold_change": .1,
                                   "p_value": .8, "fdr_q_value": q}]),
                 ArtifactType.TABLE, "review", "deg", "deg",
                 parameters={"cond_ad": "injury", "cond_ctrl": "healthy"},
                 summary_metrics={"statistical_unit": "donor_pseudobulk"})
    task = TaskContract(task_id="deg", capability="deg")
    result = TaskResult(task_id="deg", capability="deg", method_used="donor_mean_welch",
                        status=TaskStatus.SUCCESS, output_artifacts=[uri], metrics={"is_pseudobulk": True})
    return task, result


def test_nonsignificant_deg_does_not_become_evidence(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    task, result = deg(reg)
    orch = ScientificOrchestrator(artifact_registry=reg)
    audit = StatisticalValidator().audit(task, result, reg)
    assert audit.checks[0].passed  # shared donor_pseudobulk vocabulary
    assert orch.extract_evidence_from_result(task, result, audit) == []


def test_significant_gene_preserves_context_and_direction(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    task, result = deg(reg, q=.01)
    orch = ScientificOrchestrator(artifact_registry=reg)
    nodes = orch.extract_evidence_from_result(task, result, StatisticalValidator().audit(task, result, reg))
    assert len(nodes) == 1
    assert "KidneyGene" in nodes[0].summary and "injury versus healthy" in nodes[0].summary
    assert "AD microglia" not in nodes[0].summary
    for node in nodes:
        orch.evidence_graph.add_evidence(node)
    orch._synthesize_study_claims(manifest())
    claim = next(iter(orch.evidence_graph.claim_nodes.values()))
    assert claim.language_tier == LanguageTier.LEVEL_2_STATISTICAL_INFERENCE
    assert "Trem2" not in claim.statement


def test_claim_rejects_empty_unknown_and_unaudited_support():
    graph = EvidenceGraph()
    engine = ClaimEngine(graph)
    args = dict(claim_id="C", statement="Expression differs.", language_tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE)
    with pytest.raises(ValueError, match="at least one"):
        engine.create_claim(**args)
    with pytest.raises(ValueError, match="Unknown"):
        engine.create_claim(**args, support_evidence_ids=["missing"])
    graph.add_evidence(EvidenceNode(evidence_id="E", type=EvidenceType.PSEUDOBULK_DEG,
        summary="Unaudited test", source_task_id="T", metrics={"fdr_q_value": .001}))
    with pytest.raises(ValueError, match="passed audit"):
        engine.create_claim(**args, support_evidence_ids=["E"])


@pytest.mark.parametrize("q", [.9, float("nan"), -1, None])
def test_statistics_gate_requires_valid_adjusted_p(q):
    node = EvidenceNode(evidence_id="E", type=EvidenceType.PSEUDOBULK_DEG,
        summary="test", source_task_id="T", source_artifact_uris=["table://s/t/v1"],
        audit_passed=True, metrics={"fdr_q_value": q})
    assert not ClaimEngine.has_valid_statistics(node)


def test_stop_rule_blocks_dependents_and_excludes_evidence(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    orch = ScientificOrchestrator(artifact_registry=reg)
    tasks = [TaskContract(task_id="a", capability="dataset_audit"),
             TaskContract(task_id="b", capability="dataset_audit", depends_on=["a"])]
    result = TaskResult(task_id="a", capability="dataset_audit", method_used="sc_audit_v1", status=TaskStatus.SUCCESS)
    audit = ValidationReport(auditor_name="test", target_task_id="a", overall_passed=False, stop_rule_triggered=True)
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=tasks), \
         patch.object(orch.capability_registry, "execute_contract", return_value=result) as execute, \
         patch.object(orch.auditor, "audit_task", return_value=audit):
        summary = orch.run_study(manifest())
    assert execute.call_count == 1
    assert summary["status"] == "failed" and summary["tasks_blocked"] == 1
    assert summary["claims_count"] == summary["evidence_nodes_count"] == 0
    assert orch.task_history[0].status == TaskStatus.SCIENTIFIC_FAILURE


def test_resume_validates_outputs_and_does_not_reexecute(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    orch = ScientificOrchestrator(artifact_registry=reg)
    uri = "json://review/audit/v1"
    task = TaskContract(task_id="audit", capability="dataset_audit", expected_outputs=[uri])
    def execute(t, registry):
        registry.register(uri, {"ok": True}, ArtifactType.JSON, "review", "audit", "audit")
        return TaskResult(task_id=t.task_id, capability=t.capability, method_used="sc_audit_v1",
            status=TaskStatus.SUCCESS, output_artifacts=[uri], metrics={"n_cells": 2, "min_replicates": 1})
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[task]), \
         patch.object(orch.capability_registry, "execute_contract", side_effect=execute) as call:
        first = orch.run_study(manifest())
        second = orch.run_study(manifest(), {"resume": True})
        assert first["status"] == second["status"] == "success"
        assert call.call_count == 1
        assert second["evidence_nodes_count"] == first["evidence_nodes_count"]
        with pytest.raises(ValueError, match="resume"):
            orch.run_study(manifest())


def test_state_and_evidence_are_reset_between_studies(tmp_path):
    orch = ScientificOrchestrator(artifact_registry=ArtifactRegistry(str(tmp_path)))
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[]):
        orch.run_study(manifest("one"), {"target_gene": "OldGene"})
        result = orch.run_study(manifest("two"))
    assert "target_gene" not in orch.current_state
    assert result["claims_count"] == 0


def test_task_dependency_validation_and_order():
    child = TaskContract(task_id="b", capability="qc", depends_on=["a"])
    parent = TaskContract(task_id="a", capability="dataset_audit")
    assert [t.task_id for t in ComputationalDAGPlanner.order_tasks([child, parent])] == ["a", "b"]
    parent.depends_on = ["b"]
    with pytest.raises(ValueError, match="cycle"):
        ComputationalDAGPlanner.order_tasks([parent, child])


def test_resume_rejects_changed_input(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    orch = ScientificOrchestrator(artifact_registry=reg)
    uri = "json://review/input/v1"
    meta = reg.register(uri, {"value": 1}, ArtifactType.JSON, "review", "ingest", "ingest")
    task = TaskContract(task_id="audit", capability="dataset_audit", input_artifacts=[uri])
    result = TaskResult(task_id="audit", capability="dataset_audit", method_used="sc_audit_v1",
                        status=TaskStatus.SUCCESS, output_artifacts=[uri])
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[task]), \
         patch.object(orch.capability_registry, "execute_contract", return_value=result) as call:
        assert orch.run_study(manifest())["status"] == "success"
        from pathlib import Path
        Path(meta.storage_path).write_text('{"value": 999}', encoding="utf-8")
        resumed = orch.run_study(manifest(), {"resume": True})
        assert resumed["status"] == "failed"
        assert resumed["claims_count"] == 0
        assert call.call_count == 1


def test_auditor_checks_secondary_output(tmp_path):
    import numpy as np
    from eacbp.auditor import ScientificAuditor
    reg = ArtifactRegistry(str(tmp_path))
    good = "json://review/primary/v1"
    bad = "table://review/secondary/v1"
    reg.register(good, {"value": 1}, ArtifactType.JSON, "review", "compute", "compute")
    reg.register(bad, pd.DataFrame({"value": [np.inf]}), ArtifactType.TABLE, "review", "compute", "compute")
    result = TaskResult(task_id="compute", capability="other", method_used="test", status=TaskStatus.SUCCESS, output_artifacts=[good, bad])
    audit = ScientificAuditor().audit_task(TaskContract(task_id="compute", capability="other"), result, reg)
    assert not audit.overall_passed


def test_journal_excludes_concurrent_run_and_releases_lock(tmp_path):
    from eacbp.orchestrator.checkpoint import StudyJournal
    with StudyJournal(tmp_path, "s"):
        with pytest.raises(RuntimeError, match="already running"):
            with StudyJournal(tmp_path, "s"):
                pass
    with StudyJournal(tmp_path, "s"):
        pass
