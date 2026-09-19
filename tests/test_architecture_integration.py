"""Adversarial integration checks for the typed runtime and audit boundary.

These tests exercise the state and artifact contracts through their public
interfaces.  A published computation is intentionally usable for recovery,
but it is not evidence until a matching audit receipt has been persisted.
"""

from pathlib import Path

import pytest

from eacbp.artifact.registry import (
    ArtifactAuditError,
    ArtifactAuditAccessError,
    ArtifactRegistry,
)
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationCheck, ValidationReport, ValidationSeverity
from eacbp.capabilities.base import BaseCapability
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.runtime import ExecutionState, RunConfig
from eacbp.schemas.study import AnalysisPolicy, BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


class _FixtureCapability(BaseCapability):
    """Small real transaction producer for the orchestrator integration test."""

    def __init__(self):
        super().__init__(
            capability_name="fixture",
            implementation_id="fixture_v1",
            accepts_types=[],
            output_types=[ArtifactType.JSON],
        )

    def execute(self, contract, registry):
        uri = contract.expected_outputs[0]
        registry.register(
            uri,
            {"value": "committed"},
            ArtifactType.JSON,
            "study",
            contract.task_id,
            self.capability_name,
        )
        return TaskResult(
            task_id=contract.task_id,
            capability=contract.capability,
            method_used=self.implementation_id,
            status=TaskStatus.SUCCESS,
            output_artifacts=[uri],
            metrics={"observed_metric": 42},
        )


class _PassingFixtureAuditor:
    auditor_name = "fixture-loop-auditor"

    def audit_task(self, contract, result, registry):
        return ValidationReport(
            auditor_name=self.auditor_name,
            target_task_id=contract.task_id,
            target_artifact_uri=result.output_artifacts[0],
            checks=[
                ValidationCheck(
                    check_name="fixture_integrity",
                    passed=True,
                    severity=ValidationSeverity.INFO,
                    message="fixture payload is valid",
                )
            ],
            overall_passed=True,
        )


def _contract_and_result(*uris: str):
    contract = TaskContract(
        task_id="task-audit",
        capability="fixture",
        method="fixture_v1",
        input_artifacts=[],
        expected_outputs=list(uris),
        parameters={"context": "stable"},
    )
    result = TaskResult(
        task_id=contract.task_id,
        capability=contract.capability,
        method_used=contract.method or "fixture_v1",
        status=TaskStatus.SUCCESS,
        output_artifacts=list(uris),
        metrics={"observed": 1},
    )
    return contract, result


def _passed_report(contract: TaskContract, *uris: str) -> ValidationReport:
    return ValidationReport(
        auditor_name="fixture-auditor",
        target_task_id=contract.task_id,
        target_artifact_uri=uris[0] if uris else None,
        checks=[
            ValidationCheck(
                check_name="fixture_integrity",
                passed=True,
                severity=ValidationSeverity.INFO,
                message="fixture payload is valid",
            )
        ],
        overall_passed=True,
    )


def _failed_report(contract: TaskContract, *uris: str) -> ValidationReport:
    return ValidationReport(
        auditor_name="fixture-auditor",
        target_task_id=contract.task_id,
        target_artifact_uri=uris[0] if uris else None,
        checks=[
            ValidationCheck(
                check_name="fixture_integrity",
                passed=False,
                severity=ValidationSeverity.ERROR,
                message="fixture payload rejected",
            )
        ],
        overall_passed=False,
    )


def _publish_outputs(registry: ArtifactRegistry, signature: str, *uris: str) -> None:
    """Publish fixture outputs with the computation receipt audit requires."""

    tx = TaskArtifactTransaction(registry)
    for index, uri in enumerate(uris):
        tx.register(uri, {"value": index}, ArtifactType.JSON, "study", "task-audit", "fixture")
    result = TaskResult(
        task_id="task-audit",
        capability="fixture",
        method_used="fixture_v1",
        status=TaskStatus.SUCCESS,
        output_artifacts=list(uris),
    )
    registry.commit_task(tx, signature, result)


def _record_pass(registry: ArtifactRegistry, signature: str, contract: TaskContract, result: TaskResult):
    report = _passed_report(contract, *result.output_artifacts)
    registry.begin_audit(
        signature,
        contract,
        result,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )
    return registry.record_audit(
        signature,
        contract,
        result,
        report,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )


def test_metrics_cannot_inject_routing_configuration():
    """Observed metrics adapt planning, but never create routing options."""

    state = ExecutionState.from_input({"mode": "real"})
    state.record_metrics(
        {
            "method_overrides": {"deg": "untrusted_method"},
            "analysis_extensions": {"functional_activity": True},
            "advanced_analysis": True,
            "n_cells": 17,
        }
    )

    view = state.current_state
    assert view["mode"] == "real"
    assert view["n_cells"] == 17
    assert "method_overrides" not in view
    assert "analysis_extensions" not in view
    assert "advanced_analysis" not in view
    assert state.config.method_overrides is None
    assert state.config.analysis_extensions is None
    assert state.config.advanced_analysis is None


def test_explicit_config_wins_over_conflicting_observations():
    state = ExecutionState.from_input(
        {
            "method_profile": "baseline",
            "advanced_analysis": False,
            "analysis_extensions": {"functional_activity": False},
        }
    )
    state.record_metrics(
        {
            "method_profile": "standard",
            "advanced_analysis": True,
            "analysis_extensions": {"functional_activity": True},
            "n_genes": 9,
        }
    )

    assert state.current_state["method_profile"] == "baseline"
    assert state.current_state["advanced_analysis"] is False
    assert state.current_state["analysis_extensions"] == {"functional_activity": False}
    assert state.current_state["n_genes"] == 9


def test_legacy_observational_batches_are_refreshed_by_audit_metrics():
    """Legacy ``batches`` input is an observation, rather than frozen config."""

    state = ExecutionState.from_input({"batches": ["stale-before-audit"]})
    state.record_metrics({"batches": ["batch-a", "batch-b"]})

    assert state.current_state["batches"] == ["batch-a", "batch-b"]


def test_runconfig_and_mapping_inputs_have_identical_resume_semantics():
    mapping_state = ExecutionState.from_input({"resume": True, "mode": "demo"})
    config_state = ExecutionState.from_input(RunConfig(resume=True, mode="demo"))

    assert mapping_state.resume_requested is True
    assert config_state.resume_requested is True
    assert mapping_state.current_state == config_state.current_state


def test_orchestrator_does_not_admit_state_or_graph_when_audit_persistence_fails(
    tmp_path, monkeypatch
):
    """Staged evidence cannot enter state or the graph if receipt persistence fails."""

    uri = "json://study/orchestrated/v1"
    task = TaskContract(
        task_id="fixture-task",
        capability="fixture",
        method="fixture_v1",
        expected_outputs=[uri],
    )
    manifest = StudyManifest(
        study_id="fixture-study",
        biological_design=BiologicalDesign(species="mouse", tissue="cortex"),
        analysis_policy=AnalysisPolicy(strict_reproducibility=False),
    )
    capabilities = CapabilityRegistry()
    capabilities.register(_FixtureCapability())
    registry = ArtifactRegistry(str(tmp_path))
    auditor = _PassingFixtureAuditor()
    orchestrator = ScientificOrchestrator(
        artifact_registry=registry,
        capability_registry=capabilities,
        auditor=auditor,
    )
    monkeypatch.setattr(ComputationalDAGPlanner, "build_study_plan", lambda *_args: [task])
    monkeypatch.setattr(ComputationalDAGPlanner, "order_tasks", lambda tasks: tasks)

    def failing_audit(*_args, **_kwargs):
        raise ArtifactAuditError("simulated audit receipt write failure")

    registry.record_audit = failing_audit

    def extracting_evidence(contract, result, _report):
        return [
            EvidenceNode(
                evidence_id="fixture-evidence",
                type=EvidenceType.QC_METRICS,
                summary="staged fixture evidence",
                source_task_id=contract.task_id,
                source_artifact_uris=list(result.output_artifacts),
                audit_passed=True,
            )
        ]

    orchestrator.extract_evidence_from_result = extracting_evidence

    outcome = orchestrator.run_study(manifest)

    assert outcome["status"] == "failed"
    assert orchestrator.evidence_graph.evidence_nodes == {}
    assert "observed_metric" not in orchestrator.execution_state.current_state
    signature = next(iter(registry.task_commits))
    assert registry.audit_records[signature].status.value == "pending"
    with pytest.raises(ArtifactAuditAccessError):
        registry.query_audited(uri, signature=signature, contract=task)


def test_computed_checkpoint_is_not_audited_after_restart(tmp_path):
    """A committed computation receipt remains unaudited after a process restart."""

    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    tx = TaskArtifactTransaction(registry)
    tx.register(uri, {"value": 1}, ArtifactType.JSON, "study", "task-audit", "fixture")
    contract, result = _contract_and_result(uri)
    signature = "computed-only-signature"

    registry.commit_task(tx, signature, result)
    assert registry.get_task_commit(signature)["phase"] == "computed"

    restarted = ArtifactRegistry(str(tmp_path))
    assert restarted.get(uri)[1] == {"value": 1}
    with pytest.raises(ArtifactAuditAccessError):
        restarted.query_audited(uri, signature=signature, contract=contract)


def test_audit_receipt_survives_restart_but_context_mismatch_fails_closed(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    contract, result = _contract_and_result(uri)
    signature = "stable-signature"
    _publish_outputs(registry, signature, uri)
    _record_pass(registry, signature, contract, result)

    restarted = ArtifactRegistry(str(tmp_path))
    admitted = restarted.query_audited(uri, signature=signature, contract=contract)
    assert admitted[1] == {"value": 0}

    changed_contract = contract.model_copy(update={"parameters": {"context": "changed"}})
    with pytest.raises(ArtifactAuditAccessError):
        restarted.query_audited(uri, signature=signature, contract=changed_contract)
    with pytest.raises(ArtifactAuditAccessError):
        restarted.query_audited(uri, signature="different-signature", contract=contract)


def test_sibling_output_tampering_revokes_audited_access(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    uri_a = "json://study/result_a/v1"
    uri_b = "json://study/result_b/v1"
    contract, result = _contract_and_result(uri_a, uri_b)
    signature = "two-output-signature"
    _publish_outputs(registry, signature, uri_a, uri_b)
    _record_pass(registry, signature, contract, result)

    assert registry.query_audited(uri_a, signature=signature, contract=contract)[1] == {"value": 0}
    audited_many = registry.query_audited_many(signature=signature, contract=contract)
    assert {metadata.uri for metadata, _payload in audited_many} == {uri_a, uri_b}
    sibling_path = Path(registry.get_metadata(uri_b).storage_path)
    sibling_path.write_bytes(sibling_path.read_bytes() + b"tampered")

    with pytest.raises(ArtifactAuditAccessError):
        registry.query_audited(uri_a, signature=signature, contract=contract)
    with pytest.raises(ArtifactAuditAccessError):
        registry.query_audited_many(signature=signature, contract=contract)


def test_audit_report_cannot_bind_a_receipt_to_another_task(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    contract, result = _contract_and_result(uri)
    signature = "wrong-task-report-signature"
    _publish_outputs(registry, signature, uri)
    registry.begin_audit(
        signature,
        contract,
        result,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )
    wrong_task = _passed_report(contract, uri).model_copy(update={"target_task_id": "other-task"})
    record = registry.record_audit(
        signature,
        contract,
        result,
        wrong_task,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )

    assert record.status.value == "rejected"
    with pytest.raises(ArtifactAuditAccessError):
        registry.query_audited(uri, signature=signature, contract=contract)


def test_contradictory_failed_error_check_cannot_be_marked_passed(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    contract, result = _contract_and_result(uri)
    signature = "contradictory-report-signature"
    _publish_outputs(registry, signature, uri)
    registry.begin_audit(
        signature,
        contract,
        result,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )
    contradictory = ValidationReport(
        auditor_name="fixture-auditor",
        target_task_id=contract.task_id,
        target_artifact_uri=uri,
        checks=[
            ValidationCheck(
                check_name="failed_error_check",
                passed=False,
                severity=ValidationSeverity.ERROR,
                message="the independent check failed",
            )
        ],
        overall_passed=True,
    )
    record = registry.record_audit(
        signature,
        contract,
        result,
        contradictory,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )

    assert record.status.value == "rejected"
    with pytest.raises(ArtifactAuditAccessError):
        registry.query_audited(uri, signature=signature, contract=contract)


def test_failed_or_crashed_reaudit_revokes_previous_pass_after_restart(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    uri = "json://study/result/v1"
    contract, result = _contract_and_result(uri)
    signature = "re-audit-signature"
    _publish_outputs(registry, signature, uri)
    _record_pass(registry, signature, contract, result)
    assert registry.query_audited(uri, signature=signature, contract=contract)[1] == {"value": 0}

    # Starting a re-audit must revoke the previous pass before the auditor is
    # called.  Simulate a crash by leaving this pending record on disk.
    registry.begin_audit(
        signature,
        contract,
        result,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )
    crashed_restart = ArtifactRegistry(str(tmp_path))
    with pytest.raises(ArtifactAuditAccessError):
        crashed_restart.query_audited(uri, signature=signature, contract=contract)

    # A completed failed re-audit remains rejected and cannot resurrect the
    # old pass even in a fresh registry instance.
    failed = _failed_report(contract, uri)
    crashed_restart.record_audit(
        signature,
        contract,
        result,
        failed,
        auditor_version="fixture-1",
        auditor_fingerprint="fixture-fingerprint",
    )
    restarted = ArtifactRegistry(str(tmp_path))
    assert restarted.audit_records[signature].status.value == "rejected"
    with pytest.raises(ArtifactAuditAccessError):
        restarted.query_audited(uri, signature=signature, contract=contract)
