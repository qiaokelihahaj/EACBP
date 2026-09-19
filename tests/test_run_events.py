import json
from unittest.mock import patch

import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.events import RunEventJournal, read_run_events
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest, BiologicalDesign, AnalysisPolicy
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType


def test_event_log_restart_separation_and_crash_tail(tmp_path):
    first = RunEventJournal(tmp_path, "../../study")
    first.emit("run_started")
    first.emit("task_started", task_id="a")
    with first.path.open("ab") as stream:
        stream.write(b'{"truncated":')
    assert [event.kind for event in read_run_events(first.path)] == ["run_started", "task_started"]
    second = RunEventJournal(tmp_path, "../../study")
    second.emit("run_started")
    assert second.path != first.path
    assert second.path.resolve().is_relative_to(tmp_path.resolve())
    assert read_run_events(second.path)[0].sequence == 1
    with second.path.open("a") as stream:
        stream.write("not-json\n")
    with pytest.raises(ValueError):
        read_run_events(second.path)


def test_event_write_failure_is_visible_and_does_not_fail_science(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    orchestrator = ScientificOrchestrator(registry)
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="test"))
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[]), \
         patch("eacbp.orchestrator.events.os.fsync", side_effect=OSError("event disk full")):
        outcome = orchestrator.run_study(manifest)
    assert outcome["status"] == "success"
    assert "event disk full" in outcome["observability_warnings"][0]


def test_retry_resume_and_audit_rejection_are_observable(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    orch = ScientificOrchestrator(reg)
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="test"),
                             analysis_policy=AnalysisPolicy(strict_reproducibility=False))
    task = TaskContract(task_id="a", capability="dataset_audit", expected_outputs=["json://s/a/v1"])
    calls = []

    def execute(contract, staged):
        calls.append(contract.task_id)
        if len(calls) == 1:
            return TaskResult(task_id="a", capability="dataset_audit", method_used="sc_audit_v1",
                              status=TaskStatus.EXECUTION_FAILURE,
                              error_type=ExecutionFailureType.CODE_ERROR, error_message="temporary failure")
        staged.register("json://s/a/v1", {"value": 1}, ArtifactType.JSON, "s", "a", "fixture")
        return TaskResult(task_id="a", capability="dataset_audit", method_used="sc_audit_v1",
                          status=TaskStatus.SUCCESS, output_artifacts=["json://s/a/v1"])

    with patch.object(ComputationalDAGPlanner, "build_study_plan", side_effect=lambda *_: [task.model_copy(deep=True)]), \
         patch.object(orch.capability_registry, "execute_contract", side_effect=execute), \
         patch.object(orch.auditor, "audit_task", return_value=ValidationReport(auditor_name="test", target_task_id="a")):
        first = orch.run_study(manifest)
        assert first["status"] == "success", first
        resumed = orch.run_study(manifest, {"resume": True})
    assert resumed["status"] == "success", resumed
    assert len(calls) == 2
    initial_events = read_run_events(first["event_log"])
    retry = next(e for e in initial_events if e.kind == "retry_scheduled")
    assert retry.details["reason"] == "temporary failure"
    resumed_events = read_run_events(resumed["event_log"])
    assert "task_resumed" in [e.kind for e in resumed_events]
    assert "attempt_started" not in [e.kind for e in resumed_events]
    assert "audit_started" in [e.kind for e in resumed_events]
    assert first["run_id"] != resumed["run_id"]

    with patch.object(ComputationalDAGPlanner, "build_study_plan", side_effect=lambda *_: [task.model_copy(deep=True)]), \
         patch.object(orch.auditor, "audit_task", return_value=ValidationReport(
             auditor_name="test", target_task_id="a", overall_passed=False, stop_rule_triggered=True)):
        rejected = orch.run_study(manifest, {"resume": True})
    assert rejected["status"] == "failed"
    assert any(e.kind == "audit_rejected" for e in read_run_events(rejected["event_log"]))
    assert rejected["evidence_nodes_count"] == 0


def test_planning_failure_has_terminal_event(tmp_path):
    orch = ScientificOrchestrator(ArtifactRegistry(str(tmp_path)))
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="test"))
    with pytest.raises(ValueError):
        orch.run_study(manifest, {"method_profile": "invalid"})
    events = read_run_events(orch.event_journal.path)
    assert events[-1].kind == "run_failed"
