"""Resource policy must prevent work and interrupted publication, not just log it."""
import json
import subprocess
import sys
import threading
import time
from unittest.mock import Mock

import pytest

from eacbp.execution_context import (
    ExecutionContext, ExecutionCancelled, ExecutionControlError,
    ExecutionDeadlineExceeded, activate_execution, run_external,
)
from eacbp.schemas.study import Constraints
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.orchestrator.execution import TaskExecutor


def test_external_timeout_terminates_before_side_effect(tmp_path):
    marker = tmp_path / "too_late.txt"
    command = [sys.executable, "-c",
               "import time,pathlib;time.sleep(2);pathlib.Path(" + repr(str(marker)) + ").write_text('bad')"]
    with activate_execution(ExecutionContext()):
        with pytest.raises(ExecutionDeadlineExceeded):
            run_external(command, timeout=0.1, capture_output=True, text=True)
    assert not marker.exists()


def test_external_cancel_interrupts_active_command():
    context = ExecutionContext()
    timer = threading.Timer(0.3, context.cancel)
    timer.start()
    try:
        with activate_execution(context), pytest.raises(ExecutionCancelled):
            run_external([sys.executable, "-c", "import time;time.sleep(30)"], capture_output=True)
    finally:
        timer.cancel()


def test_external_environment_enforces_cpu_and_threads():
    with activate_execution(ExecutionContext(gpu_allowed=False, threads=2)):
        result = run_external([sys.executable, "-c", "import os,json;print(json.dumps([os.getenv('CUDA_VISIBLE_DEVICES'),os.getenv('OMP_NUM_THREADS')]))"],
                              capture_output=True, text=True, check=True)
    assert json.loads(result.stdout) == ["-1", "2"]


@pytest.mark.parametrize("parameters", [{"use_gpu": True}, {"extra_args": ["--cuda"]}, {"nested": {"device": "cuda:0"}}])
def test_forbidden_gpu_prevents_capability_execution(tmp_path, parameters):
    class MustNotRun:
        def execute_contract(self, *args):
            pytest.fail("GPU request must be rejected before execution")
    executor = TaskExecutor(ArtifactRegistry(str(tmp_path / "artifacts")), MustNotRun())
    executor.execution_context = ExecutionContext(gpu_allowed=False)
    result = executor.execute(TaskContract(task_id="t", capability="fake", parameters=parameters)).result
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.RESOURCE_POLICY
    assert result.metrics["execution_attempts"] == 1


def test_late_python_result_is_not_admitted(tmp_path):
    context = ExecutionContext()
    class LateResult:
        def execute_contract(self, task, staged):
            context.deadline = time.monotonic() - 1
            return TaskResult(task_id=task.task_id, capability=task.capability,
                              method_used="fake", status=TaskStatus.SUCCESS)
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    executor = TaskExecutor(registry, LateResult())
    executor.execution_context = context
    outcome = executor.execute(TaskContract(task_id="t", capability="fake"))
    assert outcome.result.error_type == ExecutionFailureType.TIMEOUT
    assert outcome.attempts == 1
    assert registry.list_artifacts() == []


def test_swallowed_external_timeout_still_fails_task(tmp_path):
    class SwallowsError:
        def execute_contract(self, task, staged):
            try:
                run_external([sys.executable, "-c", "import time;time.sleep(10)"], timeout=0.05)
            except Exception:
                return TaskResult(task_id=task.task_id, capability=task.capability,
                                  method_used="fake", status=TaskStatus.SUCCESS)
    executor = TaskExecutor(ArtifactRegistry(str(tmp_path / "artifacts")), SwallowsError())
    executor.execution_context = ExecutionContext()
    result = executor.execute(TaskContract(task_id="t", capability="fake")).result
    assert result.error_type == ExecutionFailureType.TIMEOUT


@pytest.mark.parametrize("hours", [0, -1, float("inf"), float("nan")])
def test_invalid_study_budget_rejected(hours):
    with pytest.raises(ValueError):
        Constraints(max_runtime_hours=hours)


def test_external_nonzero_preserves_diagnostics():
    with activate_execution(ExecutionContext()), pytest.raises(subprocess.CalledProcessError) as caught:
        run_external([sys.executable, "-c", "import sys;print('reason',file=sys.stderr);sys.exit(3)"],
                     capture_output=True, text=True, check=True)
    assert caught.value.returncode == 3
    assert "reason" in caught.value.stderr


def _resumable_study(tmp_path, monkeypatch):
    from eacbp.orchestrator.dag import ComputationalDAGPlanner
    from eacbp.orchestrator.loop import ScientificOrchestrator
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.study import BiologicalDesign, StudyManifest

    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    source, output = "json://control/input/v1", "json://control/audit/v1"
    registry.register(source, {"value": 1}, ArtifactType.JSON, "control", "ingest", "ingest")
    task = TaskContract(task_id="audit", capability="dataset_audit",
                        input_artifacts=[source], expected_outputs=[output])
    manifest = StudyManifest(study_id="control",
                             biological_design=BiologicalDesign(species="human", tissue="test"),
                             analysis_policy={"strict_reproducibility": False})
    orchestrator = ScientificOrchestrator(registry)

    def execute(contract, staged):
        staged.register(output, {"value": 1}, ArtifactType.JSON, "control", contract.task_id,
                        "audit", parent_uris=[source])
        return TaskResult(task_id=contract.task_id, capability=contract.capability,
                          method_used="sc_audit_v1", status=TaskStatus.SUCCESS,
                          input_artifacts=[source], output_artifacts=[output])

    execute_spy = Mock(side_effect=execute)
    audit_spy = Mock(wraps=orchestrator.auditor.audit_task)
    lookup_spy = Mock(wraps=orchestrator.resume_manager.lookup)
    monkeypatch.setattr(orchestrator.capability_registry, "execute_contract", execute_spy)
    monkeypatch.setattr(orchestrator.auditor, "audit_task", audit_spy)
    monkeypatch.setattr(orchestrator.resume_manager, "lookup", lookup_spy)
    monkeypatch.setattr(ComputationalDAGPlanner, "build_study_plan",
                        lambda *_args, **_kwargs: [task.model_copy(deep=True)])
    return orchestrator, manifest, execute_spy, audit_spy, lookup_spy


@pytest.mark.parametrize("interruption", ["cancelled", "timeout"])
def test_interrupted_resume_skips_execution_lookup_and_audit(tmp_path, monkeypatch, interruption):
    orchestrator, manifest, execute_spy, audit_spy, lookup_spy = _resumable_study(tmp_path, monkeypatch)
    assert orchestrator.run_study(manifest)["status"] == "success"
    receipts = orchestrator.artifact_registry.list_audit_records()
    context = ExecutionContext(deadline=-1) if interruption == "timeout" else ExecutionContext()
    if interruption == "cancelled":
        context.cancel()

    interrupted = orchestrator.run_study(manifest, {"resume": True}, execution_context=context)

    assert interrupted["status"] == "failed"
    assert interrupted["claims_count"] == 0
    assert orchestrator.task_history[0].error_type == ExecutionFailureType(interruption)
    assert execute_spy.call_count == audit_spy.call_count == lookup_spy.call_count == 1
    assert orchestrator.artifact_registry.list_audit_records() == receipts
    # A cancelled invocation must not overwrite the previous recovery state.
    assert orchestrator.run_study(manifest, {"resume": True})["status"] == "success"
    assert execute_spy.call_count == 1
    assert audit_spy.call_count == 2


@pytest.mark.parametrize("phase", ["audit", "evidence"])
def test_cancel_during_admission_retains_pending_receipt_for_resume(tmp_path, monkeypatch, phase):
    from eacbp.orchestrator.checkpoint import StudyJournal

    orchestrator, manifest, execute_spy, audit_spy, _ = _resumable_study(tmp_path, monkeypatch)
    context = ExecutionContext()
    target = orchestrator.auditor if phase == "audit" else orchestrator
    method_name = "audit_task" if phase == "audit" else "extract_evidence_from_result"
    original = getattr(target, method_name)

    def cancel_after_work(*args, **kwargs):
        result = original(*args, **kwargs)
        context.cancel()
        return result

    monkeypatch.setattr(target, method_name, cancel_after_work)
    interrupted = orchestrator.run_study(manifest, execution_context=context)

    assert interrupted["status"] == "failed"
    assert orchestrator.task_history[0].error_type == ExecutionFailureType.CANCELLED
    assert interrupted["claims_count"] == interrupted["evidence_nodes_count"] == 0
    assert orchestrator.audit_reports == []
    assert execute_spy.call_count == audit_spy.call_count == 1
    receipts = orchestrator.artifact_registry.list_audit_records()
    assert len(receipts) == 1 and receipts[0].status.value == "pending"
    with StudyJournal(orchestrator.artifact_registry.storage.base_dir, manifest.study_id) as journal:
        assert journal.entries["audit"]["phase"] == "computed"
        assert journal.entries["audit"]["result"]["status"] == "success"

    monkeypatch.setattr(target, method_name, original)
    recovered = orchestrator.run_study(manifest, {"resume": True})
    assert recovered["status"] == "success"
    assert recovered["claims_count"] > 0
    assert execute_spy.call_count == 1
    assert audit_spy.call_count == 2
    assert orchestrator.artifact_registry.list_audit_records()[0].status.value == "passed"
