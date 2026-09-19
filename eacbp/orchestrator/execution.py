"""Execution of one task contract, isolated from orchestration policy."""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Optional

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.capabilities import CapabilityRegistry
from eacbp.schemas.task import ExecutionFailureType, TaskContract, TaskResult, TaskStatus


@dataclass
class TaskExecutionOutcome:
    """Result of an execution attempt sequence.

    ``staged`` is the private transaction belonging to the final attempt.  A
    caller commits it only after receiving a successful result; failed
    attempts remain private and are never published by this class.
    """

    result: TaskResult
    staged: TaskArtifactTransaction
    attempts: int
    execution_time_sec: float


class TaskExecutor:
    """Run a contract with the same retry/fallback policy as the legacy loop."""

    def __init__(self, artifact_registry: ArtifactRegistry, capability_registry: CapabilityRegistry):
        self.artifact_registry = artifact_registry
        self.capability_registry = capability_registry
        self.event_journal = None

    def _emit(self, kind, task, **details):
        if self.event_journal is not None:
            self.event_journal.emit(kind, task_id=task.task_id, **details)

    def _attempt(self, task, attempt):
        staged = TaskArtifactTransaction(self.artifact_registry)
        started = time.monotonic()
        self._emit("attempt_started", task, attempt=attempt, method=task.method)
        try:
            result = self.capability_registry.execute_contract(task, staged)
        except BaseException as exc:
            staged.close()
            self._emit("attempt_finished", task, attempt=attempt, status="exception",
                       error_type=type(exc).__name__, error=str(exc),
                       duration_seconds=max(0.0, time.monotonic() - started))
            raise
        self._emit("attempt_finished", task, attempt=attempt, method=result.method_used,
                   status=result.status.value,
                   error_type=result.error_type.value if result.error_type else None,
                   error=result.error_message,
                   duration_seconds=max(0.0, time.monotonic() - started))
        if result.status != TaskStatus.SUCCESS:
            staged.close()
        return staged, result

    def execute(self, task: TaskContract) -> TaskExecutionOutcome:
        started = time.time()
        # Each attempt receives a fresh private namespace.  This preserves the
        # old public ``execute_contract`` seam used by tests and plugin code.
        attempt_limit = min(
            task.retry_policy.max_execution_retry + 1,
            task.retry_policy.require_human_after,
        )
        staged: Optional[TaskArtifactTransaction] = None
        result: Optional[TaskResult] = None
        attempt = 0
        for attempt in range(attempt_limit):
            if attempt:
                self._emit("retry_scheduled", task, attempt=attempt + 1,
                           reason=result.error_message, retry_kind="execution")
            staged, result = self._attempt(task, attempt + 1)
            if result.status != TaskStatus.EXECUTION_FAILURE:
                break
            if result.error_type not in (ExecutionFailureType.CODE_ERROR,):
                break

        # A method-level failure may switch to the configured fallback methods.
        # Keep the exact attempt cap and ordering of the previous inline loop.
        if result is not None and result.status == TaskStatus.METHOD_FAILURE:
            for fallback in task.retry_policy.fallback_methods[:task.retry_policy.max_method_retry]:
                if attempt + 1 >= task.retry_policy.require_human_after:
                    break
                attempt += 1
                fallback_task = task.model_copy(update={"method": fallback})
                self._emit("retry_scheduled", fallback_task, attempt=attempt + 1,
                           reason=result.error_message, retry_kind="method_fallback")
                staged, result = self._attempt(fallback_task, attempt + 1)
                if result.status != TaskStatus.METHOD_FAILURE:
                    break

        # ``attempt_limit`` is always >= 1 under the schema, therefore both
        # values are initialized for static type checkers and defensive use by
        # callers constructing non-Pydantic contracts.
        if result is None or staged is None:
            raise RuntimeError(f"Task executor produced no result for {task.task_id}")
        result.metrics["execution_attempts"] = attempt + 1
        elapsed = round(time.time() - started, 3)
        result.execution_time_sec = elapsed
        return TaskExecutionOutcome(
            result=result,
            staged=staged,
            attempts=attempt + 1,
            execution_time_sec=elapsed,
        )
