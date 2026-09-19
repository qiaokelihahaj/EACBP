"""Per-invocation execution events, independent of scientific checkpoints."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
from typing import Any, Literal, Optional
from uuid import uuid4

from pydantic import BaseModel, Field

from eacbp.orchestrator.checkpoint import fingerprint


EventKind = Literal[
    "run_started", "plan_created", "task_started", "method_resolved",
    "task_resumed", "attempt_started", "attempt_finished", "retry_scheduled",
    "artifacts_committed", "audit_started", "audit_finished", "audit_rejected",
    "task_blocked", "task_finished", "task_failed", "plan_adapted",
    "run_finished", "run_failed",
]


class RunEvent(BaseModel):
    schema_version: int = 1
    run_id: str
    study_id: str
    sequence: int = Field(ge=1)
    timestamp: datetime
    elapsed_seconds: float = Field(ge=0)
    kind: EventKind
    task_id: Optional[str] = None
    details: dict[str, Any] = Field(default_factory=dict)


class RunEventJournal:
    """Write one JSONL file per invocation, including resumed invocations.

    Events describe execution; only the task journal and audit receipts govern
    recovery/admission. Diagnostic I/O failure is exposed to the caller without
    changing a scientific result. Each successful event is flushed to disk.
    """

    def __init__(self, storage_dir, study_id: str):
        self.run_id = uuid4().hex
        self.study_id = study_id
        self.path = Path(storage_dir) / "_runs" / "events" / fingerprint(study_id) / f"{self.run_id}.jsonl"
        self.warnings: list[str] = []
        self._started = time.monotonic()
        self._sequence = 0
        self._disabled = False
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            # Exclusive creation prevents accidental reuse of a prior log.
            self.path.open("x", encoding="utf-8").close()
        except OSError as exc:
            self._disable(exc)

    def _disable(self, exc: Exception) -> None:
        self._disabled = True
        self.warnings.append(f"Execution event log unavailable: {type(exc).__name__}: {exc}")

    def emit(self, kind: EventKind, task_id: Optional[str] = None, **details) -> None:
        if self._disabled:
            return
        try:
            event = RunEvent(
                run_id=self.run_id, study_id=self.study_id,
                sequence=self._sequence + 1,
                timestamp=datetime.now(timezone.utc),
                elapsed_seconds=max(0.0, time.monotonic() - self._started),
                kind=kind, task_id=task_id, details=details,
            )
            line = json.dumps(event.model_dump(mode="json"), ensure_ascii=False, allow_nan=False)
            with self.path.open("a", encoding="utf-8") as stream:
                stream.write(line + "\n")
                stream.flush()
                os.fsync(stream.fileno())
            self._sequence += 1
        except (OSError, ValueError, TypeError) as exc:
            self._disable(exc)


def read_run_events(path) -> list[RunEvent]:
    """Read a log, tolerating only a crash-truncated final line.

    A malformed complete event is an error, rather than a silently omitted
    failure. Sequence/run identity checks prevent accidental log concatenation.
    """
    events = []
    data = Path(path).read_bytes()
    for index, line in enumerate(data.splitlines(keepends=True)):
        if not line.endswith(b"\n"):
            break
        event = RunEvent.model_validate_json(line)
        if event.sequence != index + 1:
            raise ValueError("Run event sequence is not contiguous")
        if events and (event.run_id, event.study_id) != (events[0].run_id, events[0].study_id):
            raise ValueError("Run event log contains different invocation identities")
        events.append(event)
    return events
