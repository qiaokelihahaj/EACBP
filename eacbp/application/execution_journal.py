"""Durable invocation diagnostics for frozen clustering sensitivity runs.

The journal records what an executor attempted. It does not admit scientific
artifacts: a resume caller must re-read native outputs through its collector.
Published events are immutable local records with a consistency hash chain,
not independently signed attestations. Elapsed times never cross processes.
"""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time
from typing import Any, Mapping
from uuid import uuid4

from eacbp._atomic_json import atomic_write_json


_DIGEST = re.compile(r"[0-9a-f]{64}")
_IDENTIFIER = re.compile(r"[a-z][a-z0-9_]{0,79}")
_UUID = re.compile(r"[0-9a-f]{32}")
_STATUSES = {"completed", "failed", "interrupted", "postprocessing_failed"}


def _encode(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False).encode("utf-8")


def _fingerprint(value: Any) -> str:
    return hashlib.sha256(_encode(value)).hexdigest()


def _timestamp() -> str:
    return datetime.now(timezone.utc).isoformat()


def _validate_timestamp(value: Any) -> None:
    if not isinstance(value, str):
        raise ValueError("Journal timestamp is not a UTC string")
    timestamp = datetime.fromisoformat(value)
    if timestamp.tzinfo is None or timestamp.utcoffset().total_seconds() != 0:
        raise ValueError("Journal timestamp must retain its UTC offset")


def _validate_seconds(value: Any) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
        raise ValueError("Journal monotonic time must be a finite nonnegative number")


def _validate_identity(identity: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(dict(identity))
    for field in ("bindings_sha256", "plan_file_sha256", "plan_sha256", "source_sha256"):
        if not isinstance(result.get(field), str) or not _DIGEST.fullmatch(result[field]):
            raise ValueError(f"Journal requires its frozen {field}")
    environment = result.get("execution_environment")
    if (not isinstance(environment, dict) or not isinstance(environment.get("source_sha256"), str)
            or not _DIGEST.fullmatch(environment["source_sha256"])):
        raise ValueError("Journal requires its frozen execution environment and package source hash")
    _encode(result)
    return result


def _publish_exclusive(path: Path, value: Any) -> None:
    """Flush a complete same-filesystem document before exclusive publication."""
    data = _encode(value) + b"\n"
    descriptor, name = tempfile.mkstemp(prefix=".journal-", suffix=".tmp", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        # A hard link publishes an already complete file and cannot replace an
        # old event. The executor lock serializes sequence selection, not data.
        os.link(temporary, path)
        if os.name != "nt":
            directory = os.open(path.parent, os.O_RDONLY)
            try:
                os.fsync(directory)
            finally:
                os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


class ExecutionJournal:
    """One locked journal per prepared output directory, with explicit resume.

    ``identity`` binds the already pinned plan, bindings, source and installed
    package environment. Status is diagnostic: even ``completed`` attempts
    require fresh native artifact/admission verification before reuse.
    """

    def __init__(self, output_root: str | Path, *, identity: Mapping[str, Any], resume: bool = False):
        self.root = Path(output_root).expanduser().resolve()
        self.path = self.root / "execution_journal"
        self.identity = _validate_identity(identity)
        if not isinstance(resume, bool):
            raise ValueError("Journal resume must be an explicit boolean")
        self.resume = resume
        self.session_id = uuid4().hex
        self.identity_sha256: str | None = None
        self._events: list[dict[str, Any]] = []
        self._active = False
        self._clock_starts: dict[str, float] = {}

    def __enter__(self):
        if self._active:
            raise RuntimeError("Journal already holds its executor lock")
        if not self.root.is_dir():
            raise FileNotFoundError("Journal requires the existing prepared output directory")
        self._lock_handle = (self.root / ".execution-journal.lock").open("a+b")
        if self._lock_handle.tell() == 0:
            self._lock_handle.write(b"0")
            self._lock_handle.flush()
        self._lock_handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(self._lock_handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(self._lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            self._lock_handle.close()
            raise RuntimeError("This prepared experiment already has an active executor") from exc
        self._active = True
        self._events, self._clock_starts = [], {}
        try:
            if self.resume:
                self._load()
                # An unclosed start proves only that its old process did not
                # persist a finish. Preserve it and append an interrupted end;
                # the native collector determines whether outputs survived.
                for attempt in self.attempts:
                    if attempt["status"] == "started":
                        self.finish_attempt(attempt["attempt_id"], status="interrupted", details={
                            "reason": "Previous executor did not persist an attempt finish",
                            "elapsed_unavailable_reason": "No elapsed time is inferred across executor processes",
                        })
            else:
                self.path.mkdir(exist_ok=False)
                (self.path / "events").mkdir()
                payload = {"schema_version": 1, "identity": self.identity, "created_utc": _timestamp()}
                self.identity_sha256 = _fingerprint(payload)
                payload["identity_sha256"] = self.identity_sha256
                _publish_exclusive(self.path / "identity.json", payload)
                _publish_exclusive(self.path / "head.json", self._head())
            return self
        except BaseException:
            self.__exit__()
            raise

    def __exit__(self, *args):
        if not self._active:
            return
        try:
            self._lock_handle.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(self._lock_handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(self._lock_handle.fileno(), fcntl.LOCK_UN)
        finally:
            self._lock_handle.close()
            self._active = False

    @property
    def events(self) -> list[dict[str, Any]]:
        return deepcopy(self._events)

    @property
    def attempts(self) -> list[dict[str, Any]]:
        attempts: dict[str, dict[str, Any]] = {}
        for event in self._events:
            attempt_id = event["attempt_id"]
            if event["kind"] == "attempt_started":
                attempts[attempt_id] = {"attempt_id": attempt_id, "stage": event["stage"],
                    "branch_id": event["branch_id"], "status": "started", "start": event, "finish": None}
            else:
                attempts[attempt_id].update(status=event["status"], finish=event)
        return deepcopy(list(attempts.values()))

    def latest_attempt(self, stage: str, branch_id: str | None = None) -> dict[str, Any] | None:
        return next((item for item in reversed(self.attempts)
                     if item["stage"] == stage and item["branch_id"] == branch_id), None)

    def _require_active(self) -> None:
        if not self._active:
            raise RuntimeError("Journal operations require the executor lock")

    def _load(self, *, repair_head: bool = True) -> None:
        payload = json.loads((self.path / "identity.json").read_text(encoding="utf-8"))
        if not isinstance(payload, dict) or set(payload) != {"schema_version", "identity", "created_utc", "identity_sha256"}:
            raise ValueError("Journal identity envelope is invalid")
        digest = payload.pop("identity_sha256")
        if payload["schema_version"] != 1 or _fingerprint(payload) != digest or payload["identity"] != self.identity:
            raise ValueError("Journal frozen plan/source/environment identity changed")
        _validate_timestamp(payload["created_utc"])
        self.identity_sha256 = digest
        files = sorted((self.path / "events").iterdir())
        files = [path for path in files if not (path.name.startswith(".journal-") and path.suffix == ".tmp")]
        for sequence, path in enumerate(files, 1):
            if path.name != f"{sequence:08d}.json" or not path.is_file():
                raise ValueError("Journal event file sequence is incomplete or contains unexpected entries")
            event = json.loads(path.read_text(encoding="utf-8"))
            self._validate_event(event)
            self._events.append(event)
        head = json.loads((self.path / "head.json").read_text(encoding="utf-8"))
        if (not isinstance(head, dict) or set(head) != {"sequence", "event_sha256", "identity_sha256"}
                or type(head["sequence"]) is not int or not 0 <= head["sequence"] <= len(self._events)
                or head["identity_sha256"] != self.identity_sha256
                or head["event_sha256"] != (self._events[head["sequence"] - 1]["event_sha256"]
                    if head["sequence"] else self.identity_sha256)):
            raise ValueError("Journal head differs from its persisted event history")
        if repair_head and head["sequence"] < len(self._events):
            # The last complete event may survive a process exit before the
            # replaceable head was advanced. Never rewrite that event.
            self._save_head()

    def _head(self) -> dict[str, Any]:
        return {"sequence": len(self._events), "identity_sha256": self.identity_sha256,
                "event_sha256": self._events[-1]["event_sha256"] if self._events else self.identity_sha256}

    def _save_head(self) -> None:
        atomic_write_json(self.path / "head.json", self._head(),
                          dump_kwargs={"sort_keys": True, "allow_nan": False})

    def _validate_event(self, event: dict[str, Any]) -> None:
        if not isinstance(event, dict):
            raise ValueError("Journal event must be an object")
        digest = event.get("event_sha256")
        unsigned = {key: value for key, value in event.items() if key != "event_sha256"}
        if (not isinstance(digest, str) or not _DIGEST.fullmatch(digest) or _fingerprint(unsigned) != digest
                or type(event.get("schema_version")) is not int or event["schema_version"] != 1
                or type(event.get("sequence")) is not int or event["sequence"] != len(self._events) + 1
                or event.get("identity_sha256") != self.identity_sha256
                or event.get("previous_event_sha256") != (self._events[-1]["event_sha256"] if self._events else self.identity_sha256)):
            raise ValueError("Journal event identity, sequence or hash chain changed")
        _validate_timestamp(event.get("timestamp_utc"))
        for field in ("session_id", "attempt_id"):
            if not isinstance(event.get(field), str) or not _UUID.fullmatch(event[field]):
                raise ValueError("Journal event invocation identity is invalid")
        if (not isinstance(event.get("stage"), str) or not _IDENTIFIER.fullmatch(event["stage"])
                or event.get("branch_id") is not None and (not isinstance(event["branch_id"], str)
                    or not _IDENTIFIER.fullmatch(event["branch_id"]))
                or not isinstance(event.get("details"), dict)):
            raise ValueError("Journal event stage, branch or details are invalid")
        earlier = next((item for item in self.attempts if item["attempt_id"] == event["attempt_id"]), None)
        if event.get("kind") == "attempt_started":
            if earlier is not None or event.get("status") != "started":
                raise ValueError("Journal attempt start is duplicated or invalid")
            _validate_seconds(event.get("monotonic_started_seconds"))
        elif event.get("kind") == "attempt_finished":
            if (earlier is None or earlier["finish"] is not None or event.get("status") not in _STATUSES
                    or (event["stage"], event["branch_id"]) != (earlier["stage"], earlier["branch_id"])):
                raise ValueError("Journal attempt finish has no matching unique start")
            if event["session_id"] == earlier["start"]["session_id"]:
                _validate_seconds(event.get("monotonic_finished_seconds"))
                _validate_seconds(event.get("elapsed_seconds"))
                expected = event["monotonic_finished_seconds"] - earlier["start"]["monotonic_started_seconds"]
                if expected < 0 or not math.isclose(event["elapsed_seconds"], expected, abs_tol=1e-12, rel_tol=1e-12):
                    raise ValueError("Journal elapsed time differs from its same-process monotonic clocks")
            elif (event["status"] != "interrupted" or event.get("elapsed_seconds") is not None
                    or event.get("monotonic_finished_seconds") is not None):
                raise ValueError("Journal cannot infer elapsed time or success across executor processes")
        else:
            raise ValueError("Unknown journal event kind")

    def _append(self, event: dict[str, Any]) -> None:
        self._require_active()
        event.update(schema_version=1, sequence=len(self._events) + 1,
                     identity_sha256=self.identity_sha256, session_id=self.session_id,
                     timestamp_utc=_timestamp(), previous_event_sha256=(
                         self._events[-1]["event_sha256"] if self._events else self.identity_sha256))
        event["event_sha256"] = _fingerprint(event)
        self._validate_event(event)
        _publish_exclusive(self.path / "events" / f"{event['sequence']:08d}.json", event)
        self._events.append(event)
        self._save_head()

    def start_attempt(self, stage: str, branch_id: str | None = None, *, details: Mapping[str, Any] | None = None) -> str:
        self._require_active()
        if any(item["stage"] == stage and item["branch_id"] == branch_id and item["status"] == "started"
               for item in self.attempts):
            raise ValueError("This journal stage and branch already have an unfinished attempt")
        attempt_id, started = uuid4().hex, time.monotonic()
        self._append({"kind": "attempt_started", "attempt_id": attempt_id, "stage": stage,
                      "branch_id": branch_id, "status": "started", "monotonic_started_seconds": started,
                      "details": deepcopy(dict(details or {}))})
        self._clock_starts[attempt_id] = started
        return attempt_id

    def finish_attempt(self, attempt_id: str, *, status: str, details: Mapping[str, Any] | None = None) -> None:
        self._require_active()
        attempt = next((item for item in self.attempts if item["attempt_id"] == attempt_id), None)
        if attempt is None or attempt["finish"] is not None:
            raise ValueError("Attempt finish requires a known unfinished start")
        finished = time.monotonic() if attempt_id in self._clock_starts else None
        elapsed = finished - self._clock_starts[attempt_id] if finished is not None else None
        self._append({"kind": "attempt_finished", "attempt_id": attempt_id, "stage": attempt["stage"],
                      "branch_id": attempt["branch_id"], "status": status,
                      "monotonic_finished_seconds": finished, "elapsed_seconds": elapsed,
                      "details": deepcopy(dict(details or {}))})
        self._clock_starts.pop(attempt_id, None)


__all__ = ["ExecutionJournal"]
