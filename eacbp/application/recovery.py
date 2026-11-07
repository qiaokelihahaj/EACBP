"""Read-only recovery inspection and explicit frozen-experiment continuation.

Inspection verifies persisted diagnostics and pinned evidence. It never opens
an artifact registry or grants scientific admission; the native collector must
re-admit every existing branch when the executor actually resumes.
"""

from __future__ import annotations

from contextlib import contextmanager
import json
import os
from pathlib import Path
import re
from typing import Any, Mapping

from eacbp.application import clustering_sensitivity_execution as execution
from eacbp.application.execution_journal import ExecutionJournal


_SHA256 = re.compile(r"[0-9a-f]{64}")


def _require_digest(value: Any) -> None:
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise ValueError("Recovery requires an explicit 64-character lowercase bindings SHA256")


@contextmanager
def _inspection_lock(root: Path):
    # Opening an existing file without append/create keeps inspection read-only.
    # Match the executor's byte-lock access without append/create. No file data
    # is written, and a missing executor lock is never created.
    with (root / ".execution-journal.lock").open("r+b") as stream:
        stream.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError("This prepared experiment already has an active executor") from exc
        try:
            if stream.read(1) != b"0":
                raise ValueError("Persisted executor lock is empty or invalid")
            yield
        finally:
            stream.seek(0)
            if os.name == "nt":
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def _pinned_file(path: Any, digest: Any, run_dir: Path, *, artifact: bool = False) -> None:
    from eacbp.application import study_service

    # ArtifactMetadata retains the storage checksum's algorithm tag. Frozen
    # bindings and native file inventories use bare SHA256 and stay strict.
    if artifact and isinstance(digest, str) and digest.startswith("sha256:"):
        digest = digest[len("sha256:"):]
    if not isinstance(path, str) or not isinstance(digest, str) or not _SHA256.fullmatch(digest):
        raise ValueError("Previous native admission has an invalid file SHA pin")
    resolved = Path(path).expanduser().resolve()
    resolved.relative_to(run_dir)
    if not resolved.is_file() or study_service._hash_file(resolved) != digest:
        raise ValueError("Previously admitted native evidence is missing or changed")


def _branch_state(binding: Mapping[str, Any], journal: ExecutionJournal) -> dict[str, Any]:
    state = {"branch_id": binding["branch_id"], "planned_status": binding["planned_status"],
             "reuse_authorized": False, "requires_native_readmission": True,
             "prior_pin_available": False, "native_metadata_pins_verified": False,
             "native_payload_pins_verified": 0}
    if binding["planned_status"] != "planned":
        return {**state, "state": "not_planned", "reason": binding.get("reason")}
    run_dir = Path(binding["run_dir"]).resolve()
    exists = run_dir.exists()
    if exists and not run_dir.is_dir():
        raise ValueError("Native run path is not a directory")
    state.update(run_dir=str(run_dir), native_exists=exists,
                 state="pending_re_admission" if exists else "new_execution")
    admitted = [attempt for attempt in journal.attempts
                if attempt["branch_id"] == binding["branch_id"]
                and attempt["stage"] in {"native_execution", "native_admission"}
                and attempt["status"] == "completed"
                and attempt["finish"]["details"].get("admitted") is True]
    if not admitted:
        return state
    if not exists:
        raise ValueError("Previously admitted native directory is missing; recovery cannot recompute it")
    previous = admitted[-1]
    provenance = previous["finish"]["details"]["provenance"]
    if (not isinstance(provenance, dict) or not isinstance(provenance.get("run_dir"), str)
            or Path(provenance["run_dir"]).resolve() != run_dir):
        raise ValueError("Previous native admission belongs to another directory")
    files = provenance.get("files")
    if not isinstance(files, dict) or not files:
        raise ValueError("Previous native admission file pins are invalid")
    for path, digest in files.items():
        _pinned_file(path, digest, run_dir)
    if execution._native_file_hashes(run_dir) != files:
        raise ValueError("Previously admitted native file inventory changed")
    artifacts = provenance.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ValueError("Previous native admission payload pins are invalid")
    for artifact in artifacts:
        if not isinstance(artifact, dict):
            raise ValueError("Previous native admission payload pin is invalid")
        _pinned_file(artifact["path"], artifact["sha256"], run_dir, artifact=True)
    state.update(prior_pin_available=True, prior_native_attempt_id=previous["attempt_id"],
                 native_metadata_pins_verified=True, native_payload_pins_verified=len(artifacts))
    return state


def inspect_recovery(prepared, *, expected_bindings_sha256: str) -> dict[str, Any]:
    """Inspect an already executed preparation without changing any files.

    ``ready`` means the frozen inputs and persisted recovery evidence permit a
    continuation attempt. Existing native directories remain pending admission,
    including those whose old journal attempt says ``completed``. A surviving
    event ahead of the head and an unfinished start are reported without repair.
    """
    report: dict[str, Any] = {
        "schema_version": 1, "status": "blocked", "can_resume": False,
        "output_dir": None, "expected_bindings_sha256": expected_bindings_sha256,
        "blockers": [], "identity": {"verified": False}, "journal": None,
        "branches": [], "latest_comparison": None, "unfinished_attempts": [], "failed_attempts": [],
        "reuse_authorized": False, "requires_native_readmission": True,
        "integrity_scope": "Frozen preparation and journal identity/hash chain; prior admitted native file and payload pins. Detached historical comparisons/reports are diagnostic, not native admission.",
        "biological_stability_claim_supported": False,
    }
    check = "bindings_sha256"
    try:
        _require_digest(expected_bindings_sha256)
        check = "preparation_path"
        path = Path(prepared["bindings_path"] if isinstance(prepared, Mapping) else prepared).expanduser().resolve()
        candidate = path if path.is_dir() else path.parent
        report["output_dir"] = str(candidate)
        if not candidate.is_dir():
            raise FileNotFoundError("Recovery requires an existing prepared output directory")
        check = "executor_lock"
        with _inspection_lock(candidate):
            check = "preparation_integrity"
            root, bindings, plan, bindings_sha256 = execution._load_preparation(prepared, expected_bindings_sha256)
            if root != candidate:
                raise ValueError("Recovery lock and prepared output directory disagree")
            check = "frozen_input_identity"
            execution._validate_frozen_inputs(root, bindings, plan)
            identity = {key: bindings[key] for key in ("plan_file_sha256", "plan_sha256", "source_sha256", "execution_environment")}
            identity["bindings_sha256"] = bindings_sha256
            check = "journal_integrity"
            journal = ExecutionJournal(root, identity=identity, resume=True)
            journal._load(repair_head=False)
            head = json.loads((journal.path / "head.json").read_text(encoding="utf-8"))
            attempts = journal.attempts
            report.update(identity={"verified": True, **identity},
                          journal={"path": str(journal.path), "identity_sha256": journal.identity_sha256,
                                   "event_count": len(journal.events), "head_sequence": head["sequence"],
                                   "head_needs_repair": head["sequence"] < len(journal.events), "attempts": attempts},
                          unfinished_attempts=[item for item in attempts if item["status"] == "started"],
                          failed_attempts=[item for item in attempts if item["status"] in {"failed", "interrupted", "postprocessing_failed"}])
            comparison = next((item for item in reversed(attempts)
                               if item["stage"] == "comparison" and item["status"] == "completed"), None)
            if comparison is not None:
                details = comparison["finish"]["details"]
                report["latest_comparison"] = {"attempt_id": comparison["attempt_id"], "path": details.get("path"),
                                               "sha256": details.get("sha256"), "status": details.get("status"),
                                               "integrity_checked": False, "reuse_authorized": False}
            check = "native_evidence_integrity"
            for binding in bindings["branches"]:
                report["branches"].append(_branch_state(binding, journal))
            report.update(status="ready", can_resume=True)
    except Exception as exc:
        report["blockers"].append({"code": check, "reason": str(exc), "error_type": type(exc).__name__})
    return report


def resume_recovery(prepared, *, expected_bindings_sha256: str) -> dict[str, Any]:
    """Continue only after inspection; the executor locks and verifies again."""
    report = inspect_recovery(prepared, expected_bindings_sha256=expected_bindings_sha256)
    if not report["can_resume"]:
        return report
    # Do not hold the inspection lock across delegation. The existing executor
    # acquires the same lock and repeats identity checks, so a competing executor
    # between inspection and continuation is rejected rather than authorized.
    return execution.execute_clustering_sensitivity(prepared, expected_bindings_sha256=expected_bindings_sha256, resume=True)


__all__ = ["inspect_recovery", "resume_recovery"]
