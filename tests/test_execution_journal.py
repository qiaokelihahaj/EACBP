"""Crash/restart diagnostics preserve native completion without admitting it."""

from copy import deepcopy
import json
import subprocess
import sys

import pytest

from eacbp.application import execution_journal as module
from eacbp.application.execution_journal import ExecutionJournal


@pytest.fixture
def identity():
    return {"bindings_sha256": "a" * 64, "plan_file_sha256": "b" * 64,
            "plan_sha256": "c" * 64, "source_sha256": "d" * 64,
            "execution_environment": {"python": "frozen", "source_sha256": "e" * 64}}


def _bytes(root):
    return {str(path.relative_to(root)): path.read_bytes()
            for path in (root / "execution_journal").rglob("*.json") if path.name != "head.json"}


def _rewrite_event(path, mutate):
    event = json.loads(path.read_text(encoding="utf-8"))
    mutate(event)
    event.pop("event_sha256")
    event["event_sha256"] = module._fingerprint(event)
    path.write_text(json.dumps(event), encoding="utf-8")


def test_postprocessing_failure_retains_native_completion_and_precise_elapsed(tmp_path, identity, monkeypatch):
    clocks = iter([100., 104.5, 105., 106.])
    monkeypatch.setattr(module.time, "monotonic", lambda: next(clocks))
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        native = journal.start_attempt("native", "primary", details={"source": "frozen"})
        journal.finish_attempt(native, status="completed", details={"admitted": True})
        report = journal.start_attempt("report")
        journal.finish_attempt(report, status="postprocessing_failed", details={"error": "AttributeError", "log": "failed-v1.log"})
        assert journal.latest_attempt("native", "primary")["status"] == "completed"
        assert journal.latest_attempt("report")["status"] == "postprocessing_failed"
        assert journal.events[1]["elapsed_seconds"] == 4.5
        assert journal.events[3]["elapsed_seconds"] == 1.
        assert all(event["timestamp_utc"].endswith("+00:00") for event in journal.events)
        assert all(event["identity_sha256"] == journal.identity_sha256 for event in journal.events)
    before = _bytes(tmp_path)
    monkeypatch.undo()
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["finish"]["details"]["admitted"] is True
        attempt = journal.start_attempt("report", details={"recovery": True})
        journal.finish_attempt(attempt, status="completed")
        assert journal.latest_attempt("report")["status"] == "completed"
    assert all(_bytes(tmp_path)[name] == data for name, data in before.items())


def test_pending_attempt_becomes_interrupted_with_no_invented_elapsed(tmp_path, identity):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
    before = _bytes(tmp_path)
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        recovered = journal.latest_attempt("native", "primary")
        assert recovered["attempt_id"] == attempt
        assert recovered["status"] == "interrupted"
        assert recovered["finish"]["elapsed_seconds"] is None
        assert recovered["finish"]["monotonic_finished_seconds"] is None
        assert recovered["start"]["session_id"] != recovered["finish"]["session_id"]
        assert "across executor processes" in recovered["finish"]["details"]["elapsed_unavailable_reason"]
    assert all(_bytes(tmp_path)[name] == data for name, data in before.items())
    final = _bytes(tmp_path)
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert len(journal.events) == 2
    assert _bytes(tmp_path) == final


def test_real_process_exit_releases_lock_and_preserves_started_event(tmp_path, identity):
    program = """
import json, os, sys
from eacbp.application.execution_journal import ExecutionJournal
with ExecutionJournal(sys.argv[1], identity=json.loads(sys.argv[2])) as journal:
    journal.start_attempt('native', 'primary', details={'fault': 'process_exit'})
    os._exit(23)
"""
    result = subprocess.run([sys.executable, "-B", "-c", program, str(tmp_path), json.dumps(identity)],
                            capture_output=True, text=True)
    assert result.returncode == 23, result.stderr
    original = _bytes(tmp_path)
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["status"] == "interrupted"
        assert journal.latest_attempt("native", "primary")["finish"]["elapsed_seconds"] is None
    assert all(_bytes(tmp_path)[name] == data for name, data in original.items())


def test_running_executor_excludes_second_process(tmp_path, identity):
    program = """
import json, sys
from eacbp.application.execution_journal import ExecutionJournal
try:
    with ExecutionJournal(sys.argv[1], identity=json.loads(sys.argv[2]), resume=True):
        pass
except RuntimeError as error:
    print(error)
    sys.exit(19)
"""
    with ExecutionJournal(tmp_path, identity=identity):
        result = subprocess.run([sys.executable, "-B", "-c", program, str(tmp_path), json.dumps(identity)],
                                capture_output=True, text=True)
        assert result.returncode == 19, result.stderr
        assert "active executor" in result.stdout


@pytest.mark.parametrize("field", ["bindings_sha256", "plan_file_sha256", "plan_sha256", "source_sha256", "execution_environment"])
def test_resume_rejects_every_changed_frozen_identity(tmp_path, identity, field):
    with ExecutionJournal(tmp_path, identity=identity):
        pass
    changed = deepcopy(identity)
    changed[field] = {"source_sha256": "f" * 64} if field == "execution_environment" else "f" * 64
    original = _bytes(tmp_path)
    with pytest.raises(ValueError, match="identity changed"):
        with ExecutionJournal(tmp_path, identity=changed, resume=True):
            pass
    assert _bytes(tmp_path) == original


def test_completed_status_is_only_diagnostic_and_never_invokes_native_code(tmp_path, identity):
    # A claimed admission in an event cannot become a scientific artifact. The
    # package integration must use collect_clustering_run again before reuse.
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
        journal.finish_attempt(attempt, status="completed", details={"admitted": True})
    assert not (tmp_path / "runs").exists()
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["status"] == "completed"
    assert not (tmp_path / "runs").exists()


def test_partial_publication_failure_cannot_replace_old_event(tmp_path, identity, monkeypatch):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
        original = _bytes(tmp_path)
        def fail_link(*args):
            raise OSError("injected publication failure")
        monkeypatch.setattr(module.os, "link", fail_link)
        with pytest.raises(OSError, match="publication failure"):
            journal.finish_attempt(attempt, status="completed")
        assert journal.latest_attempt("native", "primary")["status"] == "started"
        assert _bytes(tmp_path) == original
        assert not list((tmp_path / "execution_journal/events").glob("*.tmp"))
    monkeypatch.undo()
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["status"] == "interrupted"


def test_unpublished_crash_temporary_is_preserved_and_ignored(tmp_path, identity):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
        journal.finish_attempt(attempt, status="completed")
    leftover = tmp_path / "execution_journal/events/.journal-crashed.tmp"
    leftover.write_bytes(b'{"incomplete":')
    original = _bytes(tmp_path)
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["status"] == "completed"
    assert leftover.read_bytes() == b'{"incomplete":'
    assert _bytes(tmp_path) == original


@pytest.mark.parametrize("attack", ["bytes", "chain", "sequence_gap", "malformed_complete", "elapsed", "session_success", "trailing_deletion"])
def test_tampered_published_events_fail_closed(tmp_path, identity, attack):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
        journal.finish_attempt(attempt, status="completed")
    path = tmp_path / "execution_journal/events/00000002.json"
    if attack == "bytes":
        path.write_bytes(path.read_bytes().replace(b"completed", b"failed"))
    elif attack == "chain":
        _rewrite_event(path, lambda event: event.update(previous_event_sha256="f" * 64))
    elif attack == "sequence_gap":
        path.rename(path.with_name("00000003.json"))
    elif attack == "malformed_complete":
        path.write_bytes(b'{"broken":\n')
    elif attack == "elapsed":
        _rewrite_event(path, lambda event: event.update(elapsed_seconds=123.))
    elif attack == "trailing_deletion":
        path.unlink()
    else:
        _rewrite_event(path, lambda event: event.update(session_id="f" * 32,
                                                       elapsed_seconds=None, monotonic_finished_seconds=None))
    original = _bytes(tmp_path)
    with pytest.raises(ValueError):
        with ExecutionJournal(tmp_path, identity=identity, resume=True):
            pass
    assert _bytes(tmp_path) == original


def test_process_exit_between_event_publication_and_head_update_rolls_forward(tmp_path, identity):
    program = """
import json, os, sys
from eacbp.application.execution_journal import ExecutionJournal
with ExecutionJournal(sys.argv[1], identity=json.loads(sys.argv[2])) as journal:
    journal._save_head = lambda: os._exit(24)
    journal.start_attempt('native', 'primary')
"""
    result = subprocess.run([sys.executable, "-B", "-c", program, str(tmp_path), json.dumps(identity)],
                            capture_output=True, text=True)
    assert result.returncode == 24, result.stderr
    head_path = tmp_path / "execution_journal/head.json"
    assert json.loads(head_path.read_text())["sequence"] == 0
    original = _bytes(tmp_path)
    with ExecutionJournal(tmp_path, identity=identity, resume=True) as journal:
        assert journal.latest_attempt("native", "primary")["status"] == "interrupted"
    assert json.loads(head_path.read_text())["sequence"] == 2
    assert all(_bytes(tmp_path)[name] == data for name, data in original.items())


def test_exclusive_initialization_and_no_implicit_resume(tmp_path, identity):
    with ExecutionJournal(tmp_path, identity=identity):
        pass
    original = _bytes(tmp_path)
    with pytest.raises(FileExistsError):
        with ExecutionJournal(tmp_path, identity=identity):
            pass
    assert _bytes(tmp_path) == original
    with pytest.raises(FileNotFoundError):
        with ExecutionJournal(tmp_path / "missing", identity=identity, resume=True):
            pass


def test_invalid_or_duplicate_attempt_cannot_publish_event(tmp_path, identity):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary")
        original = _bytes(tmp_path)
        with pytest.raises(ValueError, match="unfinished attempt"):
            journal.start_attempt("native", "primary")
        with pytest.raises(ValueError):
            journal.start_attempt("../escape")
        with pytest.raises(ValueError):
            journal.finish_attempt(attempt, status="fabricated")
        with pytest.raises(ValueError):
            journal.finish_attempt("f" * 32, status="completed")
        with pytest.raises(ValueError):
            journal.start_attempt("report", details={"invalid": float("nan")})
        assert _bytes(tmp_path) == original
        journal.finish_attempt(attempt, status="failed", details={"reason": "original native failure"})
        final = _bytes(tmp_path)
        with pytest.raises(ValueError):
            journal.finish_attempt(attempt, status="completed")
        assert _bytes(tmp_path) == final
    with pytest.raises(RuntimeError, match="executor lock"):
        journal.start_attempt("report")


def test_returned_events_and_attempts_cannot_mutate_journal(tmp_path, identity):
    with ExecutionJournal(tmp_path, identity=identity) as journal:
        attempt = journal.start_attempt("native", "primary", details={"reason": "preserved"})
        journal.events[0]["details"]["reason"] = "tampered"
        journal.attempts[0]["start"]["details"]["reason"] = "tampered"
        assert journal.events[0]["details"]["reason"] == "preserved"
        journal.finish_attempt(attempt, status="completed")
