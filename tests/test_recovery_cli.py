"""Recovery commands retain explicit pins and actionable failure results."""
import json
from pathlib import Path

import pytest

from eacbp.cli import main


@pytest.mark.parametrize("command", ["recovery-status", "resume-experiment"])
def test_recovery_cli_requires_independently_supplied_pin(tmp_path, capsys, command):
    with pytest.raises(SystemExit) as exc:
        main([command, "--prepared", str(tmp_path / "missing")])
    assert exc.value.code == 2
    assert "--bindings-sha256" in capsys.readouterr().err
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("command,status,can_resume,exit_code", [
    ("recovery-status", "ready", True, 0),
    ("recovery-status", "ready", False, 1),
    ("recovery-status", "blocked", False, 1),
    ("resume-experiment", "blocked", False, 1),
    ("resume-experiment", "partial", False, 1),
    ("resume-experiment", "failed", False, 1),
    ("resume-experiment", "completed", False, 0),
])
def test_recovery_cli_preserves_failures_and_completion_scope(
    tmp_path, monkeypatch, capsys, command, status, can_resume, exit_code
):
    import eacbp.application.recovery as recovery
    payload = {"status": status, "can_resume": can_resume,
               "reuse_authorized": False, "external_scientific_validation": False,
               "blockers": [{"reason": "source bytes changed"}] if status == "blocked" else [],
               "failed_attempts": [{"stage": "postprocessing", "reason": "presentation unavailable"}]}
    calls = []

    def invoked(prepared, *, expected_bindings_sha256):
        calls.append((prepared, expected_bindings_sha256))
        return payload

    def forbidden(*args, **kwargs):
        pytest.fail("The other operation must not run")

    selected = "inspect_recovery" if command == "recovery-status" else "resume_recovery"
    other = "resume_recovery" if command == "recovery-status" else "inspect_recovery"
    monkeypatch.setattr(recovery, selected, invoked)
    monkeypatch.setattr(recovery, other, forbidden)
    prepared = tmp_path / "existing-experiment"
    pin = "a" * 64
    assert main([command, "--prepared", str(prepared), "--bindings-sha256", pin]) == exit_code
    captured = capsys.readouterr()
    assert json.loads(captured.out) == payload
    assert captured.err == ""
    assert calls == [(Path(prepared), pin)]
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("command", ["recovery-status", "resume-experiment"])
def test_recovery_cli_reports_invalid_pin_without_creating_target(tmp_path, capsys, command):
    target = tmp_path / "never-created"
    assert main([command, "--prepared", str(target), "--bindings-sha256", "not-a-sha256"]) == 1
    captured = capsys.readouterr()
    assert not target.exists()
    # The shared gate may return a structured refusal or raise a concise error.
    if captured.out:
        result = json.loads(captured.out)
        assert result["status"] == "blocked" and not result["can_resume"]
    else:
        assert "sha256" in captured.err.casefold() or "sha-256" in captured.err.casefold()
