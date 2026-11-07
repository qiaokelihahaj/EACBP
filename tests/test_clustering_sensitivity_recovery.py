"""Explicit recovery re-admits native outputs and preserves prior failures."""

from copy import deepcopy
import json
from pathlib import Path

import pytest

from eacbp.application import clustering_sensitivity_execution as execution
from eacbp.application import study_service

from test_clustering_sensitivity_execution import _fake_receipt, _read, _source


BRANCHES = ["primary", "resolution_0_5", "resolution_2_0", "clustering_seed_17", "uniform_sample_80pct"]
ENVIRONMENT = {"source_sha256": "e" * 64, "python": "frozen-test-python"}


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "execution_environment", lambda: deepcopy(ENVIRONMENT))
    source, config, manifest = _source(tmp_path)
    record = execution.prepare_clustering_sensitivity(manifest=manifest, config=config, data=source,
                                                     output_dir=tmp_path / "experiment")
    return record, source


def _stub_native(monkeypatch, *, interrupted=None, incomplete=None, collect_error=None):
    started, collected = [], []

    def run(**kwargs):
        branch = kwargs["config"]["autonomous_analysis"]["clustering_sensitivity_branch"]["branch_id"]
        started.append(branch)
        path = Path(kwargs["run_dir"])
        path.mkdir(parents=True, exist_ok=False)
        (path / "summary.json").write_text(json.dumps({"branch_id": branch, "status": "success"}), encoding="utf-8")
        (path / "report.md").write_text(f"Original native report: {branch}\n", encoding="utf-8")
        if branch == interrupted:
            raise SystemExit("Injected exit after native files were persisted")
        return {"status": "success"}

    def collect(path, *, plan, binding):
        branch = binding["branch_id"]
        collected.append(branch)
        if branch == collect_error:
            raise ValueError("Current native admission rejected")
        receipt = _fake_receipt(plan, binding, path)
        if branch == incomplete:
            receipt.update(admitted=False, reason="Current native request is incomplete")
        return receipt

    monkeypatch.setattr(study_service, "run_study", run)
    monkeypatch.setattr(execution, "collect_clustering_run", collect)
    monkeypatch.setattr("eacbp.evidence.snapshot.load_study_snapshot", lambda *args, **kwargs: None)
    monkeypatch.setattr(study_service, "_make_registry", lambda path: None)
    return started, collected


def _attempts(root):
    attempts = {}
    for path in sorted((root / "execution_journal/events").glob("*.json")):
        event = _read(path)
        if event["kind"] == "attempt_started":
            attempts[event["attempt_id"]] = {"stage": event["stage"], "branch_id": event["branch_id"],
                                             "start": event, "finish": None, "status": "started"}
        else:
            attempts[event["attempt_id"]].update(finish=event, status=event["status"])
    return list(attempts.values())


def _immutable_bytes(root):
    paths = []
    for directory in ("runs", "comparisons", "execution_journal/events", "execution_journal/resources"):
        paths.extend(path for path in (root / directory).rglob("*") if path.is_file())
    paths.extend(path for name in ("comparison.json", "report.md", "execution_journal/identity.json")
                 if (path := root / name).is_file())
    return {path: path.read_bytes() for path in paths}


def _assert_preserved(files):
    assert all(path.is_file() and path.read_bytes() == content for path, content in files.items())


def test_failed_postprocessing_after_all_native_runs_recovers_without_reexecution(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch)
    render = execution._markdown
    monkeypatch.setattr(execution, "_markdown", lambda result: (_ for _ in ()).throw(RuntimeError("Injected report failure")))
    with pytest.raises(RuntimeError, match="report failure"):
        execution.execute_clustering_sensitivity(record)
    assert started == BRANCHES and collected == BRANCHES
    original = _immutable_bytes(root)
    old_comparison = (root / "comparison.json").read_bytes()
    prior = _attempts(root)
    assert len([item for item in prior if item["stage"] == "native_execution" and item["status"] == "completed"]) == 5
    assert any(item["stage"] == "postprocessing" and item["status"] == "postprocessing_failed" for item in prior)
    monkeypatch.setattr(execution, "_markdown", render)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] == "completed"
    assert recovered["native_reused_branch_ids"] == BRANCHES
    assert recovered["native_started_branch_ids"] == []
    assert started == BRANCHES and collected == BRANCHES * 2
    assert (root / "comparison.json").read_bytes() == old_comparison
    original_science, recovered_science = json.loads(old_comparison), _read(recovered["comparison"])
    for field in ("plan_sha256", "branches", "donor_perturbation", "native_execution_receipts", "run_provenance"):
        assert recovered_science[field] == original_science[field]
    assert recovered_science["resume_prior_native_file_pins"] == {branch: True for branch in BRANCHES}
    assert _read(root / "latest_comparison.json")["invocation_id"] == recovered["invocation_id"]
    _assert_preserved(original)
    attempts = _attempts(root)
    assert len([item for item in attempts if item["stage"] == "native_execution"]) == 5
    admitted = [item for item in attempts if item["stage"] == "native_admission"]
    assert len(admitted) == 5 and all(item["status"] == "completed" for item in admitted)
    assert all(item["finish"]["elapsed_seconds"] >= 0 for item in admitted)
    assert len([item for item in attempts if item["stage"] == "postprocessing"]) == 2


def test_exit_after_first_native_resumes_by_admission_before_remaining_computation(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch, interrupted="primary")
    with pytest.raises(SystemExit, match="persisted"):
        execution.execute_clustering_sensitivity(record)
    assert started == ["primary"] and collected == []
    original = _immutable_bytes(root)
    initial = [item for item in _attempts(root) if item["stage"] == "native_execution"]
    assert len(initial) == 1 and initial[0]["status"] == "started"
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] == "completed"
    assert started == BRANCHES and collected == BRANCHES
    assert recovered["native_reused_branch_ids"] == ["primary"]
    assert recovered["native_started_branch_ids"] == BRANCHES[1:]
    assert _read(recovered["comparison"])["resumed_without_previous_native_file_pin"] == ["primary"]
    _assert_preserved(original)
    interrupted = [item for item in _attempts(root) if item["stage"] == "native_execution" and item["status"] == "interrupted"]
    assert len(interrupted) == 1
    assert interrupted[0]["finish"]["elapsed_seconds"] is None
    assert interrupted[0]["finish"]["monotonic_finished_seconds"] is None


def test_failed_comparison_is_preserved_when_all_native_outputs_are_readmitted(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch)
    compare = execution.compare_clustering_sensitivity_runs
    monkeypatch.setattr(execution, "compare_clustering_sensitivity_runs",
                        lambda *args, **kwargs: (_ for _ in ()).throw(ValueError("Injected comparison failure")))
    failed = execution.execute_clustering_sensitivity(record)
    assert failed["status"] == "failed" and started == BRANCHES and collected == BRANCHES
    original = _immutable_bytes(root)
    assert "comparison failure" in _read(root / "comparison.json")["reason"]
    assert any(item["stage"] == "comparison" and item["status"] == "failed" for item in _attempts(root))
    monkeypatch.setattr(execution, "compare_clustering_sensitivity_runs", compare)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] == "completed" and started == BRANCHES and collected == BRANCHES * 2
    assert recovered["native_reused_branch_ids"] == BRANCHES and recovered["native_started_branch_ids"] == []
    assert _read(root / "comparison.json")["status"] == "failed"
    assert _read(root / "latest_comparison.json")["status"] == "completed"
    _assert_preserved(original)


def test_known_comparison_serialization_failure_has_an_end_and_preserves_native_completion(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch)
    write = execution._write

    def fail_comparison(path, payload):
        if Path(path).name == "comparison.json":
            raise TypeError("Injected comparison JSON format failure")
        return write(path, payload)

    monkeypatch.setattr(execution, "_write", fail_comparison)
    with pytest.raises(TypeError, match="JSON format failure"):
        execution.execute_clustering_sensitivity(record)
    assert started == BRANCHES and collected == BRANCHES
    original = _immutable_bytes(root)
    assert len([item for item in _attempts(root) if item["stage"] == "native_execution" and item["status"] == "completed"]) == 5
    failed = [item for item in _attempts(root) if item["stage"] == "comparison"]
    assert len(failed) == 1 and failed[0]["status"] == "failed"
    assert failed[0]["finish"]["details"]["error_type"] == "TypeError"
    monkeypatch.setattr(execution, "_write", write)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] == "completed" and started == BRANCHES and collected == BRANCHES * 2
    assert recovered["native_reused_branch_ids"] == BRANCHES and recovered["native_started_branch_ids"] == []
    _assert_preserved(original)


@pytest.mark.parametrize("mode", ["incomplete", "collector_error"])
def test_existing_unadmitted_primary_is_preserved_without_reuse_or_alternate_rescue(prepared, monkeypatch, mode):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch, incomplete="primary" if mode == "incomplete" else None,
                                      collect_error="primary" if mode == "collector_error" else None)
    first = execution.execute_clustering_sensitivity(record)
    assert first["status"] == "partial" and started == ["primary"]
    original = _immutable_bytes(root)
    second = execution.execute_clustering_sensitivity(record, resume=True)
    assert second["status"] == "partial" and started == ["primary"]
    assert collected == ["primary", "primary"]
    assert second["native_reused_branch_ids"] == second["native_started_branch_ids"] == []
    comparison = _read(second["comparison"])
    assert comparison["primary_native_summary"]["admitted"] is False
    assert all(item["status"] == "not_available" for item in comparison["branches"])
    assert len([item for item in comparison["failures"] if item["status"] == "not_run"]) == 4
    _assert_preserved(original)


@pytest.mark.parametrize("mode", ["incomplete", "collector_error"])
def test_old_completed_journal_never_overrides_current_failed_admission(prepared, monkeypatch, mode):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, _ = _stub_native(monkeypatch)
    assert execution.execute_clustering_sensitivity(record)["status"] == "completed"
    original = _immutable_bytes(root)
    current_collect = execution.collect_clustering_run

    def reject_primary(path, *, plan, binding):
        if binding["branch_id"] == "primary":
            if mode == "collector_error":
                raise ValueError("Admission changed after original completion")
            receipt = current_collect(path, plan=plan, binding=binding)
            return {**receipt, "admitted": False, "reason": "Current admission is incomplete"}
        pytest.fail("An inadmissible primary must block alternate branches")

    monkeypatch.setattr(execution, "collect_clustering_run", reject_primary)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] == "partial" and started == BRANCHES
    assert recovered["native_reused_branch_ids"] == recovered["native_started_branch_ids"] == []
    assert _read(recovered["comparison"])["primary_native_summary"]["admitted"] is False
    assert len([item for item in _attempts(root) if item["stage"] == "native_execution" and item["status"] == "completed"]) == 5
    assert any(item["stage"] == "native_admission" and item["status"] == "failed" for item in _attempts(root))
    _assert_preserved(original)


def test_native_byte_change_is_blocked_even_when_new_collector_would_admit(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, _ = _stub_native(monkeypatch)
    assert execution.execute_clustering_sensitivity(record)["status"] == "completed"
    report = root / "runs/primary/report.md"
    report.write_bytes(report.read_bytes() + b"changed later\n")
    original = _immutable_bytes(root)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] != "completed"
    assert started == BRANCHES and recovered["native_started_branch_ids"] == []
    assert recovered["native_reused_branch_ids"] == []
    _assert_preserved(original)


def test_missing_previously_completed_native_directory_cannot_trigger_recomputation(prepared, monkeypatch, tmp_path):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, _ = _stub_native(monkeypatch)
    assert execution.execute_clustering_sensitivity(record)["status"] == "completed"
    primary, preserved = root / "runs/primary", root / "preserved-missing-primary"
    primary.resolve().relative_to(tmp_path.resolve())
    preserved.resolve().relative_to(tmp_path.resolve())
    primary.rename(preserved)
    recovered = execution.execute_clustering_sensitivity(record, resume=True)
    assert recovered["status"] != "completed"
    assert started == BRANCHES
    assert not primary.exists()
    assert (preserved / "summary.json").is_file()


@pytest.mark.parametrize("target", ["original_source", "derived_source", "config", "manifest", "environment", "bindings", "plan"])
def test_frozen_identity_mutation_on_resume_precedes_any_new_native_execution(prepared, monkeypatch, target):
    record, source = prepared
    root = Path(record["output_dir"])
    started, _ = _stub_native(monkeypatch)
    assert execution.execute_clustering_sensitivity(record)["status"] == "completed"
    original = _immutable_bytes(root)
    bindings = _read(record["bindings_path"])
    if target == "environment":
        monkeypatch.setattr(execution, "execution_environment", lambda: {**ENVIRONMENT, "source_sha256": "f" * 64})
    else:
        path = {"original_source": source, "derived_source": Path(bindings["branches"][-1]["data_path"]),
                "config": Path(bindings["branches"][1]["config_path"]),
                "manifest": Path(bindings["branches"][1]["manifest_path"]),
                "bindings": Path(record["bindings_path"]), "plan": Path(record["plan_path"])}[target]
        # Even JSON whitespace is a changed frozen file, without relying on a
        # parse error to reject semantically equivalent modified configuration.
        path.write_bytes(path.read_bytes() + b" ")
    if target in {"bindings", "plan"}:
        with pytest.raises(ValueError):
            execution.execute_clustering_sensitivity(record, resume=True)
    else:
        recovered = execution.execute_clustering_sensitivity(record, resume=True)
        assert recovered["status"] == "failed"
        assert recovered["native_started_branch_ids"] == recovered["native_reused_branch_ids"] == []
    assert started == BRANCHES
    _assert_preserved(original)


@pytest.mark.parametrize("resume", [0, 1, "true"])
def test_resume_accepts_only_explicit_boolean_before_execution(prepared, monkeypatch, resume):
    record, _ = prepared
    started, _ = _stub_native(monkeypatch)
    with pytest.raises(ValueError, match="explicit boolean"):
        execution.execute_clustering_sensitivity(record, resume=resume)
    assert started == []


def test_fresh_directory_cannot_be_adopted_as_a_resume_without_a_journal(prepared, monkeypatch):
    record, _ = prepared
    started, _ = _stub_native(monkeypatch)
    with pytest.raises(FileNotFoundError):
        execution.execute_clustering_sensitivity(record, resume=True)
    assert started == []


def test_existing_results_are_never_reused_without_explicit_resume(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, _ = _stub_native(monkeypatch)
    assert execution.execute_clustering_sensitivity(record)["status"] == "completed"
    original = _immutable_bytes(root)
    with pytest.raises(ValueError, match="already executed"):
        execution.execute_clustering_sensitivity(record)
    assert started == BRANCHES
    _assert_preserved(original)
