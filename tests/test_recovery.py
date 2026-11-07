"""Inspection is read-only; continuation still requires current admission."""

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from eacbp.application import clustering_sensitivity_execution as execution
from eacbp.application import execution_journal as journal_module
from eacbp.application import recovery
from eacbp.application.execution_journal import ExecutionJournal
from eacbp.application import study_service

from test_clustering_sensitivity_execution import _read, _source
from test_clustering_sensitivity_recovery import BRANCHES, ENVIRONMENT, _stub_native as _original_stub_native


@pytest.fixture
def prepared(tmp_path, monkeypatch):
    monkeypatch.setattr(execution, "execution_environment", lambda: deepcopy(ENVIRONMENT))
    source, config, manifest = _source(tmp_path)
    record = execution.prepare_clustering_sensitivity(manifest=manifest, config=config, data=source,
                                                     output_dir=tmp_path / "experiment")
    return record, source


def _identity(record):
    bindings = _read(record["bindings_path"])
    result = {key: bindings[key] for key in ("plan_file_sha256", "plan_sha256", "source_sha256", "execution_environment")}
    result["bindings_sha256"] = record["bindings_sha256"]
    return result


def _initialize(record, *, stage="invocation", finish=None):
    with ExecutionJournal(record["output_dir"], identity=_identity(record)) as journal:
        attempt = journal.start_attempt(stage, "primary" if stage == "native_execution" else None)
        if finish:
            journal.finish_attempt(attempt, status=finish)
    return attempt


def _files(root, *, exclude_lock=False):
    return {path.relative_to(root).as_posix(): path.read_bytes() for path in root.rglob("*") if path.is_file()
            and not (exclude_lock and path.name == ".execution-journal.lock")}


def _inspect(record, **kwargs):
    return recovery.inspect_recovery(record, expected_bindings_sha256=record["bindings_sha256"], **kwargs)


def _stub_native(monkeypatch, **kwargs):
    started, collected = _original_stub_native(monkeypatch, **kwargs)
    run, collect = study_service.run_study, execution.collect_clustering_run

    def write_payload(**arguments):
        try:
            return run(**arguments)
        finally:
            native = Path(arguments["run_dir"])
            if native.is_dir():
                (native / "artifacts").mkdir(exist_ok=True)
                (native / "artifacts" / "fixture.bin").write_bytes(b"native-artifact-fixture")

    def collect_payload(path, **arguments):
        receipt = collect(path, **arguments)
        payload = Path(path) / "artifacts" / "fixture.bin"
        receipt["provenance"]["artifacts"] = [{"path": str(payload), "sha256": "sha256:" + study_service._hash_file(payload)}]
        return receipt

    monkeypatch.setattr(study_service, "run_study", write_payload)
    monkeypatch.setattr(execution, "collect_clustering_run", collect_payload)
    return started, collected


@pytest.mark.parametrize("pin", [None, "", "a" * 63, "A" * 64, "g" * 64, True])
def test_strong_pin_is_required_even_when_prepared_record_contains_its_own_pin(prepared, pin, monkeypatch):
    record, _ = prepared
    monkeypatch.setattr(execution, "execute_clustering_sensitivity", lambda *args, **kwargs: pytest.fail("Blocked inspection executed"))
    before = _files(Path(record["output_dir"]))
    result = recovery.resume_recovery(record, expected_bindings_sha256=pin)
    assert result["status"] == "blocked" and result["can_resume"] is False
    assert result["blockers"][0]["code"] == "bindings_sha256"
    assert _files(Path(record["output_dir"])) == before


def test_a_fresh_preparation_cannot_implicitly_create_a_recovery_journal(prepared):
    record, _ = prepared
    root = Path(record["output_dir"])
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked" and result["blockers"][0]["code"] == "executor_lock"
    assert not (root / "execution_journal").exists()
    assert not (root / ".execution-journal.lock").exists()
    assert _files(root) == before


def test_missing_preparation_does_not_create_a_directory(tmp_path):
    target = tmp_path / "missing" / "experiment"
    result = recovery.inspect_recovery(target, expected_bindings_sha256="a" * 64)
    assert result["status"] == "blocked"
    assert not target.parent.exists()


@pytest.mark.parametrize("kind", ["record", "directory", "bindings_file"])
def test_inspect_preserves_unfinished_events_without_imagined_interruption_or_elapsed(prepared, monkeypatch, kind):
    record, _ = prepared
    attempt = _initialize(record, stage="native_execution")
    root = Path(record["output_dir"])
    before = _files(root)
    monkeypatch.setattr(execution, "collect_clustering_run", lambda *args, **kwargs: pytest.fail("Inspect collected native outputs"))
    monkeypatch.setattr(study_service, "_make_registry", lambda *args, **kwargs: pytest.fail("Inspect opened an artifact registry"))
    monkeypatch.setattr(study_service, "_import_h5ad", lambda *args, **kwargs: pytest.fail("Inspect imported a matrix"))
    prepared_arg = {"record": record, "directory": root, "bindings_file": root / "bindings.json"}[kind]
    result = recovery.inspect_recovery(prepared_arg, expected_bindings_sha256=record["bindings_sha256"])
    assert result["status"] == "ready" and result["can_resume"] is True
    pending = result["unfinished_attempts"]
    assert len(pending) == 1 and pending[0]["attempt_id"] == attempt
    assert pending[0]["status"] == "started" and pending[0]["finish"] is None
    assert result["failed_attempts"] == []
    assert result["reuse_authorized"] is False and result["requires_native_readmission"] is True
    assert all(branch["reuse_authorized"] is False for branch in result["branches"])
    assert result["identity"]["verified"] is True
    assert _files(root) == before


def test_valid_event_ahead_of_head_is_reported_without_repair(prepared):
    record, _ = prepared
    root = Path(record["output_dir"])
    with ExecutionJournal(root, identity=_identity(record)) as journal:
        initial_head = (journal.path / "head.json").read_bytes()
        journal.start_attempt("invocation")
    (root / "execution_journal/head.json").write_bytes(initial_head)
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "ready"
    assert result["journal"]["head_needs_repair"] is True
    assert result["journal"]["head_sequence"] == 0 and result["journal"]["event_count"] == 1
    assert _files(root) == before
    # Default resume loading still repairs the head and closes the old start.
    with ExecutionJournal(root, identity=_identity(record), resume=True) as journal:
        assert journal.attempts[0]["status"] == "interrupted"
        assert journal.attempts[0]["finish"]["elapsed_seconds"] is None
    assert _read(root / "execution_journal/head.json")["sequence"] == 2


@pytest.mark.parametrize("target", ["bindings", "plan", "source", "derived_source", "config", "manifest", "environment"])
def test_frozen_identity_change_blocks_before_any_continuation(prepared, monkeypatch, target):
    record, source = prepared
    _initialize(record)
    bindings = _read(record["bindings_path"])
    if target == "environment":
        monkeypatch.setattr(execution, "execution_environment", lambda: {**ENVIRONMENT, "source_sha256": "f" * 64})
    else:
        path = {"bindings": Path(record["bindings_path"]), "plan": Path(record["plan_path"]), "source": source,
                "derived_source": Path(bindings["branches"][-1]["data_path"]),
                "config": Path(bindings["branches"][1]["config_path"]),
                "manifest": Path(bindings["branches"][1]["manifest_path"])}[target]
        path.write_bytes(path.read_bytes() + b"\n")
    monkeypatch.setattr(execution, "execute_clustering_sensitivity", lambda *args, **kwargs: pytest.fail("Changed identity executed"))
    before = _files(Path(record["output_dir"]))
    result = recovery.resume_recovery(record, expected_bindings_sha256=record["bindings_sha256"])
    assert result["status"] == "blocked" and result["identity"]["verified"] is False
    assert _files(Path(record["output_dir"])) == before


@pytest.mark.parametrize("target", ["identity", "event", "missing_event", "head", "missing_lock"])
def test_journal_damage_is_blocked_without_repair(prepared, target):
    record, _ = prepared
    _initialize(record, finish="failed")
    root = Path(record["output_dir"])
    path = {"identity": root / "execution_journal/identity.json", "event": root / "execution_journal/events/00000001.json",
            "missing_event": root / "execution_journal/events/00000001.json", "head": root / "execution_journal/head.json",
            "missing_lock": root / ".execution-journal.lock"}[target]
    if target.startswith("missing"):
        path.unlink()
    else:
        value = _read(path)
        value["unexpected_change"] = True
        path.write_text(json.dumps(value), encoding="utf-8")
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked"
    assert _files(root) == before


def test_consistently_rehashed_false_elapsed_is_rejected(prepared):
    record, _ = prepared
    _initialize(record, finish="failed")
    root = Path(record["output_dir"])
    path = root / "execution_journal/events/00000002.json"
    event = _read(path)
    event["elapsed_seconds"] += 10
    event.pop("event_sha256")
    event["event_sha256"] = journal_module._fingerprint(event)
    path.write_text(json.dumps(event), encoding="utf-8")
    head = _read(root / "execution_journal/head.json")
    head["event_sha256"] = event["event_sha256"]
    (root / "execution_journal/head.json").write_text(json.dumps(head), encoding="utf-8")
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked" and "elapsed" in result["blockers"][0]["reason"]
    assert _files(root) == before


def test_concurrent_executor_is_rejected_without_touching_events(prepared):
    record, _ = prepared
    root = Path(record["output_dir"])
    with ExecutionJournal(root, identity=_identity(record)) as journal:
        journal.start_attempt("invocation")
        before = _files(root, exclude_lock=True)
        result = _inspect(record)
        assert result["status"] == "blocked" and result["blockers"][0]["code"] == "executor_lock"
        assert result["blockers"][0]["error_type"] == "RuntimeError"
        assert _files(root, exclude_lock=True) == before


def test_actual_child_process_rejects_the_live_parent_executor_lock(prepared):
    record, _ = prepared
    root = Path(record["output_dir"])
    program = """
import json, sys
from eacbp.application.recovery import inspect_recovery
print(json.dumps(inspect_recovery(json.loads(sys.argv[1]), expected_bindings_sha256=sys.argv[2])))
"""
    with ExecutionJournal(root, identity=_identity(record)) as journal:
        journal.start_attempt("invocation")
        before = _files(root, exclude_lock=True)
        child = subprocess.run([sys.executable, "-B", "-c", program, json.dumps(record), record["bindings_sha256"]],
                               capture_output=True, text=True, timeout=30)
        assert child.returncode == 0, child.stderr
        result = json.loads(child.stdout)
        assert result["status"] == "blocked" and result["blockers"][0]["code"] == "executor_lock"
        assert result["blockers"][0]["error_type"] == "RuntimeError"
        assert "active executor" in result["blockers"][0]["reason"]
        assert _files(root, exclude_lock=True) == before


def _admitted_fixture(record):
    root = Path(record["output_dir"])
    binding = _read(record["bindings_path"])["branches"][0]
    native = Path(binding["run_dir"])
    native.mkdir(parents=True)
    (native / "summary.json").write_text('{"status":"success"}', encoding="utf-8")
    payload = native / "artifacts" / "payload.bin"
    payload.parent.mkdir()
    payload.write_bytes(b"original-scientific-values")
    with ExecutionJournal(root, identity=_identity(record)) as journal:
        attempt = journal.start_attempt("native_execution", "primary")
        journal.finish_attempt(attempt, status="completed", details={"admitted": True,
            "provenance": {"run_dir": str(native), "files": execution._native_file_hashes(native),
                           "artifacts": [{"path": str(payload), "sha256": "sha256:" + study_service._hash_file(payload)}]}})
    return native, payload


@pytest.mark.parametrize("target", ["metadata", "extra_metadata", "payload", "missing_payload", "missing_native"])
def test_previous_admission_metadata_and_binary_payload_pins_are_checked(prepared, target):
    record, _ = prepared
    native, payload = _admitted_fixture(record)
    if target == "metadata":
        (native / "summary.json").write_text('{"status":"failed"}', encoding="utf-8")
    elif target == "extra_metadata":
        (native / "invented.json").write_text("{}", encoding="utf-8")
    elif target == "payload":
        payload.write_bytes(b"changed-scientific-values")
    elif target == "missing_payload":
        payload.unlink()
    else:
        native.rename(native.with_name("preserved-old-primary"))
    before = _files(Path(record["output_dir"]))
    result = _inspect(record)
    assert result["status"] == "blocked" and result["blockers"][0]["code"] == "native_evidence_integrity"
    assert _files(Path(record["output_dir"])) == before


def test_completed_admission_pins_never_authorize_reuse_at_inspection(prepared):
    record, _ = prepared
    _admitted_fixture(record)
    result = _inspect(record)
    branch = result["branches"][0]
    assert result["status"] == "ready" and branch["state"] == "pending_re_admission"
    assert branch["prior_pin_available"] is True and branch["native_metadata_pins_verified"] is True
    assert branch["native_payload_pins_verified"] == 1
    assert branch["reuse_authorized"] is False and branch["requires_native_readmission"] is True


def _serialized_admission_fixture(record, payload_type):
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.schemas.artifact import ArtifactType
    import pandas as pd

    root = Path(record["output_dir"])
    native = root / "runs" / "primary"
    native.mkdir(parents=True)
    registry = ArtifactRegistry(str(native / "artifacts"))
    if payload_type == "table":
        payload = pd.DataFrame({"mean_raw_count": [36775107 / 2933, float("nan")]}, index=["gene_1", "gene_2"])
        artifact_type = ArtifactType.TABLE
    else:
        payload = {"fixture_values": [1, 2, 3]}
        artifact_type = ArtifactType.JSON
    metadata = registry.register(f"{payload_type}://fixture/published/v1", payload, artifact_type,
                                 "fixture", "task_fixture", "serializer_format_regression")
    assert metadata.sha256_hash == "sha256:" + study_service._hash_file(Path(metadata.storage_path))
    (native / "summary.json").write_text('{"status":"success"}', encoding="utf-8")
    with ExecutionJournal(root, identity=_identity(record)) as journal:
        attempt = journal.start_attempt("native_execution", "primary")
        journal.finish_attempt(attempt, status="completed", details={"admitted": True,
            "provenance": {"run_dir": str(native), "files": execution._native_file_hashes(native),
                           "artifacts": [{"uri": metadata.uri, "path": metadata.storage_path,
                                          "sha256": metadata.sha256_hash}]}})
    return Path(metadata.storage_path), metadata.sha256_hash


@pytest.mark.parametrize("payload_type", ["table", "json"])
@pytest.mark.parametrize("changed", [False, True])
def test_actual_serializer_tagged_payload_is_checked_without_registry_or_matrix_loading(prepared, monkeypatch, payload_type, changed):
    record, _ = prepared
    payload, tagged_digest = _serialized_admission_fixture(record, payload_type)
    if changed:
        payload.write_bytes(payload.read_bytes() + b"\nchanged-after-admission\n")
    monkeypatch.setattr(execution, "collect_clustering_run", lambda *args, **kwargs: pytest.fail("Inspection performed native admission"))
    monkeypatch.setattr(study_service, "_make_registry", lambda *args, **kwargs: pytest.fail("Inspection opened a registry"))
    monkeypatch.setattr(study_service, "_import_h5ad", lambda *args, **kwargs: pytest.fail("Inspection loaded a matrix"))
    root = Path(record["output_dir"])
    before = _files(root)
    result = _inspect(record)
    assert _files(root) == before
    if changed:
        assert result["status"] == "blocked" and result["blockers"][0]["code"] == "native_evidence_integrity"
    else:
        assert result["status"] == "ready"
        assert result["branches"][0]["native_payload_pins_verified"] == 1
        assert result["branches"][0]["reuse_authorized"] is False
        native = result["journal"]["attempts"][0]
        assert native["finish"]["details"]["provenance"]["artifacts"][0]["sha256"] == tagged_digest


def _rewrite_finished_pin(root, mutate):
    event_path = root / "execution_journal/events/00000002.json"
    event = _read(event_path)
    mutate(event["details"]["provenance"])
    event.pop("event_sha256")
    event["event_sha256"] = journal_module._fingerprint(event)
    event_path.write_text(json.dumps(event), encoding="utf-8")
    head = _read(root / "execution_journal/head.json")
    head["event_sha256"] = event["event_sha256"]
    (root / "execution_journal/head.json").write_text(json.dumps(head), encoding="utf-8")


def test_legacy_bare_payload_digest_remains_compatible(prepared):
    record, _ = prepared
    _, tagged_digest = _serialized_admission_fixture(record, "table")
    root = Path(record["output_dir"])
    _rewrite_finished_pin(root, lambda provenance: provenance["artifacts"][0].update(sha256=tagged_digest[len("sha256:"):]))
    before = _files(root)
    assert _inspect(record)["status"] == "ready"
    assert _files(root) == before


@pytest.mark.parametrize("digest", ["md5:" + "a" * 64, "sha256:sha256:" + "a" * 64,
                                   "SHA256:" + "a" * 64, "sha256:" + "A" * 64, "sha256:" + "a" * 63])
def test_algorithm_tag_does_not_admit_unknown_or_malformed_payload_checksums(prepared, digest):
    record, _ = prepared
    _admitted_fixture(record)
    root = Path(record["output_dir"])
    _rewrite_finished_pin(root, lambda provenance: provenance["artifacts"][0].update(sha256=digest))
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked" and result["blockers"][0]["code"] == "native_evidence_integrity"
    assert "invalid file SHA pin" in result["blockers"][0]["reason"]
    assert _files(root) == before


def test_native_metadata_file_map_does_not_accept_artifact_algorithm_tags(prepared):
    record, _ = prepared
    _admitted_fixture(record)
    root = Path(record["output_dir"])

    def tag_file_pin(provenance):
        path = next(iter(provenance["files"]))
        provenance["files"][path] = "sha256:" + provenance["files"][path]

    _rewrite_finished_pin(root, tag_file_pin)
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked" and result["blockers"][0]["code"] == "native_evidence_integrity"
    assert "invalid file SHA pin" in result["blockers"][0]["reason"]
    assert _files(root) == before


@pytest.mark.parametrize("damage", ["missing_artifacts", "empty_artifacts", "missing_files", "foreign_payload"])
def test_consistent_journal_with_incomplete_or_foreign_admission_pins_is_blocked(prepared, damage):
    record, _ = prepared
    native, payload = _admitted_fixture(record)
    root = Path(record["output_dir"])
    path = root / "execution_journal/events/00000002.json"
    event = _read(path)
    provenance = event["details"]["provenance"]
    if damage == "missing_artifacts":
        provenance.pop("artifacts")
    elif damage == "empty_artifacts":
        provenance["artifacts"] = []
    elif damage == "missing_files":
        provenance.pop("files")
    else:
        foreign = root / "foreign-payload.bin"
        foreign.write_bytes(payload.read_bytes())
        provenance["artifacts"][0]["path"] = str(foreign)
    event.pop("event_sha256")
    event["event_sha256"] = journal_module._fingerprint(event)
    path.write_text(json.dumps(event), encoding="utf-8")
    head = _read(root / "execution_journal/head.json")
    head["event_sha256"] = event["event_sha256"]
    (root / "execution_journal/head.json").write_text(json.dumps(head), encoding="utf-8")
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "blocked" and result["blockers"][0]["code"] == "native_evidence_integrity"
    assert _files(root) == before


def test_native_directory_without_prior_pin_is_pending_admission(prepared):
    record, _ = prepared
    _initialize(record, stage="native_execution")
    root = Path(record["output_dir"])
    native = root / "runs" / "primary"
    native.mkdir(parents=True)
    (native / "summary.json").write_text('{"status":"success"}', encoding="utf-8")
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "ready" and result["branches"][0]["state"] == "pending_re_admission"
    assert result["branches"][0]["prior_pin_available"] is False
    assert result["branches"][0]["reuse_authorized"] is False
    assert _files(root) == before


def test_postprocessing_failure_resumes_via_real_collector_without_new_fits(prepared, monkeypatch):
    record, _ = prepared
    root = Path(record["output_dir"])
    started, collected = _stub_native(monkeypatch)
    render = execution._markdown
    monkeypatch.setattr(execution, "_markdown", lambda result: (_ for _ in ()).throw(RuntimeError("Report rendering failed")))
    with pytest.raises(RuntimeError, match="rendering failed"):
        execution.execute_clustering_sensitivity(record)
    before = _files(root)
    inspected = _inspect(record)
    assert inspected["status"] == "ready"
    assert any(item["status"] == "postprocessing_failed" for item in inspected["failed_attempts"])
    assert inspected["latest_comparison"]["integrity_checked"] is False
    assert _files(root) == before
    monkeypatch.setattr(execution, "_markdown", render)
    result = recovery.resume_recovery(record, expected_bindings_sha256=record["bindings_sha256"])
    assert result["status"] == "completed" and result["native_reused_branch_ids"] == BRANCHES
    assert result["native_started_branch_ids"] == [] and started == BRANCHES and collected == BRANCHES * 2
    assert result["biological_stability_claim_supported"] is False
    assert all(path.exists() and path.read_bytes() == data for relative, data in before.items()
               if relative != "execution_journal/head.json" for path in [root / relative])


def test_crash_after_native_write_is_readmitted_before_other_branches_run(prepared, monkeypatch):
    record, _ = prepared
    started, collected = _stub_native(monkeypatch, interrupted="primary")
    with pytest.raises(SystemExit):
        execution.execute_clustering_sensitivity(record)
    result = _inspect(record)
    assert result["status"] == "ready" and result["branches"][0]["prior_pin_available"] is False
    assert result["branches"][0]["state"] == "pending_re_admission"
    assert any(item["stage"] == "native_execution" and item["status"] == "started" for item in result["unfinished_attempts"])
    resumed = recovery.resume_recovery(record, expected_bindings_sha256=record["bindings_sha256"])
    assert resumed["status"] == "completed" and resumed["native_reused_branch_ids"] == ["primary"]
    assert resumed["native_started_branch_ids"] == BRANCHES[1:]
    assert started == BRANCHES and collected == BRANCHES


def test_old_completed_journal_cannot_override_a_current_collector_rejection(prepared, monkeypatch):
    record, _ = prepared
    started, collected = _stub_native(monkeypatch)
    execution.execute_clustering_sensitivity(record)
    assert _inspect(record)["status"] == "ready"
    original_collect = execution.collect_clustering_run

    def reject_primary(path, *, plan, binding):
        receipt = original_collect(path, plan=plan, binding=binding)
        if binding["branch_id"] == "primary":
            receipt.update(admitted=False, reason="Current audit admission rejected")
        return receipt

    monkeypatch.setattr(execution, "collect_clustering_run", reject_primary)
    result = recovery.resume_recovery(record, expected_bindings_sha256=record["bindings_sha256"])
    assert result["status"] != "completed" and result["native_reused_branch_ids"] == []
    assert result["primary_native_summary"]["admitted"] is False
    assert result["native_started_branch_ids"] == [] and started == BRANCHES
    assert collected == BRANCHES + ["primary"]


def test_executor_still_rejects_a_lock_acquired_after_inspection(prepared, monkeypatch):
    record, _ = prepared
    _initialize(record, finish="failed")
    execute = execution.execute_clustering_sensitivity

    def compete(*args, **kwargs):
        with ExecutionJournal(record["output_dir"], identity=_identity(record), resume=True):
            return execute(*args, **kwargs)

    monkeypatch.setattr(execution, "execute_clustering_sensitivity", compete)
    with pytest.raises(RuntimeError, match="active executor"):
        recovery.resume_recovery(record, expected_bindings_sha256=record["bindings_sha256"])


def test_unbound_detached_report_is_not_a_native_reuse_gate(prepared, monkeypatch):
    record, _ = prepared
    _stub_native(monkeypatch)
    execution.execute_clustering_sensitivity(record)
    root = Path(record["output_dir"])
    (root / "comparison.json").write_text("An independently edited historical alias", encoding="utf-8")
    (root / "report.md").write_text("An independently edited detached report", encoding="utf-8")
    before = _files(root)
    result = _inspect(record)
    assert result["status"] == "ready" and result["latest_comparison"]["integrity_checked"] is False
    assert "Detached historical" in result["integrity_scope"]
    assert _files(root) == before
