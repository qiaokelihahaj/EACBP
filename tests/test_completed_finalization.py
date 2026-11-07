"""Presentation recovery must not compute, re-audit, or forgive missing science."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

import numpy as np
import pandas as pd
import pytest

from eacbp.application import finalization
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport, ValidationSeverity
from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.sc_data import SCData
from eacbp.evidence.snapshot import load_study_snapshot
from eacbp.orchestrator.checkpoint import StudyJournal, fingerprint
from eacbp.orchestrator.planning import resolve_task_contract
from eacbp.orchestrator.resume import ResumeManager
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.runtime import ExecutionState
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, default=str), encoding="utf-8")


@pytest.fixture
def completed(tmp_path, monkeypatch):
    run = tmp_path / "runs" / "completed"
    run.mkdir(parents=True)
    manifest = StudyManifest.model_validate({"study_id": "recovery", "title": "Describe observed data",
        "biological_design": {"species": "human", "tissue": "PBMC"},
        "data": {"raw_artifact_uri": "adata://recovery/raw/v1"}})
    config = {"method_profile": "baseline", "threads": 1}
    source_hash = hashlib.sha256(b"original public counts").hexdigest()
    source = {"path": "/input/data.h5ad", "sha256": source_hash,
              "sha256_before_import": source_hash, "sha256_after_import": source_hash}
    registry = ArtifactRegistry(str(run / "artifacts"))
    registry.register(manifest.data.raw_artifact_uri,
        SCData(X=np.array([[1, 2], [3, 4]], dtype=np.float32), obs=pd.DataFrame(index=["a", "b"]),
               var=pd.DataFrame(index=["g1", "g2"])), ArtifactType.ANNDATA, manifest.study_id,
        "task_000_ingest", "import", parameters={"source_sha256": source_hash})
    task = TaskContract(task_id="task_001_audit", capability="dataset_audit",
                        input_artifacts=[manifest.data.raw_artifact_uri], expected_outputs=["table://recovery/dataset_audit/v1"])
    capabilities = create_default_capability_registry()
    resolved = resolve_task_contract(task.model_copy(deep=True), manifest, config, capabilities, registry)
    files = {"core.py": "a" * 64}
    environment = {"python": "original-python", "platform": "original-platform", "distributions": {"numpy": "original"},
                   "source_sha256": fingerprint(files)}
    proof = tmp_path / "installed-source.json"
    _write(proof, {"status": "passed", "source_matches_installed_python_files": True, "n_python_files": 1,
                   "environment": environment, "source_files_sha256": files})
    signature = ResumeManager.signature(resolved, {manifest.data.raw_artifact_uri: registry.get_metadata(manifest.data.raw_artifact_uri).sha256_hash}, manifest, environment)
    result = TaskResult(task_id=task.task_id, capability=task.capability, method_used="sc_audit_v1", status=TaskStatus.SUCCESS,
                        input_artifacts=task.input_artifacts, output_artifacts=task.expected_outputs,
                        metrics={"n_cells": 2, "min_replicates": 3, "biological_replication_sufficient": True})
    with TaskArtifactTransaction(registry) as transaction:
        transaction.register(task.expected_outputs[0], pd.DataFrame({"field": ["n_cells"], "value": [2]}),
                             ArtifactType.TABLE, manifest.study_id, task.task_id, "audit_metadata", parent_uris=task.input_artifacts)
        registry.commit_task(transaction, signature, result)
    report = ValidationReport(auditor_name="original", target_task_id=task.task_id, overall_passed=True)
    report.add_check("clustering_separation_silhouette", False, ValidationSeverity.WARNING, "A retained warning")
    report.add_check("canonical_marker_representation", False, ValidationSeverity.INFO, "A retained diagnostic")
    registry.begin_audit(signature, resolved, result, auditor_name="original", auditor_version="1", auditor_fingerprint="original-auditor")
    registry.record_audit(signature, resolved, result, report, auditor_name="original", auditor_version="1", auditor_fingerprint="original-auditor")
    _, decisions = finalization.ComputationalDAGPlanner.adapt_after_audit([task.model_copy(deep=True)],
                                ExecutionState.from_input(config).preview_with_metrics(result.metrics))
    result.metrics["planning_decisions"] = decisions
    with StudyJournal(registry.storage.base_dir, manifest.study_id) as journal:
        journal.entries[task.task_id] = ResumeManager(registry).journal_entry(signature=signature, result=result,
                                                                            phase="audited", environment=environment)
        journal.save()
    journal_path = run / "artifacts/_runs" / f"{fingerprint(manifest.study_id)}.json"
    run_id = uuid4().hex
    history = {"study_id": manifest.study_id, "status": "success", "tasks_executed": 1, "tasks_blocked": 0,
               "run_id": run_id, "artifacts_created": 2, "evidence_nodes_count": 0, "claims_count": 0,
               "claims": [], "planning_decisions": decisions, "failures": [], "is_simulated": False}
    envelope = {"schema_version": 1, "status": "failed", "import_completed": True, "last_run_id": run_id,
                "manifest": manifest.model_dump(mode="json"), "config": config, "source": source,
                "error_type": "SnapshotValidationError", "error": "presentation serialization failed"}
    failure = {"status": "failed", "error_type": "SnapshotValidationError", "error": "presentation serialization failed"}
    for name, value in (("manifest.json", envelope["manifest"]), ("config.json", config), ("source.json", source),
                        ("run_config.json", envelope), ("summary.json", failure), ("error.json", failure),
                        (f"summaries/{run_id}.json", history)):
        _write(run / name, value)
    image = "sha256:" + "b" * 64
    invocation_path = run.parent / f"{run.name}-execution-command.json"
    _write(invocation_path, {"answer_reference_mounted": False, "command": ["docker", "run", image, "analyze",
            "--question", manifest.title, "--data", source["path"], "--species", "human", "--tissue", "PBMC",
            "--threads", "1", "--run-dir", "/results/completed"]})
    event_path = run / "artifacts/_runs/events" / journal_path.stem / f"{run_id}.jsonl"
    events = [("run_started", None, {}), ("plan_created", None, {"task_count": 1}),
              ("audit_finished", task.task_id, {"signature": signature, "duration_seconds": .5, "reason": None,
                                               "checks_failed": [check.check_name for check in report.checks if not check.passed]}),
              ("task_finished", task.task_id, {"status": "success", "method": result.method_used}),
              ("run_finished", None, {"status": "success", "tasks_executed": 1, "tasks_blocked": 0})]
    event_path.parent.mkdir(parents=True)
    event_path.write_text("\n".join(json.dumps({"sequence": i, "run_id": run_id, "kind": kind,
                          "task_id": task_id, "details": details}) for i, (kind, task_id, details) in enumerate(events, 1)), encoding="utf-8")
    monkeypatch.setattr(finalization, "build_study_tasks", lambda *args: [task.model_copy(deep=True)])
    monkeypatch.setattr(finalization, "execution_environment", lambda: {"source_sha256": "different-finalizer-source", "python": "current-python"})
    return SimpleNamespace(run=run, registry=registry, task=resolved, result=result, journal=journal_path, signature=signature,
                           proof=proof, proof_hash=finalization._hash(proof), image=image, events=event_path, invocation=invocation_path,
                           failure_bytes=(run / "summary.json").read_bytes())


def _finalize(case):
    return finalization.finalize_completed_journal(run_dir=case.run, compute_image=case.image,
                                                   source_proof=case.proof, source_proof_sha256=case.proof_hash)


def test_completed_recovery_preserves_science_and_calls_no_executor_or_auditor(completed, monkeypatch):
    from eacbp.orchestrator.execution import TaskExecutor
    from eacbp.auditor import ScientificAuditor
    def forbidden(*args, **kwargs):
        pytest.fail("Presentation finalization invoked computation or scientific mutation")
    monkeypatch.setattr(finalization.ScientificOrchestrator, "run_study", forbidden)
    monkeypatch.setattr(finalization.ScientificOrchestrator, "__init__", forbidden)
    monkeypatch.setattr(TaskExecutor, "execute", forbidden)
    monkeypatch.setattr(ScientificAuditor, "audit_task", forbidden)
    for method in ("register", "commit_task", "begin_audit", "record_audit"):
        monkeypatch.setattr(ArtifactRegistry, method, forbidden)
    protected = {str(path): finalization._hash(path) for path in completed.run.rglob("*") if path.is_file() and path.name != "summary.json"}
    summary = _finalize(completed)
    assert summary["status"] == "success" and summary["run_id"] == json.loads((completed.run / "run_config.json").read_text())["last_run_id"]
    assert all(finalization._hash(path) == digest for path, digest in protected.items())
    snapshot = load_study_snapshot(completed.run / "snapshot.json")
    assert snapshot.summary == summary == json.loads((completed.run / "summary.json").read_text())
    assert len(snapshot.task_results) == 1 and snapshot.evidence_graph.evidence_nodes
    binding = summary["finalization"]
    receipt_path = completed.run / binding["receipt"]
    receipt = json.loads(receipt_path.read_text())
    assert finalization._hash(receipt_path) == binding["receipt_sha256"]
    assert receipt["original_compute_environment"]["python"] == "original-python"
    assert receipt["current_finalizer_environment"]["python"] == "current-python"
    assert receipt["scientific_execution_calls"] == receipt["scientific_audit_calls"] == 0
    assert (receipt_path.parent / "previous-summary.json").read_bytes() == completed.failure_bytes
    assert (receipt_path.parent / "previous-run-config.json").read_bytes() == (completed.run / "run_config.json").read_bytes()
    assert (completed.run / "report.md").is_file()
    assert [(check.check_name, check.severity.value, check.passed) for check in snapshot.audit_reports[0].checks] == [
        ("clustering_separation_silhouette", "warning", False), ("canonical_marker_representation", "info", False)]


def test_complete_compute_journal_can_publish_partial_scientific_request(completed, monkeypatch):
    from eacbp.application.analysis_requirements import freeze_question_requirements

    envelope = json.loads((completed.run / "run_config.json").read_text())
    question = envelope["manifest"]["title"]
    source_hash = envelope["source"]["sha256"]
    config = envelope["config"]
    config["autonomous_analysis"] = {"question": question, "source_sha256": source_hash,
        "scientific_gaps": [], "requirement_contract": freeze_question_requirements(
            question=question, observations={}, source_sha256=source_hash)}
    _write(completed.run / "config.json", config)
    _write(completed.run / "run_config.json", envelope)
    protected = {str(path): finalization._hash(path) for path in completed.run.rglob("*")
                 if path.is_file() and path.name != "summary.json"}

    def forbidden(*args, **kwargs):
        pytest.fail("Partial presentation recovery must not compute or re-audit")
    from eacbp.orchestrator.execution import TaskExecutor
    from eacbp.auditor import ScientificAuditor
    monkeypatch.setattr(TaskExecutor, "execute", forbidden)
    monkeypatch.setattr(ScientificAuditor, "audit_task", forbidden)
    summary = _finalize(completed)
    assert summary["status"] == "success" and summary["analysis_request_status"] == "partial"
    assert all(finalization._hash(path) == digest for path, digest in protected.items())
    snapshot = load_study_snapshot(completed.run / "snapshot.json")
    assert snapshot.summary["requirement_completion_checks"] == summary["requirement_completion_checks"]
    assert "Requested analysis: `partial`" in (completed.run / "report.md").read_text(encoding="utf-8")
    receipt = json.loads((completed.run / summary["finalization"]["receipt"]).read_text())
    assert receipt["analysis_request_status"] == "partial"
    assert receipt["scientific_execution_calls"] == receipt["scientific_audit_calls"] == 0


@pytest.mark.parametrize("fault", ["missing_check", "extra_check", "reordered_checks", "signature", "reason", "error", "stop_rule"])
def test_original_warning_events_require_exact_admitted_receipt_match(completed, fault):
    events = [json.loads(line) for line in completed.events.read_text().splitlines()]
    details = events[2]["details"]
    if fault == "missing_check":
        details["checks_failed"].pop()
    elif fault == "extra_check":
        details["checks_failed"].append("unrecorded_check")
    elif fault == "reordered_checks":
        details["checks_failed"].reverse()
    elif fault == "signature":
        details["signature"] = "another-task-signature"
    elif fault == "reason":
        details["reason"] = "Original admission rejected"
    else:
        index_path = completed.run / "artifacts/.artifact_registry.json"
        index = json.loads(index_path.read_text())
        # Even inconsistent global passed flags cannot forgive a rejecting
        # original check severity.
        index["audit_records"][completed.signature]["report"]["checks"][0]["severity"] = fault
        _write(index_path, index)
    completed.events.write_text("\n".join(map(json.dumps, events)), encoding="utf-8")
    with pytest.raises(ValueError, match="Original audit event contradicts its saved admission"):
        _finalize(completed)
    assert (completed.run / "summary.json").read_bytes() == completed.failure_bytes
    assert not (completed.run / "snapshot.json").exists()


@pytest.mark.parametrize("fault", ["missing_task", "phase", "status", "environment", "signature", "result_metric", "output_hash",
                                   "contract", "audit_pending", "audit_stop", "payload", "config", "manifest", "source",
                                   "proof_hash", "proof_source", "compute_image", "event_coverage"])
def test_incomplete_or_tampered_original_science_is_never_finalized(completed, fault):
    case = completed
    journal = json.loads(case.journal.read_text())
    index_path = case.run / "artifacts/.artifact_registry.json"
    index = json.loads(index_path.read_text())
    if fault == "missing_task":
        journal.clear()
    elif fault == "phase":
        journal[case.task.task_id]["phase"] = "computed"
    elif fault == "status":
        journal[case.task.task_id]["result"]["status"] = "execution_failure"
    elif fault == "environment":
        journal[case.task.task_id]["environment"]["python"] = "another-compute-env"
    elif fault == "signature":
        journal[case.task.task_id]["signature"] = "forged"
    elif fault == "result_metric":
        journal[case.task.task_id]["result"]["metrics"]["n_cells"] = 999
    elif fault == "output_hash":
        journal[case.task.task_id]["output_hashes"].clear()
    elif fault == "contract":
        index["audit_records"][case.signature]["contract_fingerprint"] = "forged"
    elif fault == "audit_pending":
        index["audit_records"][case.signature]["status"] = "pending"
    elif fault == "audit_stop":
        index["audit_records"][case.signature]["stop_rule_triggered"] = True
    elif fault == "payload":
        meta = case.registry.get_metadata(case.task.expected_outputs[0])
        path = case.registry.storage._validated_path(Path(meta.storage_path))
        path.write_bytes(path.read_bytes() + b"tampered")
    elif fault in {"config", "manifest", "source"}:
        path = case.run / f"{fault}.json"
        value = json.loads(path.read_text())
        value["changed"] = True
        if fault == "config":
            value["threads"] = 99
        _write(path, value)
    elif fault == "proof_hash":
        case.proof_hash = "0" * 64
    elif fault == "proof_source":
        proof = json.loads(case.proof.read_text())
        proof["source_files_sha256"]["core.py"] = "c" * 64
        _write(case.proof, proof)
        case.proof_hash = finalization._hash(case.proof)
    elif fault == "compute_image":
        case.image = "sha256:" + "c" * 64
    else:
        events = [json.loads(line) for line in case.events.read_text().splitlines()]
        events[2]["task_id"] = "unknown"
        case.events.write_text("\n".join(map(json.dumps, events)), encoding="utf-8")
    if fault in {"missing_task", "phase", "status", "environment", "signature", "result_metric", "output_hash"}:
        _write(case.journal, journal)
    if fault in {"contract", "audit_pending", "audit_stop"}:
        _write(index_path, index)
    with pytest.raises(Exception):
        _finalize(case)
    assert (case.run / "summary.json").read_bytes() == case.failure_bytes
    assert not (case.run / "snapshot.json").exists()
    assert not (case.run / "report.md").exists()


def test_current_scientific_contract_change_is_rejected(completed, monkeypatch):
    changed = completed.task.model_copy(deep=True)
    changed.parameters["a_new_scientific_choice"] = True
    monkeypatch.setattr(finalization, "build_study_tasks", lambda *args: [changed.model_copy(deep=True)])
    with pytest.raises(ValueError, match="signature changed"):
        _finalize(completed)
    assert (completed.run / "summary.json").read_bytes() == completed.failure_bytes


def test_cli_finalization_is_a_distinct_presentation_entrypoint(completed, monkeypatch, capsys):
    from eacbp.cli import main
    assert main(["finalize-completed-journal", "--run-dir", str(completed.run), "--compute-image", completed.image,
                 "--source-proof", str(completed.proof), "--source-proof-sha256", completed.proof_hash]) == 0
    assert '"presentation_only_completed_journal"' in capsys.readouterr().out


@pytest.mark.parametrize("changes", [{"skipped": True}, {"all_requested_fits_completed": False},
                                   {"n_successful_fits": 7}, {"n_failed_fits": 1}, {"complete_fit_coverage": .875},
                                   {"n_requested_fits": None, "n_successful_fits": None},
                                   {"n_requested_fits": 0, "n_successful_fits": 0}])
def test_sensitivity_coverage_cannot_be_inferred_from_a_success_status(changes):
    original = {"all_requested_fits_completed": True, "n_requested_fits": 8, "n_successful_fits": 8,
                "n_failed_fits": 0, "complete_fit_coverage": 1}
    assert finalization._completed_sensitivity(original)
    assert not finalization._completed_sensitivity({**original, **changes})
