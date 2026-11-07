"""Rebuild presentation from a complete, immutable, previously audited journal.

This is intentionally separate from resume: it never executes a task, runs an
auditor, updates scientific receipts, or changes the original compute environment.
Local source/image records provide consistency checks, not signed attestation.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from types import SimpleNamespace
from uuid import UUID, uuid4

from eacbp._atomic_json import atomic_write_json
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport, ValidationSeverity
from eacbp.capabilities import create_default_capability_registry
from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.graph import EvidenceGraph
from eacbp.orchestrator.checkpoint import StudyJournal, execution_environment, fingerprint
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.orchestrator.planning import build_study_tasks, resolve_task_contract
from eacbp.orchestrator.resume import ResumeManager
from eacbp.schemas.audit import stable_fingerprint
from eacbp.schemas.runtime import ExecutionState, RunConfig
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskResult, TaskStatus


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _object(path):
    value = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return value


def _write_json(path, value):
    return atomic_write_json(path, value, dump_kwargs={"allow_nan": False, "indent": 2})


def _write_report(path, markdown):
    descriptor, temporary = tempfile.mkstemp(prefix=".finalization-report-", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(markdown)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _original_environment(entries, proof_path, proof_sha256):
    if not re.fullmatch(r"[0-9a-f]{64}", proof_sha256) or _hash(proof_path) != proof_sha256:
        raise ValueError("Original installed-source proof hash does not match its independent pin")
    proof = _object(proof_path)
    files = proof.get("source_files_sha256", {})
    if (proof.get("status") != "passed" or proof.get("source_matches_installed_python_files") is not True
            or not files or proof.get("n_python_files") != len(files)
            or any(not isinstance(key, str) or Path(key).is_absolute() or ".." in Path(key).parts
                   or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest)
                   for key, digest in files.items())):
        raise ValueError("Original source proof is incomplete or invalid")
    environment = proof.get("environment")
    if not isinstance(environment, dict) or environment.get("source_sha256") != fingerprint(files):
        raise ValueError("Original source proof does not reconstruct its environment source fingerprint")
    if not entries or any(entry.get("environment") != environment for entry in entries.values()):
        raise ValueError("Every completed journal task must retain the exact original proved environment")
    return environment


def _original_invocation(run_dir, envelope, manifest, config, compute_image):
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", compute_image):
        raise ValueError("Original compute image must be an independently checked sha256 digest")
    path = run_dir.parent / f"{run_dir.name}-execution-command.json"
    receipt = _object(path)
    command = receipt.get("command", [])
    if not isinstance(command, list) or any(not isinstance(item, str) for item in command) or "analyze" not in command:
        raise ValueError("An original native analyze execution receipt is required")
    index = command.index("analyze")
    if index < 1 or command[index - 1] != compute_image or receipt.get("answer_reference_mounted") is not False:
        raise ValueError("Original compute image/reference declaration does not match the pinned invocation")
    tokens = command[index + 1:]
    if len(tokens) % 2 or len(set(tokens[::2])) != len(tokens[::2]):
        raise ValueError("Original analyze invocation has an unsupported shape")
    flags = dict(zip(tokens[::2], tokens[1::2]))
    request = config.get("autonomous_analysis", {})
    if (flags.get("--question") != request.get("question", manifest.title)
            or flags.get("--data") != envelope.get("source", {}).get("path")
            or flags.get("--species") != manifest.biological_design.species
            or flags.get("--tissue") != manifest.biological_design.tissue
            or flags.get("--threads") != str(config.get("threads", 1))
            or Path(flags.get("--run-dir", "")).name != run_dir.name):
        raise ValueError("Original invocation differs from the saved scientific question/source/run/config")
    return path


def _event_proof(run_dir, journal_path, run_id, entries, registry):
    path = run_dir / "artifacts" / "_runs" / "events" / journal_path.stem / f"{UUID(run_id).hex}.jsonl"
    events = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if (not events or [event.get("sequence") for event in events] != list(range(1, len(events) + 1))
            or any(event.get("run_id") != run_id for event in events)
            or events[-1].get("kind") != "run_finished"
            or events[-1].get("details", {}).get("status") != "success"
            or events[-1].get("details", {}).get("tasks_executed") != len(entries)
            or events[-1].get("details", {}).get("tasks_blocked") != 0
            or any(event.get("kind") in {"task_failed", "task_blocked", "audit_rejected", "run_failed"} for event in events)):
        raise ValueError("Original compute event log does not certify a complete successful native run")
    for kind in ("task_finished", "audit_finished"):
        selected = [event for event in events if event.get("kind") == kind]
        if len(selected) != len(entries) or {event.get("task_id") for event in selected} != set(entries):
            raise ValueError("Original event task/audit coverage differs from the completed journal")
        for event in selected:
            entry = entries[event["task_id"]]
            detail = event.get("details", {})
            if kind == "task_finished" and (detail.get("status") != "success" or detail.get("method") != entry["result"].get("method_used")):
                raise ValueError("Original task event contradicts its saved result")
            if kind == "audit_finished":
                audit = registry.get_audit_record(entry.get("signature"))
                report = ValidationReport.model_validate(audit.report) if audit is not None else None
                # Events retain every failed diagnostic, including admissible
                # warning/info checks. Compare them to the original receipt;
                # only rejecting severities prevent scientific admission.
                failed = [check for check in report.checks if not check.passed] if report is not None else []
                if (audit is None or audit.task_id != event["task_id"] or audit.status.value != "passed"
                        or not audit.overall_passed or audit.stop_rule_triggered
                        or report.target_task_id != event["task_id"] or not report.overall_passed or report.stop_rule_triggered
                        or any(check.severity in {ValidationSeverity.ERROR, ValidationSeverity.STOP_RULE} for check in failed)
                        or detail.get("signature") != entry.get("signature")
                        or detail.get("checks_failed") != [check.check_name for check in failed]
                        or detail.get("reason") != audit.rejection_reason):
                    raise ValueError("Original audit event contradicts its saved admission")
    return path


def _completed_sensitivity(metrics):
    requested = metrics.get("n_requested_fits")
    return (not metrics.get("skipped") and metrics.get("all_requested_fits_completed") is True
            and isinstance(requested, int) and not isinstance(requested, bool) and requested > 0
            and requested == metrics.get("n_successful_fits") and metrics.get("n_failed_fits") == 0
            and metrics.get("complete_fit_coverage") == 1)


def _validate_completed_plan(manifest, config, entries, registry, capabilities, environment):
    state = ExecutionState.from_input(config)
    if "method_profile" not in state.config:
        raw = registry.get_metadata(manifest.data.raw_artifact_uri)
        state.set_derived("method_profile", "baseline" if state.current_state.get("mode") == "demo" or raw.summary_metrics.get("is_simulated") else "standard")
    tasks = build_study_tasks(manifest, state.current_state, capabilities)
    records, results, reports = [], [], []
    index = 0
    while index < len(tasks):
        task = tasks[index]
        entry = entries.get(task.task_id)
        if not isinstance(entry, dict) or entry.get("phase") != "audited":
            raise ValueError(f"Planned task is missing completed audit phase: {task.task_id}")
        if any(state.completed.get(dep) != "success" for dep in task.depends_on):
            raise ValueError("Original task dependency coverage is incomplete")
        branch = task.parameters.get("target_branch")
        task = resolve_task_contract(task, manifest, state.state_for_branch(branch), capabilities, registry)
        result = TaskResult.model_validate(entry.get("result", {}))
        if (result.status != TaskStatus.SUCCESS or result.task_id != task.task_id or result.capability != task.capability
                or result.input_artifacts != task.input_artifacts or result.output_artifacts != task.expected_outputs
                or capabilities.get(task.capability, task.method).implementation_id != result.method_used):
            raise ValueError(f"Saved result does not cover its complete planned contract: {task.task_id}")
        input_hashes = {uri: meta.sha256_hash for uri, meta in registry.verify_many(task.input_artifacts).items()}
        signature = ResumeManager.signature(task, input_hashes, manifest, environment)
        if entry.get("signature") != signature:
            raise ValueError(f"Saved scientific configuration/contract signature changed: {task.task_id}")
        commit = registry.get_task_commit(signature)
        audit = registry.get_audit_record(signature, task)
        if not commit or audit is None or audit.status.value != "passed" or not audit.overall_passed or audit.stop_rule_triggered:
            raise ValueError(f"Original computation/audit admission is incomplete: {task.task_id}")
        if (audit.contract != registry._contract_payload(task) or audit.contract_fingerprint != stable_fingerprint(audit.contract)
                or audit.result_fingerprint != registry._result_identity(entry["result"])
                or audit.input_artifacts != result.input_artifacts or audit.output_artifacts != result.output_artifacts):
            raise ValueError(f"Original audit contract/result fingerprints do not match: {task.task_id}")
        saved_result, committed_result = deepcopy(entry["result"]), deepcopy(commit.get("result", {}))
        if task.capability == "dataset_audit":
            saved_result.get("metrics", {}).pop("planning_decisions", None)
            committed_result.get("metrics", {}).pop("planning_decisions", None)
        if stable_fingerprint(saved_result) != stable_fingerprint(committed_result):
            raise ValueError(f"Full scientific journal result differs from its immutable computation receipt: {task.task_id}")
        hashes = {uri: meta.sha256_hash for uri, meta in registry.verify_many(result.output_artifacts).items()}
        if entry.get("output_hashes") != hashes or commit.get("output_hashes") != hashes or audit.artifact_hashes != hashes:
            raise ValueError(f"Original complete output hash coverage changed: {task.task_id}")
        # Enforce the original persisted auditor identity; this is read-only
        # admission access, never a fresh scientific audit.
        registry.get_audited_many(signature, task, expected_auditor_version=audit.auditor_version,
                                  expected_auditor_fingerprint=audit.auditor_fingerprint)
        report = ValidationReport.model_validate(audit.report)
        if not report.overall_passed or report.stop_rule_triggered:
            raise ValueError("Original report does not admit scientific evidence")
        if task.capability == "donor_sensitivity":
            if not _completed_sensitivity(result.metrics):
                raise ValueError("Original donor sensitivity fit coverage is incomplete")
        if task.capability == "dataset_audit":
            adapted, decisions = ComputationalDAGPlanner.adapt_after_audit(deepcopy(tasks), state.preview_with_metrics(result.metrics))
            if result.metrics.get("planning_decisions", []) != decisions:
                raise ValueError("Saved planning decisions differ from the original observed contract policy")
            tasks[:] = adapted
        state_metrics = deepcopy(result.metrics)
        state_metrics.pop("planning_decisions", None)
        state.record_metrics(state_metrics, branch=branch)
        state.mark_completed(task.task_id, result.status)
        results.append(result)
        reports.append(report)
        records.append((task, result, report))
        index += 1
    if {task.task_id for task in tasks} != set(entries) or len(tasks) != len(entries):
        raise ValueError("Completed journal task coverage differs from the complete declared plan")
    return records, results, reports


def finalize_completed_journal(*, run_dir, compute_image, source_proof, source_proof_sha256):
    """Recover delivery after presentation failure, without scientific execution."""
    from eacbp.application.study_service import _run_lock, _autonomous_request_summary, _jsonable
    from eacbp.evidence.snapshot import write_study_snapshot, load_study_snapshot, render_snapshot_report

    run_dir, proof_path = Path(run_dir).resolve(strict=True), Path(source_proof).resolve(strict=True)
    with _run_lock(run_dir):
        envelope = _object(run_dir / "run_config.json")
        if envelope.get("schema_version") != 1 or envelope.get("import_completed") is not True:
            raise ValueError("Finalization requires fully imported, persisted native inputs")
        manifest = StudyManifest.model_validate(envelope.get("manifest", {}))
        config = RunConfig.from_mapping(envelope.get("config", {})).to_mapping()
        if not manifest.analysis_policy.strict_reproducibility:
            raise ValueError("Completed-journal finalization requires an original strict environment identity")
        if (_object(run_dir / "manifest.json") != envelope["manifest"]
                or RunConfig.from_mapping(_object(run_dir / "config.json")).to_mapping() != config
                or _object(run_dir / "source.json") != envelope.get("source")):
            raise ValueError("Mirrored saved manifest/config/source differs from its original envelope")
        failure_bytes = (run_dir / "summary.json").read_bytes()
        failure = json.loads(failure_bytes)
        if failure.get("status") != "failed":
            raise ValueError("This entry point requires a preserved presentation failure summary")
        run_id = UUID(envelope["last_run_id"]).hex
        historical_path = run_dir / "summaries" / f"{run_id}.json"
        historical = _object(historical_path)
        if historical.get("status") != "success" or historical.get("tasks_blocked") != 0 or historical.get("run_id") != run_id:
            raise ValueError("Original completed compute summary is missing or inconsistent")
        registry = ArtifactRegistry(str(run_dir / "artifacts"))
        journal_path = registry.storage.base_dir / "_runs" / f"{fingerprint(manifest.study_id)}.json"
        if not journal_path.is_file():
            raise ValueError("Original task journal is missing")
        with StudyJournal(registry.storage.base_dir, manifest.study_id) as journal:
            entries = journal.entries
            environment = _original_environment(entries, proof_path, source_proof_sha256)
            invocation_path = _original_invocation(run_dir, envelope, manifest, config, compute_image)
            events_path = _event_proof(run_dir, journal_path, run_id, entries, registry)
            source = envelope.get("source", {})
            digest = source.get("sha256")
            raw = registry.get_metadata(manifest.data.raw_artifact_uri)
            if (not digest or any(source.get(key) != digest for key in ("sha256_before_import", "sha256_after_import"))
                    or raw.created_by_task != "task_000_ingest" or raw.parameters.get("source_sha256") != digest
                    or config.get("autonomous_analysis", {}).get("source_sha256", digest) != digest):
                raise ValueError("Original imported source/count provenance is inconsistent")
            protected = [run_dir / name for name in ("run_config.json", "manifest.json", "config.json", "source.json",
                         "autonomy.json", "input_preparation.json", "error.json", "artifacts/.artifact_registry.json")]
            protected += [journal_path, historical_path, invocation_path, events_path, proof_path]
            before = {str(path): _hash(path) for path in protected if path.is_file()}
            capabilities = create_default_capability_registry()
            records, results, reports = _validate_completed_plan(manifest, config, entries, registry, capabilities, environment)
            if historical.get("tasks_executed") != len(results):
                raise ValueError("Original compute summary task coverage differs from the declared journal")
            context = SimpleNamespace(artifact_registry=registry, capability_registry=capabilities,
                                      task_history=results, audit_reports=reports, evidence_graph=EvidenceGraph())
            context.claim_engine = ClaimEngine(context.evidence_graph)
            for task, result, report in records:
                for node in ScientificOrchestrator.extract_evidence_from_result(context, task, result, report):
                    context.evidence_graph.add_evidence(node)
            ScientificOrchestrator._synthesize_study_claims(context, manifest)
            summary = deepcopy(historical)
            summary.update(evidence_nodes_count=len(context.evidence_graph.evidence_nodes),
                           claims_count=len(context.evidence_graph.claim_nodes),
                           claims=[claim.model_dump() for claim in context.evidence_graph.claim_nodes.values()])
            summary = _autonomous_request_summary(summary, manifest=manifest, config=config, orchestrator=context)
            # A fully computed and admitted plan can leave unsupported parts of
            # its scientific question unanswered. Publish that partial scope;
            # journal/plan admission above still rejects missing computations.
            finalization_id = uuid4().hex
            stage = run_dir / "finalizations" / finalization_id
            stage.mkdir(parents=True, exist_ok=False)
            (stage / "previous-summary.json").write_bytes(failure_bytes)
            (stage / "previous-run-config.json").write_bytes((run_dir / "run_config.json").read_bytes())
            if (run_dir / "error.json").is_file():
                (stage / "previous-error.json").write_bytes((run_dir / "error.json").read_bytes())
            receipt = {
                "schema_version": 1, "finalization_id": finalization_id, "kind": "presentation_only_completed_journal",
                "created_at_utc": datetime.now(timezone.utc).isoformat(), "original_compute_run_id": run_id,
                "original_compute_image": compute_image, "original_compute_environment": environment,
                "current_finalizer_environment": execution_environment(), "source_proof_sha256": source_proof_sha256,
                "source_proof_path": str(proof_path), "journal_sha256": before[str(journal_path)],
                "compute_event_log_sha256": before[str(events_path)], "execution_receipt_sha256": before[str(invocation_path)],
                "saved_config_fingerprint": stable_fingerprint(config), "saved_config_file_sha256": before[str(run_dir / "config.json")],
                "config_binding": "Mirrored saved configuration and reconstructed declared DAG/contracts; the original journal did not store a separate full-config signature",
                "original_failure_summary_sha256": hashlib.sha256(failure_bytes).hexdigest(),
                "protected_native_file_hashes": before, "task_signatures": {task: entry["signature"] for task, entry in entries.items()},
                "scientific_execution_calls": 0, "scientific_audit_calls": 0,
                "analysis_request_status": summary.get("analysis_request_status"),
                "requirement_contract_id": summary.get("requirement_contract_id"),
                "attestation_scope": "Local original-image/source/journal consistency; not a signed remote execution proof",
            }
            _write_json(stage / "receipt.json", receipt)
            summary["finalization"] = {"kind": receipt["kind"], "finalization_id": finalization_id,
                                       "receipt": f"finalizations/{finalization_id}/receipt.json", "receipt_sha256": _hash(stage / "receipt.json")}
            summary["run_dir"], summary["report"], summary["snapshot"] = str(run_dir), str(run_dir / "report.md"), str(run_dir / "snapshot.json")
            write_study_snapshot(stage / "snapshot.json", manifest=manifest, config=config, summary=_jsonable(summary),
                                 evidence_graph=context.evidence_graph, artifact_registry=registry,
                                 task_history=results, audit_reports=reports)
            load_study_snapshot(stage / "snapshot.json")
            markdown = render_snapshot_report(stage / "snapshot.json", registry)
            (stage / "report.md").write_text(markdown, encoding="utf-8")
            snapshot_payload = _object(stage / "snapshot.json")
            summary = snapshot_payload["summary"]
            _write_json(stage / "summary.json", summary)
            if before != {str(path): _hash(path) for path in protected if path.is_file()}:
                raise ValueError("Original compute files changed during presentation finalization")
            # The successful summary is the final publication marker. A crash
            # before it leaves the old failure and the complete staged files.
            _write_json(run_dir / "snapshot.json", snapshot_payload)
            _write_report(run_dir / "report.md", markdown)
            _write_json(run_dir / "summary.json", summary)
            if before != {str(path): _hash(path) for path in protected if path.is_file()}:
                raise ValueError("Original compute files changed during publication")
            return summary


__all__ = ["finalize_completed_journal"]
