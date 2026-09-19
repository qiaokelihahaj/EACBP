"""
Scientific Orchestration Engine: Full execution loop tying dynamic DAG execution,
independent auditing, structured multi-plane evidence extraction, and 5-pillar claim synthesis.
"""

from typing import Dict, Any, List, Optional
from copy import deepcopy
import time

from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
    ClaimNode,
    ClaimType,
    LanguageTier,
)
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities import CapabilityRegistry, create_default_capability_registry
from eacbp.auditor import ScientificAuditor, ValidationReport
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.claim import ClaimEngine
from eacbp.orchestrator.router import CapabilityRouter
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.schemas.runtime import ExecutionState
from eacbp.orchestrator.execution import TaskExecutor
from eacbp.orchestrator.resume import ResumeManager
from eacbp.orchestrator.admission import EvidenceAdmission
from eacbp.orchestrator.events import RunEventJournal
from eacbp.orchestrator.planning import build_study_tasks, resolve_task_contract

class ScientificOrchestrator:
    """The central scientific orchestrator coordinating computation, independent validation, and evidence synthesis."""

    def __init__(
        self,
        artifact_registry: Optional[ArtifactRegistry] = None,
        capability_registry: Optional[CapabilityRegistry] = None,
        auditor: Optional[ScientificAuditor] = None,
    ):
        self.artifact_registry = artifact_registry if artifact_registry is not None else ArtifactRegistry()
        # ``None`` selects the complete built-in assembly.  An explicitly
        # supplied registry is a dependency injection boundary: preserve it
        # exactly so callers can provide replacements, test doubles, or a
        # deliberately restricted capability surface.
        self.capability_registry = (
            create_default_capability_registry()
            if capability_registry is None else capability_registry
        )

        descriptor_validators = self.capability_registry.audit_validators()
        self.auditor = auditor or ScientificAuditor(additional_validators=descriptor_validators)
        if auditor is not None and descriptor_validators:
            if not hasattr(auditor, "additional_validators"):
                raise ValueError("Custom auditor must support additional_validators for declared capability audits")
            auditor.additional_validators.extend(descriptor_validators)
        self.router = CapabilityRouter(self.capability_registry)
        self.evidence_graph = EvidenceGraph()
        self.claim_engine = ClaimEngine(self.evidence_graph)
        
        self.task_history: List[TaskResult] = []
        self.audit_reports: List[ValidationReport] = []
        self.current_state: Dict[str, Any] = {}
        self.execution_state = ExecutionState()
        self.task_executor = TaskExecutor(self.artifact_registry, self.capability_registry)
        self.resume_manager = ResumeManager(self.artifact_registry)
        # Resolve the method at call time so tests and plugin users can still
        # replace ``extract_evidence_from_result`` on an orchestrator instance.
        self.evidence_admission = EvidenceAdmission(
            self.auditor,
            lambda contract, result, report: self.extract_evidence_from_result(contract, result, report),
        )

    def extract_evidence_from_result(self, contract, result, report):
        from eacbp.evidence.extraction import extract_evidence
        descriptor = (self.capability_registry.describe(contract.capability, result.method_used or contract.method)
                      if self.capability_registry.has(contract.capability, result.method_used or contract.method) else None)
        nodes = (self.capability_registry.extract_evidence(contract, result, report, self.artifact_registry)
                 if descriptor is not None and descriptor.evidence_extractor is not None
                 else extract_evidence(contract, result, report, self.artifact_registry))
        if any(node.evidence_id in self.evidence_graph.evidence_nodes for node in nodes):
            raise ValueError("Evidence IDs must be unique across study tasks")
        for node in nodes:
            for key in ("target_cell_type", "target_branch", "fdr_family"):
                if key in contract.parameters:
                    node.biological_context[key] = contract.parameters[key]
        return nodes

    def run_study(self, manifest: StudyManifest, current_state=None):
        """Run a study with a separate diagnostic event log for this invocation."""
        self.event_journal = RunEventJournal(self.artifact_registry.storage.base_dir, manifest.study_id)
        self.task_executor.event_journal = self.event_journal
        self.event_journal.emit("run_started")
        try:
            summary = self._run_study(manifest, current_state)
        except BaseException as exc:
            self.event_journal.emit("run_failed", error_type=type(exc).__name__, error=str(exc))
            raise
        self.event_journal.emit("run_finished", status=summary["status"],
                                tasks_executed=summary["tasks_executed"],
                                tasks_blocked=summary["tasks_blocked"])
        summary["run_id"] = self.event_journal.run_id
        summary["event_log"] = str(self.event_journal.path.resolve())
        summary["observability_warnings"] = list(self.event_journal.warnings)
        return summary

    def _run_study(self, manifest: StudyManifest, current_state=None):
        """Run dependency-ordered tasks; only audited outputs enter evidence synthesis."""
        from eacbp.orchestrator.checkpoint import StudyJournal, execution_environment
        from eacbp.artifact.transaction import TaskArtifactTransaction
        self.execution_state = ExecutionState.from_input(current_state)
        self.current_state = self.execution_state.current_state
        self.task_history = []
        self.audit_reports = []
        self.evidence_graph = EvidenceGraph()
        self.claim_engine = ClaimEngine(self.evidence_graph)
        self.manifest = manifest
        if "method_profile" not in self.execution_state.config:
            raw_uri = manifest.data.raw_artifact_uri
            simulated = self.execution_state.current_state.get("mode") == "demo"
            if raw_uri and self.artifact_registry.exists(raw_uri):
                simulated = simulated or bool(self.artifact_registry.get_metadata(raw_uri).summary_metrics.get("is_simulated"))
            self.execution_state.set_derived("method_profile", "baseline" if simulated else "standard")
            self.current_state = self.execution_state.current_state
        resume = self.execution_state.resume_requested
        planned_tasks = build_study_tasks(manifest, self.execution_state.current_state, self.capability_registry)
        self.event_journal.emit("plan_created", task_count=len(planned_tasks), resume=resume)
        completed = {}
        environment = execution_environment() if manifest.analysis_policy.strict_reproducibility else {}
        planning_decisions = []
        with StudyJournal(self.artifact_registry.storage.base_dir, manifest.study_id) as journal:
            if journal.entries and not resume:
                raise ValueError("Study already has a task journal. Use resume=True or a new study/run directory.")
            for task in planned_tasks:
                task_started = time.monotonic()
                task.parameters["study_id"] = manifest.study_id
                blocked = [dep for dep in task.depends_on if completed.get(dep) != TaskStatus.SUCCESS]
                if blocked:
                    result = TaskResult(task_id=task.task_id, capability=task.capability,
                        method_used=task.method or "unresolved", status=TaskStatus.BLOCKED,
                        error_type=ExecutionFailureType.INSUFFICIENT_EVIDENCE,
                        error_message=f"Dependencies did not pass audit: {blocked}")
                    self.task_history.append(result)
                    completed[task.task_id] = result.status
                    self.execution_state.mark_completed(task.task_id, result.status)
                    self.event_journal.emit("task_blocked", task_id=task.task_id, dependencies=blocked)
                    continue
                self.event_journal.emit("task_started", task_id=task.task_id, capability=task.capability)
                try:
                    pending_nodes = []
                    admission = None
                    branch = task.parameters.get("target_branch")
                    resolve_task_contract(task, manifest, self.execution_state.state_for_branch(branch),
                                          self.capability_registry, self.artifact_registry, router=self.router)
                    self.event_journal.emit("method_resolved", task_id=task.task_id, method=task.method)
                    input_hashes = self.resume_manager.input_hashes(task)
                    signature = self.resume_manager.signature(task, input_hashes, manifest, environment)
                    lookup = self.resume_manager.lookup(journal, task, signature)
                    if lookup.reused:
                        result = lookup.result
                        self.event_journal.emit("task_resumed", task_id=task.task_id, signature=signature)
                    else:
                        if any(self.artifact_registry.exists(uri) for uri in task.expected_outputs):
                            raise ValueError("Uncommitted or rejected outputs exist; use a new run directory.")
                        execution = self.task_executor.execute(task)
                        result = execution.result
                        for key in ("target_cell_type", "target_branch", "fdr_family"):
                            if key in task.parameters:
                                result.metrics[key] = task.parameters[key]
                        staged = execution.staged
                        if result.status == TaskStatus.SUCCESS:
                            try:
                                self.artifact_registry.commit_task(staged, signature, result)
                            finally:
                                staged.close()
                            self.event_journal.emit("artifacts_committed", task_id=task.task_id,
                                                    outputs=result.output_artifacts, signature=signature)
                        # Preserve completed computation before auditing. A crash
                        # during audit can resume by re-auditing these exact outputs.
                        journal.entries[task.task_id] = self.resume_manager.journal_entry(
                            signature=signature, result=result, phase="computed", environment=environment
                        )
                        journal.save()
                    if result.status == TaskStatus.SUCCESS:
                        # Invalidate any previous scientific pass before the
                        # auditor runs.  A crash during audit must leave a
                        # durable pending receipt instead of a stale pass.
                        audited_result = result.model_copy(deep=True)
                        self.artifact_registry.begin_audit(
                            signature, task, audited_result, auditor=self.auditor
                        )
                        audit_started = time.monotonic()
                        self.event_journal.emit("audit_started", task_id=task.task_id, signature=signature)
                        report = self.auditor.audit_task(task, result, self.artifact_registry)
                        admission = self.evidence_admission.stage_report(
                            task, result, manifest, self.execution_state, planned_tasks, report
                        )
                        pending_nodes = admission.pending_nodes
                    journal.entries[task.task_id] = self.resume_manager.journal_entry(
                        signature=signature, result=result, phase="audited", environment=environment
                    )
                    journal.save()
                    audit_admitted = True
                    if admission is not None:
                        audit_record = self.artifact_registry.record_audit(
                            signature,
                            task,
                            audited_result,
                            admission.audit,
                            auditor=self.auditor,
                        )
                        audit_admitted = getattr(audit_record.status, "value", audit_record.status) == "passed"
                        self.event_journal.emit("audit_finished" if audit_admitted else "audit_rejected",
                                                task_id=task.task_id, signature=signature,
                                                duration_seconds=max(0.0, time.monotonic() - audit_started),
                                                reason=audit_record.rejection_reason,
                                                checks_failed=[c.check_name for c in admission.audit.checks if not c.passed])
                        if not audit_admitted:
                            result.status = TaskStatus.SCIENTIFIC_FAILURE
                            result.error_type = ExecutionFailureType.INSUFFICIENT_EVIDENCE
                            result.error_message = (
                                getattr(audit_record, "rejection_reason", None)
                                or "Durable audit receipt rejected output."
                            )
                            pending_nodes = []
                            # Keep the journal result consistent with the
                            # final durable audit decision.  This second write
                            # is only reached after the normal audited entry is
                            # durable and is therefore safe for resume.
                            journal.entries[task.task_id] = self.resume_manager.journal_entry(
                                signature=signature, result=result, phase="audited", environment=environment
                            )
                            journal.save()
                    # Only now are audit history, planning observations, and
                    # evidence admitted.  A failed durable journal write leaves
                    # all three untouched and the next run re-audits the saved
                    # computed outputs.
                    if admission is not None:
                        self.audit_reports.append(admission.audit)
                        if admission.accepted and audit_admitted:
                            state_metrics = dict(result.metrics)
                            state_metrics.pop("planning_decisions", None)
                            self.execution_state.record_metrics(state_metrics, branch=branch)
                            if admission.adapted_tasks is not None:
                                planned_tasks[:] = admission.adapted_tasks
                                planning_decisions.extend(admission.planning_decisions)
                                if admission.planning_decisions:
                                    self.event_journal.emit("plan_adapted", task_id=task.task_id,
                                                            decisions=admission.planning_decisions)
                            self.execution_state.record_decisions(admission.planning_decisions)
                            self.current_state = self.execution_state.current_state
                    for node in pending_nodes:
                        self.evidence_graph.add_evidence(node)
                except Exception as exc:
                    result = TaskResult(task_id=task.task_id, capability=task.capability,
                        method_used=task.method or "unresolved", status=TaskStatus.EXECUTION_FAILURE,
                        error_type=ExecutionFailureType.CODE_ERROR, error_message=f"{type(exc).__name__}: {exc}")
                self.task_history.append(result)
                completed[task.task_id] = result.status
                self.execution_state.mark_completed(task.task_id, result.status)
                self.event_journal.emit("task_finished" if result.status == TaskStatus.SUCCESS else "task_failed",
                                        task_id=task.task_id, status=result.status.value,
                                        method=result.method_used, error=result.error_message,
                                        error_type=result.error_type.value if result.error_type else None,
                                        duration_seconds=max(0.0, time.monotonic() - task_started))
        self._synthesize_study_claims(manifest)
        failures = [r for r in self.task_history if r.status not in (TaskStatus.SUCCESS, TaskStatus.BLOCKED)]
        return {
            "study_id": manifest.study_id,
            "status": "failed" if failures else "success",
            "tasks_executed": sum(r.status != TaskStatus.BLOCKED for r in self.task_history),
            "tasks_blocked": sum(r.status == TaskStatus.BLOCKED for r in self.task_history),
            "artifacts_created": len(self.artifact_registry.list_artifacts(study_id=manifest.study_id)),
            "evidence_nodes_count": len(self.evidence_graph.evidence_nodes),
            "claims_count": len(self.evidence_graph.claim_nodes),
            "claims": [c.model_dump() for c in self.evidence_graph.claim_nodes.values()],
            "planning_decisions": planning_decisions,
            "failures": [{"task_id": r.task_id, "error": r.error_message} for r in failures],
            "is_simulated": any(e.is_simulated for e in self.evidence_graph.evidence_nodes.values()),
        }

    def _synthesize_study_claims(self, manifest):
        """Each sentence restates its own evidence; no disease-specific conclusions."""
        for index, node in enumerate(self.evidence_graph.evidence_nodes.values()):
            if node.polarity != EvidencePolarity.SUPPORTING:
                continue
            tier = LanguageTier.LEVEL_1_OBSERVATION
            causal = "observational"
            statement = node.summary
            if node.type == EvidenceType.PERTURBATION:
                tier = LanguageTier.LEVEL_4_HYPOTHESIS
                causal = "in_silico_perturbed"
            elif node.metrics.get("source_mode") == "local_curated_unverified":
                tier = LanguageTier.LEVEL_4_HYPOTHESIS
                statement = f"Unverified local reference context was returned for this query (item {node.evidence_id}); biological interpretation requires source verification."
                if manifest.analysis_policy.prior_guided_analysis or manifest.hypotheses.user_provided:
                    statement = "[PRIOR-GUIDED HYPOTHESIS TESTING] " + statement
            elif self.claim_engine.has_valid_statistics(node):
                tier = LanguageTier.LEVEL_2_STATISTICAL_INFERENCE
            self.claim_engine.create_claim(
                claim_id=f"C_{node.evidence_id}", statement=statement,
                language_tier=tier, claim_type=ClaimType.DESCRIPTIVE,
                causal_status=causal, support_evidence_ids=[node.evidence_id],
            )
