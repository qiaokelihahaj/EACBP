"""Independent audit and evidence-admission boundary."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Any, Callable, List, Optional

from eacbp.auditor import ScientificAuditor, ValidationReport
from eacbp.schemas.evidence import EvidenceNode
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import ExecutionFailureType, TaskContract, TaskResult, TaskStatus
from eacbp.schemas.runtime import ExecutionState
from eacbp.orchestrator.dag import ComputationalDAGPlanner


@dataclass
class AdmissionOutcome:
    """Audit result and staged planning/evidence effects.

    The outcome is intentionally side-effect free with respect to the
    orchestrator's evidence graph, audit history, and execution state.  The
    caller can persist its checkpoint first, then apply these effects.
    """

    audit: ValidationReport
    pending_nodes: List[EvidenceNode] = field(default_factory=list)
    prospective_state: Optional[dict] = None
    adapted_tasks: Optional[list] = None
    planning_decisions: List[dict] = field(default_factory=list)

    @property
    def accepted(self) -> bool:
        return self.audit.overall_passed and not self.audit.stop_rule_triggered


class EvidenceAdmission:
    """Audit a successful computation and stage admissible evidence."""

    def __init__(
        self,
        auditor: ScientificAuditor,
        evidence_extractor: Callable[[TaskContract, TaskResult, ValidationReport], List[EvidenceNode]],
    ):
        self.auditor = auditor
        self.evidence_extractor = evidence_extractor

    def stage_report(
        self,
        task: TaskContract,
        result: TaskResult,
        manifest: StudyManifest,
        state: ExecutionState,
        planned_tasks: list,
        audit: ValidationReport,
    ) -> AdmissionOutcome:
        """Stage a report that has already been produced by the auditor.

        The durable registry integration needs to mark a receipt pending before
        invoking the auditor.  This entry point lets the loop keep that single
        auditor call while reusing all admission rules here.
        """

        outcome = self._stage(task, result, manifest, state, planned_tasks, audit)
        if outcome.accepted:
            outcome.pending_nodes = self.evidence_extractor(task, result, audit)
        return outcome

    @staticmethod
    def _stage(
        task: TaskContract,
        result: TaskResult,
        manifest: StudyManifest,
        state: ExecutionState,
        planned_tasks: list,
        audit: ValidationReport,
    ) -> AdmissionOutcome:
        if not audit.overall_passed or audit.stop_rule_triggered:
            result.status = TaskStatus.SCIENTIFIC_FAILURE
            result.error_type = ExecutionFailureType.INSUFFICIENT_EVIDENCE
            result.error_message = (
                "; ".join(c.message for c in audit.checks if not c.passed)
                or "Scientific audit rejected output."
            )
            return AdmissionOutcome(audit=audit)

        # Metrics are staged in the result and the journal first.  The actual
        # ExecutionState update is performed by the loop only after its durable
        # audited checkpoint succeeds.
        prospective_state = state.preview_with_metrics(result.metrics)
        adapted_tasks = None
        decisions: List[dict] = []
        if task.capability == "dataset_audit":
            adapted_tasks, decisions = ComputationalDAGPlanner.adapt_after_audit(
                deepcopy(planned_tasks), prospective_state
            )
            result.metrics["planning_decisions"] = decisions

        return AdmissionOutcome(
            audit=audit,
            prospective_state=prospective_state,
            adapted_tasks=adapted_tasks,
            planning_decisions=decisions,
            pending_nodes=[],
        )
