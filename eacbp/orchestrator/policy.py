"""
Scientific Policies and Stop Rules for workflow orchestration.
"""

from typing import Dict, Any, Tuple, Optional
import warnings
from eacbp.schemas.study import StudyManifest
from eacbp.auditor.base import ValidationReport


class ScientificPolicy:
    """Legacy policy facade retained for import compatibility.

    Study execution now applies policy through the router and independent
    auditors.  This class is kept so older integrations fail softly while
    making its deprecated status visible to callers that still invoke it.
    """

    def __init__(self) -> None:
        warnings.warn(
            "ScientificPolicy is deprecated; use the orchestrator router and independent auditors.",
            DeprecationWarning,
            stacklevel=2,
        )

    @staticmethod
    def _warn_deprecated() -> None:
        warnings.warn(
            "ScientificPolicy methods are deprecated; use the orchestrator router and independent auditors.",
            DeprecationWarning,
            stacklevel=3,
        )

    @staticmethod
    def evaluate_deg_policy(manifest: StudyManifest, audit_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """Enforces policy on differential expression based on donor replicate count."""
        ScientificPolicy._warn_deprecated()
        min_reps = audit_metrics.get("min_replicates", 1)
        if min_reps < max(2, manifest.constraints.min_biological_replicates):
            return {
                "confirmatory_allowed": False,
                "preferred_method": "single_cell_exploratory_welch_v1",
                "policy_notice": "STOP RULE ACTIVE: Fewer than 2 biological replicates per condition. Confirmatory DEG disallowed; flagged as exploratory.",
            }
        elif min_reps >= 3 and manifest.analysis_policy.prefer_pseudobulk:
            return {
                "confirmatory_allowed": True,
                "preferred_method": "donor_pseudobulk_welch_v1",
                "policy_notice": "Pseudobulk donor aggregation enabled for statistical rigor.",
            }
        else:
            return {
                "confirmatory_allowed": True,
                "preferred_method": "donor_pseudobulk_welch_v1",
                "policy_notice": "Standard DEG enabled.",
            }

    @staticmethod
    def evaluate_trajectory_stop_rules(validation_report: ValidationReport) -> Tuple[bool, Optional[str]]:
        """Checks whether trajectory stability meets scientific thresholds."""
        ScientificPolicy._warn_deprecated()
        for check in validation_report.checks:
            if check.check_name == "trajectory_subsampling_stability" and not check.passed:
                return False, "STOP RULE: Trajectory failed stability audit (< 0.60 correlation). Mechanistic state claims suppressed."
        return True, None
