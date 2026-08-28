"""
Scientific Policies and Stop Rules for workflow orchestration.
"""

from typing import Dict, Any, Tuple, Optional
from eacbp.schemas.study import StudyManifest
from eacbp.auditor.base import ValidationReport


class ScientificPolicy:
    """Evaluates stop rules and scientific gating policies across execution stages."""

    @staticmethod
    def evaluate_deg_policy(manifest: StudyManifest, audit_metrics: Dict[str, Any]) -> Dict[str, Any]:
        """Enforces policy on differential expression based on donor replicate count."""
        min_reps = audit_metrics.get("min_replicates", 1)
        if min_reps < 2:
            return {
                "confirmatory_allowed": False,
                "preferred_method": "cell_level_mannwhitney",
                "policy_notice": "STOP RULE ACTIVE: Fewer than 2 biological replicates per condition. Confirmatory DEG disallowed; flagged as exploratory.",
            }
        elif min_reps >= 3 and manifest.analysis_policy.prefer_pseudobulk:
            return {
                "confirmatory_allowed": True,
                "preferred_method": "deg_pseudobulk_v1",
                "policy_notice": "Pseudobulk donor aggregation enabled for statistical rigor.",
            }
        else:
            return {
                "confirmatory_allowed": True,
                "preferred_method": "deg_pseudobulk_v1",
                "policy_notice": "Standard DEG enabled.",
            }

    @staticmethod
    def evaluate_trajectory_stop_rules(validation_report: ValidationReport) -> Tuple[bool, Optional[str]]:
        """Checks whether trajectory stability meets scientific thresholds."""
        for check in validation_report.checks:
            if check.check_name == "trajectory_subsampling_stability" and not check.passed:
                return False, "STOP RULE: Trajectory failed stability audit (< 0.60 correlation). Mechanistic state claims suppressed."
        return True, None
