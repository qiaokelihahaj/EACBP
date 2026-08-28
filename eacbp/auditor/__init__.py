"""
EACBP Scientific Auditor Plane: Computational, Statistical, and Biological validators.
"""

from eacbp.auditor.base import (
    BaseAuditor,
    ValidationReport,
    ValidationCheck,
    ValidationSeverity,
)
from eacbp.auditor.computational import ComputationalValidator
from eacbp.auditor.statistical import StatisticalValidator
from eacbp.auditor.biological import BiologicalValidator, CANONICAL_MARKERS


class ScientificAuditor:
    """Unified Independent Scientific Auditor orchestrating all validation checks."""

    def __init__(self):
        self.computational_validator = ComputationalValidator()
        self.statistical_validator = StatisticalValidator()
        self.biological_validator = BiologicalValidator()

    def audit_task(self, contract, result, registry):
        """Runs computational, statistical, and biological audits for a completed task."""
        comp_report = self.computational_validator.audit(contract, result, registry)
        stat_report = self.statistical_validator.audit(contract, result, registry)
        bio_report = self.biological_validator.audit(contract, result, registry)

        all_checks = comp_report.checks + stat_report.checks + bio_report.checks
        overall_passed = comp_report.overall_passed and stat_report.overall_passed and bio_report.overall_passed
        stop_rule_triggered = comp_report.stop_rule_triggered or stat_report.stop_rule_triggered or bio_report.stop_rule_triggered
        recommendations = comp_report.recommendations + stat_report.recommendations + bio_report.recommendations

        return ValidationReport(
            auditor_name="scientific_auditor_suite",
            target_task_id=contract.task_id,
            target_artifact_uri=result.output_artifacts[0] if result.output_artifacts else None,
            checks=all_checks,
            overall_passed=overall_passed,
            stop_rule_triggered=stop_rule_triggered,
            recommendations=recommendations,
        )


__all__ = [
    "BaseAuditor",
    "ValidationReport",
    "ValidationCheck",
    "ValidationSeverity",
    "ComputationalValidator",
    "StatisticalValidator",
    "BiologicalValidator",
    "CANONICAL_MARKERS",
    "ScientificAuditor",
]
