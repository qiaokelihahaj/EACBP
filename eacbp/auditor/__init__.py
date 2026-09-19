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

    def __init__(self, additional_validators=None):
        self.computational_validator = ComputationalValidator()
        self.statistical_validator = StatisticalValidator()
        self.biological_validator = BiologicalValidator()
        from eacbp.auditor.advanced_statistics import AdvancedStatisticsValidator
        from eacbp.auditor.advanced_qc import AdvancedQCValidator
        from eacbp.auditor.advanced_communication import LianaValidator
        self.additional_validators = [AdvancedStatisticsValidator(), AdvancedQCValidator(), LianaValidator(), *list(additional_validators or [])]

    def audit_task(self, contract, result, registry):
        """Runs computational, statistical, and biological audits for a completed task."""
        comp_reports = [self.computational_validator.audit(
            contract, result.model_copy(update={"output_artifacts": [uri]}), registry)
            for uri in result.output_artifacts]
        if not comp_reports:
            comp_reports = [self.computational_validator.audit(contract, result, registry)]
        stat_report = self.statistical_validator.audit(contract, result, registry)
        bio_report = self.biological_validator.audit(contract, result, registry)

        all_checks = [c for r in comp_reports for c in r.checks] + stat_report.checks + bio_report.checks
        extra_reports = [validator.audit(contract, result, registry) for validator in self.additional_validators]
        all_checks.extend(c for report in extra_reports for c in report.checks)
        from eacbp.auditor.requirements import run_payload_requirements, missing_required_checks
        all_checks.extend(run_payload_requirements(contract, result, registry, all_checks))
        all_checks.extend(missing_required_checks(contract.validation_requirements, all_checks))
        overall_passed = all(r.overall_passed for r in comp_reports) and stat_report.overall_passed and bio_report.overall_passed
        overall_passed = overall_passed and all(r.overall_passed for r in extra_reports)
        overall_passed = overall_passed and not any(not c.passed and c.severity in (ValidationSeverity.ERROR, ValidationSeverity.STOP_RULE) for c in all_checks)
        stop_rule_triggered = any(r.stop_rule_triggered for r in comp_reports) or stat_report.stop_rule_triggered or bio_report.stop_rule_triggered
        stop_rule_triggered = stop_rule_triggered or any(r.stop_rule_triggered for r in extra_reports)
        recommendations = [c for r in comp_reports for c in r.recommendations] + stat_report.recommendations + bio_report.recommendations
        recommendations.extend(item for report in extra_reports for item in report.recommendations)

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
