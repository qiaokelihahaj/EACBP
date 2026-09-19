"""Adapt a capability's independent validator into the audit suite."""
from eacbp.auditor.base import ValidationReport


class DescriptorValidator:
    def __init__(self, descriptor):
        self.capability_name = descriptor.capability_name
        self.method = descriptor.method
        self.validator = descriptor.validator

    def audit(self, contract, result, registry):
        if (contract.capability, result.method_used) != (self.capability_name, self.method):
            return ValidationReport(auditor_name="descriptor_not_applicable", target_task_id=contract.task_id)
        callback = getattr(self.validator, "audit", self.validator)
        report = callback(contract, result, registry)
        if not isinstance(report, ValidationReport) or report.target_task_id != contract.task_id:
            raise ValueError("Descriptor validator must return a ValidationReport for the current task")
        return report
