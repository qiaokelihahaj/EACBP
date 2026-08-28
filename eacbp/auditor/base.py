"""
Scientific Auditor base interfaces and validation report definitions.
Core principle: Executor cannot certify itself (Author != Reviewer).
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from enum import Enum
from pydantic import BaseModel, Field
from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.artifact.registry import ArtifactRegistry


class ValidationSeverity(str, Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    STOP_RULE = "stop_rule"


class ValidationCheck(BaseModel):
    check_name: str
    passed: bool
    severity: ValidationSeverity = Field(ValidationSeverity.INFO)
    message: str
    metrics: Dict[str, Any] = Field(default_factory=dict)


class ValidationReport(BaseModel):
    auditor_name: str
    target_task_id: str
    target_artifact_uri: Optional[str] = None
    checks: List[ValidationCheck] = Field(default_factory=list)
    overall_passed: bool = Field(True)
    stop_rule_triggered: bool = Field(False)
    recommendations: List[str] = Field(default_factory=list)

    def add_check(
        self,
        name: str,
        passed: bool,
        severity: ValidationSeverity,
        message: str,
        metrics: Optional[Dict[str, Any]] = None
    ):
        metrics = metrics or {}
        check = ValidationCheck(
            check_name=name,
            passed=passed,
            severity=severity,
            message=message,
            metrics=metrics,
        )
        self.checks.append(check)
        if not passed and severity in (ValidationSeverity.ERROR, ValidationSeverity.STOP_RULE):
            self.overall_passed = False
        if not passed and severity == ValidationSeverity.STOP_RULE:
            self.stop_rule_triggered = True


class BaseAuditor(ABC):
    """Independent auditor evaluating computational, statistical, or biological validity."""

    def __init__(self, auditor_name: str):
        self.auditor_name = auditor_name

    @abstractmethod
    def audit(
        self,
        contract: TaskContract,
        result: TaskResult,
        registry: ArtifactRegistry
    ) -> ValidationReport:
        pass
