"""Read-only acceptance closure and reported scientific conclusions."""
from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field


class VerificationFinding(BaseModel):
    model_config = ConfigDict(extra="forbid")

    code: str
    status: Literal["passed", "failed"]
    message: str
    evidence_paths: list[str] = Field(default_factory=list)


class VerifiedFile(BaseModel):
    model_config = ConfigDict(extra="forbid")

    path: str
    sha256: str
    bytes: int
    verified_path: str
    source: Literal["current_file", "bound_archive"]


class ReportedEngineering(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["passed", "failed", "not_verified"] = "not_verified"
    reported_status: str | None = None
    native_scope_count: int = 0
    independent_pair_count: int = 0
    durable_task_count: int = 0
    validation_scope: str = "Bound independent checker receipts; native computations were not rerun."


class ReportedScience(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["passed", "failed", "not_estimable", "not_verified"] = "not_verified"
    reported_status: str | None = None
    reference_field: str | None = None
    hypothesis: dict[str, Any] = Field(default_factory=dict)
    improvement: dict[str, Any] = Field(default_factory=dict)
    absolute_acceptance: dict[str, Any] = Field(default_factory=dict)
    failure_conditions: list[dict[str, Any]] = Field(default_factory=list)
    branches: list[dict[str, Any]] = Field(default_factory=list)
    reference: dict[str, Any] = Field(default_factory=dict)
    interpretation: str = ""
    validation_scope: str = "Reported scientific checker result only; reference values were not parsed and scores were not recomputed."


class VerificationResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal[1] = 1
    verification_status: Literal["passed", "failed"] = "failed"
    index: dict[str, Any]
    engineering: ReportedEngineering = Field(default_factory=ReportedEngineering)
    science: ReportedScience = Field(default_factory=ReportedScience)
    findings: list[VerificationFinding] = Field(default_factory=list)
    verified_files: list[VerifiedFile] = Field(default_factory=list)
    evidence: dict[str, Any] = Field(default_factory=dict)
    relationships: dict[str, Any] = Field(default_factory=dict)
    prior_knowledge_and_exposure: dict[str, Any] = Field(default_factory=dict)
    publication_logs_outside_index: list[str] = Field(default_factory=list)
    labels_read: Literal[False] = False
    scores_recomputed: Literal[False] = False
    scientific_result_reestablished: Literal[False] = False
    integrity_scope: str = "Sequential file-byte and receipt-binding checks; no external attestation or atomic filesystem snapshot."
    summary: str = ""
