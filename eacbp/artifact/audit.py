"""Public artifact-audit schema and identity helpers.

The durable records live alongside the schema definitions so the registry can
load old/new indexes without importing auditor implementations.  This module
provides the artifact-plane import path for integrations that naturally look
for audit records under ``eacbp.artifact``.
"""

from eacbp.schemas.audit import (
    ArtifactAuditRecord,
    AuditRecord,
    AuditRecordStatus,
    AuditStatus,
    auditor_version,
    fingerprint_auditor,
    model_to_jsonable,
    stable_fingerprint,
)

__all__ = [
    "ArtifactAuditRecord",
    "AuditRecord",
    "AuditRecordStatus",
    "AuditStatus",
    "auditor_version",
    "fingerprint_auditor",
    "model_to_jsonable",
    "stable_fingerprint",
]
