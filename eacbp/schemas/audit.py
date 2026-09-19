"""Durable artifact-audit records and audit identity helpers.

An artifact payload being present in :class:`ArtifactRegistry` means that the
computation was published.  It does not mean that the payload is suitable for
scientific evidence.  ``ArtifactAuditRecord`` is the durable boundary between
those two states.  Records are deliberately kept in a small schema module so
the registry can persist them without importing the auditor plane (which would
create a circular dependency).

The orchestration hook is intentionally two phase::

    registry.begin_audit(signature, contract, result, ...)
    report = auditor.audit_task(contract, result, registry)
    registry.record_audit(signature, contract, result, report, ...)

The first call invalidates an earlier pass.  If the process exits while the
auditor is running, the durable state therefore remains ``pending`` and an
audited query cannot admit the output on resume.
"""

from __future__ import annotations

import hashlib
import inspect
import json
from datetime import datetime, timezone
from enum import Enum
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

from pydantic import BaseModel, ConfigDict, Field


class AuditRecordStatus(str, Enum):
    """Durable lifecycle state for one task-signature audit."""

    PENDING = "pending"
    PASSED = "passed"
    REJECTED = "rejected"

    # A failed audit is represented by ``rejected`` in the durable index.  The
    # alias is useful to callers that use execution terminology while keeping a
    # single canonical on-disk value.
    FAILED = "rejected"


def model_to_jsonable(value: Any) -> Any:
    """Return a deterministic JSON-compatible representation of *value*.

    Pydantic models are intentionally handled without importing any project
    schema.  This keeps the helper usable for ``TaskContract``, ``TaskResult``
    and ``ValidationReport`` while preserving compatibility with plain mapping
    contexts supplied by integrations.
    """

    if hasattr(value, "model_dump"):
        try:
            # ``mode='json'`` eagerly serializes numpy arrays and other
            # scientific values before this helper can normalize them.  Keep
            # Python values intact and recurse below instead.
            return model_to_jsonable(value.model_dump(mode="python"))
        except TypeError:  # pragma: no cover - pydantic v1 compatibility
            return model_to_jsonable(value.model_dump())
    if isinstance(value, Mapping):
        return {str(key): model_to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [model_to_jsonable(item) for item in value]
    if isinstance(value, set):
        return sorted(model_to_jsonable(item) for item in value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, (datetime,)):
        return value.isoformat()
    # numpy scalars/arrays are optional dependencies of the package.  Avoid an
    # unconditional import here, but support them when present in a context.
    if hasattr(value, "item") and callable(value.item):
        try:
            return model_to_jsonable(value.item())
        except (TypeError, ValueError):
            pass
    if hasattr(value, "tolist") and callable(value.tolist):
        try:
            return model_to_jsonable(value.tolist())
        except (TypeError, ValueError):
            pass
    return value


def stable_fingerprint(value: Any) -> str:
    """Hash a JSON-compatible context using the same ordering every run."""

    encoded = json.dumps(
        model_to_jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _source_identity(obj: Any, seen: Optional[set[int]] = None) -> Any:
    """Collect source identities for an auditor and its configured children.

    The implementation source hash is intentionally only a deterministic local
    identity, not a cryptographic signature of a trusted package.  Runtime
    monkeypatches and mutable external configuration cannot be inferred from
    source files; callers that need that distinction should pass an explicit
    ``auditor_fingerprint``.
    """

    seen = seen or set()
    marker = id(obj)
    if marker in seen:
        return None
    seen.add(marker)

    if isinstance(obj, (str, int, float, bool, type(None))):
        return obj
    if isinstance(obj, (list, tuple)):
        return [_source_identity(item, seen) for item in obj]
    if isinstance(obj, dict):
        return {
            str(key): _source_identity(value, seen)
            for key, value in sorted(obj.items(), key=lambda item: str(item[0]))
        }

    target = obj if inspect.isclass(obj) or inspect.isfunction(obj) else obj.__class__
    identity: Dict[str, Any] = {
        "module": getattr(target, "__module__", ""),
        "qualname": getattr(target, "__qualname__", getattr(target, "__name__", "")),
    }
    try:
        source = inspect.getsource(target)
    except (OSError, TypeError):
        source = ""
    identity["source_sha256"] = hashlib.sha256(source.encode("utf-8")).hexdigest()

    # ScientificAuditor stores its validators as attributes and in
    # ``additional_validators``.  Include all configured validator classes so a
    # custom extension changes the fingerprint instead of inheriting the core
    # auditor's identity.
    if not inspect.isclass(obj) and not inspect.isfunction(obj):
        children: Dict[str, Any] = {}
        try:
            attributes = vars(obj)
        except TypeError:  # pragma: no cover - slots-based custom auditors
            attributes = {}
        for name, value in attributes.items():
            if name.startswith("_") or name in {"auditor_name", "version"}:
                continue
            if "validator" in name.lower() or name in {"additional_validators"}:
                children[name] = _source_identity(value, seen)
        if children:
            identity["children"] = children
    return identity


def fingerprint_auditor(auditor: Any = None, *, explicit: Optional[str] = None) -> str:
    """Return the stable identity of an auditor implementation.

    ``explicit`` is preferred for deployed/custom auditors.  Otherwise the
    class source and configured validator classes are hashed.  This captures
    normal source changes and additional validators while documenting the
    unavoidable limit around runtime monkeypatches and external configuration.
    """

    if explicit:
        return str(explicit)
    if auditor is None:
        return stable_fingerprint({"auditor": "unknown"})
    return stable_fingerprint(_source_identity(auditor))


def auditor_version(auditor: Any = None, explicit: Optional[str] = None) -> str:
    """Resolve a durable, human-readable auditor version label."""

    if explicit:
        return str(explicit)
    if auditor is not None:
        for name in ("auditor_version", "version", "__version__"):
            value = getattr(auditor, name, None)
            if value:
                return str(value)
        module_name = getattr(auditor.__class__, "__module__", "")
        try:
            package_name = module_name.split(".")[0]
            return str(importlib_metadata.version(package_name))
        except importlib_metadata.PackageNotFoundError:
            pass
    return "unknown"


class ArtifactAuditRecord(BaseModel):
    """Persisted audit receipt for all outputs of one task signature."""

    model_config = ConfigDict(extra="ignore")

    schema_version: int = Field(1, ge=1)
    signature: str = Field(..., min_length=1)
    task_id: str = Field(..., min_length=1)
    contract_fingerprint: str = Field(..., min_length=1)
    # The full contract context is retained to make records inspectable and to
    # allow future readers to verify a context without trusting a short hash.
    contract: Dict[str, Any] = Field(default_factory=dict)
    result_fingerprint: str = Field("", description="Stable task-result identity")
    input_artifacts: list[str] = Field(default_factory=list)
    output_artifacts: list[str] = Field(default_factory=list)
    artifact_hashes: Dict[str, str] = Field(default_factory=dict)
    auditor_name: str = Field(..., min_length=1)
    auditor_version: str = Field(..., min_length=1)
    auditor_fingerprint: str = Field(..., min_length=1)
    status: AuditRecordStatus = AuditRecordStatus.PENDING
    # Denormalized admission flags make index inspection possible without
    # trusting or reparsing the report blob.  ``stop_rule_triggered`` always
    # vetoes ``overall_passed`` in registry finalization.
    overall_passed: bool = False
    stop_rule_triggered: bool = False
    report: Dict[str, Any] = Field(default_factory=dict)
    rejection_reason: Optional[str] = None
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    updated_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))

    @property
    def task_signature(self) -> str:
        """Compatibility name for callers that use task-signature wording."""

        return self.signature

    @property
    def contract_hash(self) -> str:
        return self.contract_fingerprint

    @property
    def output_hashes(self) -> Dict[str, str]:
        return dict(self.artifact_hashes)

    @property
    def passed(self) -> bool:
        return self.status == AuditRecordStatus.PASSED


# Short aliases are useful to integrations and keep the public schema surface
# consistent with the existing ``ArtifactMetadata``/``TaskResult`` naming.
AuditRecord = ArtifactAuditRecord
AuditStatus = AuditRecordStatus


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
