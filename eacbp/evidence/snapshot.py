"""Versioned, provenance checked study snapshots.

The snapshot is a reviewable description of a completed study.  It contains
schemas, metadata, hashes and audit context, but never contains artifact
payloads.  A snapshot is deliberately not an audit receipt: loading or
rendering one never mutates an :class:`ArtifactRegistry` or creates an audit
record.

The public entry points are intentionally small because the command line
front-end uses them as its persistence boundary::

    write_study_snapshot(...)
    load_study_snapshot(path)
    render_snapshot_report(path, artifact_registry)

Validation is performed both when a snapshot is written and when it is read.
The latter matters for a snapshot copied between machines or edited by hand.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from dataclasses import dataclass, field
from datetime import date, datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from pydantic import BaseModel, ValidationError

from eacbp._atomic_json import atomic_write_json
from eacbp.artifact.uri import ArtifactURI
from eacbp.auditor.base import ValidationReport
from eacbp.evidence.graph import EvidenceGraph
from eacbp.schemas.artifact import ArtifactMetadata
from eacbp.schemas.evidence import ClaimNode, EvidenceNode
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskResult, TaskStatus


SNAPSHOT_FORMAT = "eacbp.study_snapshot"
SNAPSHOT_SCHEMA_VERSION = 1
_MAX_INLINE_ARRAY_ELEMENTS = 4096
_EXTERNAL_SOURCE_SCHEMES = {
    "doi",
    "http",
    "https",
    "pmid",
    "pubmed",
    "literature",
}


class SnapshotError(ValueError):
    """Base class for invalid or unusable study snapshots."""


class SnapshotSchemaError(SnapshotError):
    """Raised when the on-disk snapshot schema is unsupported or malformed."""


class SnapshotValidationError(SnapshotError):
    """Raised when cross-object provenance or admission checks fail."""


class SnapshotIntegrityError(SnapshotError):
    """Raised when a live artifact, source file or audit receipt changed."""


class SnapshotWriteError(RuntimeError):
    """Raised when an atomic snapshot publication fails."""


@dataclass(frozen=True)
class ArtifactSource:
    """Metadata and file hashes for one artifact and all known ancestors."""

    uri: str
    sha256: str
    size_bytes: int
    storage_path: str
    parent_uris: List[str] = field(default_factory=list)
    ancestor_uris: List[str] = field(default_factory=list)
    ancestor_hashes: Dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class AuditContextSummary:
    """Small, inspectable audit context retained beside a report."""

    key: str
    task_id: str
    signature: str = ""
    receipt_present: bool = False
    status: str = "unknown"
    overall_passed: bool = False
    stop_rule_triggered: bool = False
    contract_fingerprint: str = ""
    result_fingerprint: str = ""
    input_artifacts: List[str] = field(default_factory=list)
    output_artifacts: List[str] = field(default_factory=list)
    artifact_hashes: Dict[str, str] = field(default_factory=dict)
    auditor_name: str = ""
    auditor_version: str = ""
    auditor_fingerprint: str = ""
    rejection_reason: Optional[str] = None
    report_fingerprint: str = ""
    report: Dict[str, Any] = field(default_factory=dict)
    receipt: Dict[str, Any] = field(default_factory=dict)


@dataclass
class StudySnapshot:
    """Typed in-memory representation of one validated study snapshot."""

    schema_version: int
    manifest: StudyManifest
    config: Any
    summary: Any
    evidence_graph: EvidenceGraph
    task_results: List[TaskResult]
    audit_reports: List[ValidationReport]
    artifact_metadata: Dict[str, ArtifactMetadata]
    artifact_sources: Dict[str, ArtifactSource]
    audit_contexts: Dict[str, AuditContextSummary]
    path: Optional[Path] = None

    @property
    def task_history(self) -> List[TaskResult]:
        """Compatibility name used by the report generator and orchestrator."""

        return self.task_results

    @property
    def run_config(self) -> Any:
        """Compatibility name for callers that call the config run settings."""

        return self.config

    @property
    def artifacts(self) -> Dict[str, ArtifactMetadata]:
        return self.artifact_metadata

    @property
    def audits(self) -> Dict[str, AuditContextSummary]:
        return self.audit_contexts

    def rebuild_evidence_graph(self) -> EvidenceGraph:
        """Return a fresh graph, preserving the validated canonical ordering."""

        return _build_graph_from_models(
            list(self.evidence_graph.evidence_nodes.values()),
            list(self.evidence_graph.claim_nodes.values()),
            _graph_edges(self.evidence_graph),
        )


def _jsonable(value: Any, *, path: str = "") -> Any:
    """Convert values to JSON without silently stringifying scientific data."""

    if isinstance(value, BaseModel):
        try:
            value = value.model_dump(mode="python")
        except TypeError:  # pragma: no cover - pydantic v1 compatibility
            value = value.model_dump()
    if isinstance(value, Enum):
        return _jsonable(value.value, path=path)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise SnapshotValidationError(f"non-finite value at {path or '<root>'}")
        return value
    if isinstance(value, Mapping):
        result: Dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, (str, int, float, bool)):
                raise SnapshotValidationError(
                    f"non-JSON mapping key at {path or '<root>'}: {key!r}"
                )
            key_text = str(key)
            result[key_text] = _jsonable(item, path=f"{path}.{key_text}" if path else key_text)
        return result
    if isinstance(value, (list, tuple, set, frozenset)):
        values = list(value)
        return [
            _jsonable(item, path=f"{path}[{index}]")
            for index, item in enumerate(values)
        ]

    # Numpy is optional at import time.  This branch deliberately refuses a
    # large matrix so a caller cannot accidentally turn a snapshot into a
    # second artifact store.
    ndim = getattr(value, "ndim", None)
    size = getattr(value, "size", None)
    if ndim is not None and size is not None and hasattr(value, "tolist"):
        try:
            if int(ndim) >= 2 and int(size) > _MAX_INLINE_ARRAY_ELEMENTS:
                raise SnapshotValidationError(
                    f"large matrix at {path or '<root>'} cannot be embedded in a snapshot"
                )
            return _jsonable(value.tolist(), path=path)
        except (TypeError, ValueError) as exc:
            if isinstance(exc, SnapshotError):
                raise
            raise SnapshotValidationError(
                f"unsupported array value at {path or '<root>'}: {type(value).__name__}"
            ) from exc
    if hasattr(value, "item") and callable(value.item):
        try:
            return _jsonable(value.item(), path=path)
        except (TypeError, ValueError):
            pass
    raise SnapshotValidationError(
        f"value at {path or '<root>'} is not JSON serializable: {type(value).__name__}"
    )


def _stable_fingerprint(value: Any) -> str:
    payload = json.dumps(
        _jsonable(value),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _model_dict(value: Any, label: str) -> Dict[str, Any]:
    if isinstance(value, BaseModel):
        raw = value.model_dump(mode="python")
    elif isinstance(value, Mapping):
        raw = dict(value)
    else:
        raise SnapshotValidationError(f"{label} must be a model or object mapping")
    if not isinstance(raw, dict):
        raise SnapshotValidationError(f"{label} must serialize to an object")
    return _jsonable(raw, path=label)


def _coerce_model(value: Any, model: Any, label: str) -> Any:
    try:
        if isinstance(value, model):
            return value.model_copy(deep=True)
        if not isinstance(value, Mapping):
            raise SnapshotValidationError(f"{label} must be an object")
        return model.model_validate(value)
    except (ValidationError, TypeError, ValueError) as exc:
        if isinstance(exc, SnapshotError):
            raise
        raise SnapshotValidationError(f"invalid {label}: {exc}") from exc


def _status_value(value: Any) -> str:
    return getattr(value, "value", str(value))


def _external_source(uri: str) -> bool:
    if "://" not in uri:
        return False
    return uri.split("://", 1)[0].casefold() in _EXTERNAL_SOURCE_SCHEMES


def _artifact_uri(value: Any, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise SnapshotValidationError(f"{label} must be a non-empty artifact URI")
    try:
        return ArtifactURI.parse(value).to_string()
    except Exception as exc:
        if _external_source(value):
            return value
        raise SnapshotValidationError(f"invalid artifact URI in {label}: {value!r}") from exc


def _graph_edges(graph: EvidenceGraph) -> List[Dict[str, str]]:
    edges: List[Dict[str, str]] = []
    for source, target, attrs in graph.graph.edges(data=True):
        relation = attrs.get("relationship")
        if relation not in {"supports", "contradicts"}:
            raise SnapshotValidationError(
                f"unsupported evidence graph relationship {relation!r}"
            )
        edges.append({"source": str(source), "target": str(target), "relationship": relation})
    return sorted(edges, key=lambda item: (item["source"], item["target"], item["relationship"]))


def _build_graph_from_models(
    evidence_nodes: Sequence[EvidenceNode],
    claim_nodes: Sequence[ClaimNode],
    edges: Sequence[Mapping[str, str]],
) -> EvidenceGraph:
    graph = EvidenceGraph()
    for node in sorted(evidence_nodes, key=lambda item: item.evidence_id):
        graph.add_evidence(node.model_copy(deep=True))
    for claim in sorted(claim_nodes, key=lambda item: item.claim_id):
        graph.add_claim(claim.model_copy(deep=True))
    # ``add_claim`` only creates edges when its references are known.  The
    # semantic validator has already checked this; this final check prevents a
    # future EvidenceGraph implementation from silently changing the snapshot.
    if _graph_edges(graph) != [dict(item) for item in sorted(edges, key=lambda item: (item["source"], item["target"], item["relationship"]))]:
        raise SnapshotValidationError("evidence graph edge set is inconsistent with claims")
    return graph


def _validate_graph_shape(
    evidence_nodes: Sequence[EvidenceNode],
    claim_nodes: Sequence[ClaimNode],
    edges: Sequence[Mapping[str, str]],
) -> None:
    evidence_ids = [node.evidence_id for node in evidence_nodes]
    claim_ids = [node.claim_id for node in claim_nodes]
    if any(not isinstance(item, str) or not item for item in evidence_ids + claim_ids):
        raise SnapshotValidationError("evidence and claim IDs must be non-empty strings")
    if len(evidence_ids) != len(set(evidence_ids)):
        raise SnapshotValidationError("duplicate evidence ID in snapshot")
    if len(claim_ids) != len(set(claim_ids)):
        raise SnapshotValidationError("duplicate claim ID in snapshot")
    overlap = set(evidence_ids).intersection(claim_ids)
    if overlap:
        raise SnapshotValidationError(f"evidence/claim ID collision: {sorted(overlap)}")

    evidence_set = set(evidence_ids)
    claim_set = set(claim_ids)
    expected: set[Tuple[str, str, str]] = set()
    for claim in claim_nodes:
        support = list(claim.support_evidence_ids)
        contradiction = list(claim.contradiction_evidence_ids)
        if len(support) != len(set(support)) or len(contradiction) != len(set(contradiction)):
            raise SnapshotValidationError(f"duplicate evidence reference in claim {claim.claim_id}")
        if set(support).intersection(contradiction):
            raise SnapshotValidationError(f"evidence is both support and contradiction for {claim.claim_id}")
        if not support:
            raise SnapshotValidationError(f"claim {claim.claim_id} has no supporting evidence")
        for evidence_id in support:
            if evidence_id not in evidence_set:
                raise SnapshotValidationError(
                    f"claim {claim.claim_id} references missing evidence {evidence_id!r}"
                )
            expected.add((evidence_id, claim.claim_id, "supports"))
        for evidence_id in contradiction:
            if evidence_id not in evidence_set:
                raise SnapshotValidationError(
                    f"claim {claim.claim_id} references missing evidence {evidence_id!r}"
                )
            expected.add((evidence_id, claim.claim_id, "contradicts"))

    actual: List[Tuple[str, str, str]] = []
    for edge in edges:
        if set(edge) != {"source", "target", "relationship"}:
            raise SnapshotSchemaError("evidence edge has an invalid shape")
        source, target, relation = edge["source"], edge["target"], edge["relationship"]
        if source not in evidence_set or target not in claim_set:
            raise SnapshotValidationError("evidence graph contains a dangling edge")
        if relation not in {"supports", "contradicts"}:
            raise SnapshotValidationError(f"invalid evidence graph relationship {relation!r}")
        actual.append((source, target, relation))
    if len(actual) != len(set(actual)):
        raise SnapshotValidationError("duplicate evidence graph edge")
    if set(actual) != expected:
        raise SnapshotValidationError("evidence graph edges do not match claim references")


def _registry_artifacts(registry: Any, study_id: str) -> Dict[str, ArtifactMetadata]:
    try:
        values = registry.list_artifacts(study_id=study_id)
    except (AttributeError, TypeError) as exc:
        raise SnapshotValidationError("artifact registry cannot list study metadata") from exc
    result: Dict[str, ArtifactMetadata] = {}
    for value in values:
        metadata = _coerce_model(value, ArtifactMetadata, "artifact metadata")
        try:
            canonical = ArtifactURI.parse(metadata.uri).to_string()
        except Exception as exc:
            raise SnapshotValidationError(f"invalid registered artifact URI {metadata.uri!r}") from exc
        if canonical != metadata.uri:
            raise SnapshotValidationError(f"non-canonical registered artifact URI {metadata.uri!r}")
        if canonical in result:
            raise SnapshotValidationError(f"duplicate artifact URI {canonical!r}")
        result[canonical] = metadata
    return result


def _registry_metadata(registry: Any, uri: str) -> Optional[ArtifactMetadata]:
    try:
        value = registry.get_metadata(uri)
    except (AttributeError, KeyError, ValueError):
        return None
    return _coerce_model(value, ArtifactMetadata, f"artifact metadata {uri}")


def _artifact_refs(
    manifest: StudyManifest,
    tasks: Sequence[TaskResult],
    evidence_nodes: Sequence[EvidenceNode],
    contexts: Sequence[AuditContextSummary],
) -> Iterable[str]:
    if manifest.data.raw_artifact_uri:
        yield manifest.data.raw_artifact_uri
    for task in tasks:
        yield from task.input_artifacts
        yield from task.output_artifacts
    for node in evidence_nodes:
        yield from node.source_artifact_uris
        yield from node.data_origin_uris
    for context in contexts:
        yield from context.input_artifacts
        yield from context.output_artifacts
        yield from context.artifact_hashes


def _collect_artifact_metadata(
    registry: Any,
    manifest: StudyManifest,
    tasks: Sequence[TaskResult],
    evidence_nodes: Sequence[EvidenceNode],
    contexts: Sequence[AuditContextSummary],
) -> Dict[str, ArtifactMetadata]:
    result = _registry_artifacts(registry, manifest.study_id)
    pending = list(result.values())
    for raw_uri in _artifact_refs(manifest, tasks, evidence_nodes, contexts):
        uri = _artifact_uri(raw_uri, "artifact reference")
        if _external_source(uri) or uri in result:
            continue
        metadata = _registry_metadata(registry, uri)
        if metadata is not None:
            result[metadata.uri] = metadata
            pending.append(metadata)
    seen = set()
    while pending:
        metadata = pending.pop()
        if metadata.uri in seen:
            continue
        seen.add(metadata.uri)
        for parent_raw in metadata.parent_uris:
            parent = _artifact_uri(parent_raw, f"parent of {metadata.uri}")
            if _external_source(parent):
                raise SnapshotValidationError(
                    f"artifact {metadata.uri} has an external parent URI; source lineage cannot be verified"
                )
            if parent not in result:
                parent_metadata = _registry_metadata(registry, parent)
                if parent_metadata is None:
                    raise SnapshotValidationError(
                        f"artifact {metadata.uri} has missing parent metadata {parent!r}"
                    )
                result[parent] = parent_metadata
                pending.append(parent_metadata)
    return dict(sorted(result.items()))


def _receipt_values(registry: Any, task_ids: set[str], artifact_uris: set[str]) -> List[Dict[str, Any]]:
    try:
        records = registry.list_audit_records()
    except (AttributeError, TypeError):
        records = list(getattr(registry, "audit_records", {}).values())
    result: List[Dict[str, Any]] = []
    for record in records or []:
        raw = _model_dict(record, "audit receipt")
        task_id = str(raw.get("task_id", ""))
        outputs = set(raw.get("output_artifacts") or [])
        if task_id not in task_ids and not outputs.intersection(artifact_uris):
            continue
        result.append(raw)
    return sorted(result, key=lambda item: str(item.get("signature", "")))


def _audit_contexts(
    reports: Sequence[ValidationReport],
    receipts: Sequence[Mapping[str, Any]],
) -> Dict[str, AuditContextSummary]:
    contexts: Dict[str, AuditContextSummary] = {}
    receipt_tasks: set[str] = set()
    for raw in receipts:
        signature = str(raw.get("signature", ""))
        task_id = str(raw.get("task_id", ""))
        if not signature or not task_id:
            raise SnapshotValidationError("audit receipt requires signature and task_id")
        if signature in contexts:
            raise SnapshotValidationError(f"duplicate audit receipt signature {signature!r}")
        report_raw = raw.get("report")
        if not isinstance(report_raw, Mapping):
            report_raw = {}
        report_dict = _jsonable(dict(report_raw), path=f"audit:{signature}.report")
        context = AuditContextSummary(
            key=signature,
            task_id=task_id,
            signature=signature,
            receipt_present=True,
            status=str(raw.get("status", "unknown")),
            overall_passed=bool(raw.get("overall_passed", False)),
            stop_rule_triggered=bool(raw.get("stop_rule_triggered", False)),
            contract_fingerprint=str(raw.get("contract_fingerprint", "")),
            result_fingerprint=str(raw.get("result_fingerprint", "")),
            input_artifacts=[str(item) for item in raw.get("input_artifacts", []) or []],
            output_artifacts=[str(item) for item in raw.get("output_artifacts", []) or []],
            artifact_hashes={str(key): str(value) for key, value in (raw.get("artifact_hashes", {}) or {}).items()},
            auditor_name=str(raw.get("auditor_name", "")),
            auditor_version=str(raw.get("auditor_version", "")),
            auditor_fingerprint=str(raw.get("auditor_fingerprint", "")),
            rejection_reason=raw.get("rejection_reason"),
            report_fingerprint=_stable_fingerprint(report_dict),
            report=report_dict,
            receipt=dict(_jsonable(dict(raw), path=f"audit:{signature}")),
        )
        contexts[signature] = context
        receipt_tasks.add(task_id)

    for report in reports:
        report_dict = _model_dict(report, f"audit report {report.target_task_id}")
        target = str(report.target_task_id)
        # A durable receipt is authoritative when it exists; retaining a
        # report-only context still lets lightweight integrations use the API
        # without manufacturing a registry receipt.
        if target in receipt_tasks:
            continue
        auditor = str(report_dict.get("auditor_name", "unknown"))
        key = f"report:{target}:{auditor}"
        if key in contexts:
            raise SnapshotValidationError(f"duplicate audit report context {key!r}")
        contexts[key] = AuditContextSummary(
            key=key,
            task_id=target,
            receipt_present=False,
            status="passed" if bool(report.overall_passed) and not bool(report.stop_rule_triggered) else "rejected",
            overall_passed=bool(report.overall_passed),
            stop_rule_triggered=bool(report.stop_rule_triggered),
            auditor_name=auditor,
            report_fingerprint=_stable_fingerprint(report_dict),
            report=report_dict,
        )
    return dict(sorted(contexts.items()))


def _validate_cross_references(
    manifest: StudyManifest,
    tasks: Sequence[TaskResult],
    reports: Sequence[ValidationReport],
    graph: EvidenceGraph,
    artifact_metadata: Mapping[str, ArtifactMetadata],
    artifact_sources: Mapping[str, ArtifactSource],
    contexts: Mapping[str, AuditContextSummary],
) -> None:
    task_ids = [item.task_id for item in tasks]
    if len(task_ids) != len(set(task_ids)):
        raise SnapshotValidationError("duplicate task ID in snapshot")
    task_map = {item.task_id: item for item in tasks}
    for report in reports:
        if report.target_task_id not in task_map:
            raise SnapshotValidationError(
                f"audit report references missing task {report.target_task_id!r}"
            )
    context_list = list(contexts.values())
    context_by_task: Dict[str, List[AuditContextSummary]] = {}
    for context in context_list:
        if context.task_id not in task_map:
            raise SnapshotValidationError(
                f"audit context references missing task {context.task_id!r}"
            )
        context_by_task.setdefault(context.task_id, []).append(context)
        if context.receipt_present:
            derived = _audit_contexts([], [context.receipt]).get(context.key)
            if derived != context:
                raise SnapshotValidationError("audit context differs from its durable receipt")
            from eacbp.schemas.audit import stable_fingerprint
            if stable_fingerprint(context.receipt.get("contract", {})) != context.contract_fingerprint:
                raise SnapshotValidationError("audit contract fingerprint is inconsistent")

    local_artifacts = set(artifact_metadata)
    for metadata in artifact_metadata.values():
        if metadata.uri != metadata.artifact_id and metadata.artifact_id.startswith("://"):
            raise SnapshotValidationError(f"artifact ID is malformed for {metadata.uri!r}")
        for parent in metadata.parent_uris:
            canonical = _artifact_uri(parent, f"parent of {metadata.uri}")
            if canonical not in local_artifacts:
                raise SnapshotValidationError(
                    f"artifact {metadata.uri} has dangling parent {canonical!r}"
                )
    if set(artifact_sources) != local_artifacts:
        raise SnapshotValidationError("artifact source hash coverage differs from metadata coverage")
    for uri, source in artifact_sources.items():
        metadata = artifact_metadata[uri]
        if source.uri != uri or source.sha256 != metadata.sha256_hash:
            raise SnapshotValidationError(f"source hash does not match metadata for {uri}")
        if set(source.ancestor_uris) != set(source.ancestor_hashes):
            raise SnapshotValidationError(f"ancestor hash coverage differs for {uri}")
        if not set(source.ancestor_uris).issubset(local_artifacts):
            raise SnapshotValidationError(f"artifact {uri} has a dangling ancestor hash")
        for ancestor, digest in source.ancestor_hashes.items():
            if digest != artifact_metadata[ancestor].sha256_hash:
                raise SnapshotValidationError(f"ancestor hash mismatch for {ancestor}")

    for task in tasks:
        if _status_value(task.status) != TaskStatus.SUCCESS.value:
            # The task itself is retained for a reviewable failure, but any
            # evidence pointing at it is rejected below.
            continue
        for raw_uri in list(task.input_artifacts) + list(task.output_artifacts):
            uri = _artifact_uri(raw_uri, f"task {task.task_id} artifact")
            if not _external_source(uri) and uri not in local_artifacts:
                raise SnapshotValidationError(
                    f"task {task.task_id} references missing artifact {uri!r}"
                )

    reports_by_task: Dict[str, List[ValidationReport]] = {}
    for report in reports:
        reports_by_task.setdefault(report.target_task_id, []).append(report)

    def report_is_admitted(report: ValidationReport) -> bool:
        if not bool(report.overall_passed) or bool(report.stop_rule_triggered):
            return False
        for check in report.checks:
            severity = _status_value(check.severity)
            if not bool(check.passed) and severity in {"error", "stop_rule"}:
                return False
        return True

    def context_is_admitted(context: AuditContextSummary) -> bool:
        if not context.receipt_present:
            return False
        if context.status != "passed" or not context.overall_passed or context.stop_rule_triggered:
            return False
        report_payload = context.report
        if not bool(report_payload.get("overall_passed", False)) or bool(report_payload.get("stop_rule_triggered", False)):
            return False
        for check in report_payload.get("checks", []) or []:
            if isinstance(check, Mapping):
                severity = _status_value(check.get("severity"))
                if check.get("passed") is False and severity in {"error", "stop_rule"}:
                    return False
        return True

    def admitted(task_id: str, evidence: Optional[EvidenceNode] = None) -> bool:
        task = task_map[task_id]
        if _status_value(task.status) != TaskStatus.SUCCESS.value:
            return False
        task_reports = reports_by_task.get(task_id, [])
        if task_reports and any(not report_is_admitted(report) for report in task_reports):
            return False
        task_contexts = context_by_task.get(task_id, [])
        # A snapshot report cannot manufacture a durable scientific audit.  A
        # report-only context is useful for retaining a failure explanation,
        # but it never admits evidence.  Match a passing receipt to the task's
        # output set (and, when called for a node, to that node's sources).
        from eacbp.artifact.registry import ArtifactRegistry
        passing_contexts = [context for context in task_contexts if context_is_admitted(context)
                            and context.result_fingerprint == ArtifactRegistry._result_identity(task.model_dump(mode="python"))
                            and set(context.output_artifacts) == set(task.output_artifacts)
                            and set(context.artifact_hashes) == set(task.output_artifacts)
                            and all(uri in artifact_metadata and artifact_metadata[uri].sha256_hash == digest
                                    for uri, digest in context.artifact_hashes.items())]
        if not passing_contexts:
            return False
        if evidence is not None:
            source_uris = set(evidence.source_artifact_uris)
            for context in passing_contexts:
                if source_uris.issubset(set(context.output_artifacts)):
                    # A rejected receipt for the same output is a revocation,
                    # even if another historical receipt passed.
                    revoked = any(
                        not context_is_admitted(other)
                        and source_uris.intersection(other.output_artifacts)
                        for other in task_contexts
                    )
                    if not revoked:
                        return True
            return False
        return True

    for node in graph.evidence_nodes.values():
        if node.source_task_id not in task_map:
            raise SnapshotValidationError(
                f"evidence {node.evidence_id} references missing source task {node.source_task_id!r}"
            )
        task = task_map[node.source_task_id]
        if _status_value(task.status) != TaskStatus.SUCCESS.value:
            raise SnapshotValidationError(
                f"evidence {node.evidence_id} is sourced from failed task {node.source_task_id!r}"
            )
        if not node.source_artifact_uris:
            raise SnapshotValidationError(
                f"evidence {node.evidence_id} has no source artifact URI"
            )
        task_outputs = set()
        for raw_output in task.output_artifacts:
            task_outputs.add(_artifact_uri(raw_output, f"task {task.task_id} output"))
        for raw_uri in node.source_artifact_uris:
            uri = _artifact_uri(raw_uri, f"evidence {node.evidence_id} source artifact")
            if _external_source(uri) or uri not in task_outputs:
                raise SnapshotValidationError(
                    f"evidence {node.evidence_id} source artifact {uri!r} is not a task output"
                )
            if uri not in local_artifacts:
                raise SnapshotValidationError(
                    f"evidence {node.evidence_id} references missing artifact {uri!r}"
                )
            if artifact_metadata[uri].created_by_task != node.source_task_id:
                raise SnapshotValidationError(
                    f"evidence {node.evidence_id} source artifact {uri!r} was created by a different task"
                )
        for raw_uri in node.data_origin_uris:
            uri = _artifact_uri(raw_uri, f"evidence {node.evidence_id} artifact")
            if not _external_source(uri) and uri not in local_artifacts:
                raise SnapshotValidationError(
                    f"evidence {node.evidence_id} references missing artifact {uri!r}"
                )
        if not node.audit_passed or not admitted(node.source_task_id, node):
            raise SnapshotValidationError(
                f"evidence {node.evidence_id} is marked audited without a passing audit context"
            )

    for task in tasks:
        if task.status == TaskStatus.SUCCESS and not admitted(task.task_id):
            raise SnapshotValidationError(f"successful task {task.task_id} has no matching durable audit")

    linked_ids: set[str] = set()
    for claim in graph.claim_nodes.values():
        linked_ids.update(claim.support_evidence_ids)
        linked_ids.update(claim.contradiction_evidence_ids)
    for evidence_id in linked_ids:
        node = graph.evidence_nodes[evidence_id]
        if not node.audit_passed or not admitted(node.source_task_id, node):
            raise SnapshotValidationError(
                f"claim-linked evidence {evidence_id} is not admitted by a passing audit"
            )


def _artifact_source_records(
    registry: Any,
    metadata: Mapping[str, ArtifactMetadata],
) -> Dict[str, ArtifactSource]:
    records: Dict[str, ArtifactSource] = {}
    lineage = getattr(registry, "lineage", None)
    for uri, item in sorted(metadata.items()):
        ancestor_uris: set[str] = set()
        if lineage is not None:
            try:
                ancestor_uris.update(str(value) for value in lineage.get_ancestors(uri))
            except Exception:
                pass
        # Parent metadata is authoritative even if a legacy lineage graph had
        # no placeholder nodes.
        queue = list(item.parent_uris)
        while queue:
            parent = str(queue.pop())
            if parent in ancestor_uris:
                continue
            ancestor_uris.add(parent)
            if parent in metadata:
                queue.extend(metadata[parent].parent_uris)
        if not ancestor_uris.issubset(metadata):
            raise SnapshotValidationError(f"artifact {uri} has an unresolved ancestor")
        records[uri] = ArtifactSource(
            uri=uri,
            sha256=item.sha256_hash,
            size_bytes=int(item.size_bytes),
            storage_path=str(item.storage_path),
            parent_uris=sorted(item.parent_uris),
            ancestor_uris=sorted(ancestor_uris),
            ancestor_hashes={
                ancestor: metadata[ancestor].sha256_hash
                for ancestor in sorted(ancestor_uris)
            },
        )
    return records


def _artifact_payload_hash(registry: Any, metadata: ArtifactMetadata) -> Tuple[str, int]:
    storage = getattr(registry, "storage", None)
    raw_path = Path(metadata.storage_path)
    if not raw_path.is_absolute() and storage is not None:
        raw_path = Path(storage.base_dir) / raw_path
    try:
        if storage is not None and hasattr(storage, "_validated_path"):
            path = storage._validated_path(raw_path)
        else:
            path = raw_path.expanduser().resolve(strict=False)
    except (OSError, ValueError) as exc:
        raise SnapshotIntegrityError(
            f"artifact storage path is invalid for {metadata.uri}: {metadata.storage_path}"
        ) from exc
    if not path.exists() or not path.is_file():
        raise SnapshotIntegrityError(f"artifact payload is missing for {metadata.uri}: {path}")
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
                size += len(chunk)
    except OSError as exc:
        raise SnapshotIntegrityError(f"unable to read artifact payload for {metadata.uri}") from exc
    return f"sha256:{digest.hexdigest()}", size


def _verify_registry_scope(snapshot: StudySnapshot, registry: Any) -> None:
    live = _registry_artifacts(registry, snapshot.manifest.study_id)
    expected_study = {
        uri for uri, metadata in snapshot.artifact_metadata.items()
        if metadata.study_id == snapshot.manifest.study_id
    }
    if set(live) != expected_study:
        added = sorted(set(live) - expected_study)
        missing = sorted(expected_study - set(live))
        raise SnapshotIntegrityError(
            f"live artifact registry scope changed (added={added}, missing={missing})"
        )

    for uri, expected in sorted(snapshot.artifact_metadata.items()):
        current = _registry_metadata(registry, uri)
        if current is None:
            raise SnapshotIntegrityError(f"artifact metadata is missing for {uri}")
        if _model_dict(current, f"artifact metadata {uri}") != _model_dict(expected, f"snapshot artifact {uri}"):
            raise SnapshotIntegrityError(f"artifact metadata changed for {uri}")
        actual_hash, actual_size = _artifact_payload_hash(registry, current)
        source = snapshot.artifact_sources[uri]
        if actual_hash != source.sha256 or actual_size != source.size_bytes:
            raise SnapshotIntegrityError(f"artifact payload changed for {uri}")
        if current.sha256_hash != expected.sha256_hash:
            raise SnapshotIntegrityError(f"artifact source hash changed for {uri}")
        for ancestor, digest in source.ancestor_hashes.items():
            ancestor_meta = _registry_metadata(registry, ancestor)
            if ancestor_meta is None or ancestor_meta.sha256_hash != digest:
                raise SnapshotIntegrityError(f"artifact ancestor source changed for {ancestor}")

    for key, context in snapshot.audit_contexts.items():
        if not context.receipt_present:
            continue
        try:
            current = registry.get_audit_record(context.signature)
        except Exception as exc:
            raise SnapshotIntegrityError(
                f"audit receipt cannot be read for {context.signature}"
            ) from exc
        if current is None:
            raise SnapshotIntegrityError(f"audit receipt is missing for {context.signature}")
        current_payload = _model_dict(current, f"live audit receipt {context.signature}")
        expected_payload = context.receipt
        fields = (
            "signature", "task_id", "contract_fingerprint", "result_fingerprint",
            "input_artifacts", "output_artifacts", "artifact_hashes", "auditor_name",
            "auditor_version", "auditor_fingerprint", "status", "overall_passed",
            "stop_rule_triggered", "rejection_reason", "report",
        )
        if any(current_payload.get(field) != expected_payload.get(field) for field in fields):
            raise SnapshotIntegrityError(f"audit receipt changed for {context.signature}")
        commit = registry.get_task_commit(context.signature)
        if commit is None or commit.get("output_hashes") != context.artifact_hashes:
            raise SnapshotIntegrityError(f"computation receipt changed for {context.signature}")
        if context.status != "passed" or not context.overall_passed or context.stop_rule_triggered:
            # A previously rejected receipt is never upgraded by a later live
            # registry pass; the snapshot remains a historical failure.
            if current_payload.get("status") == "passed":
                raise SnapshotIntegrityError(
                    f"rejected audit receipt {context.signature} cannot be reactivated"
                )


def _snapshot_payload(snapshot: StudySnapshot) -> Dict[str, Any]:
    evidence_nodes = [
        _model_dict(node, f"evidence node {node.evidence_id}")
        for node in sorted(snapshot.evidence_graph.evidence_nodes.values(), key=lambda item: item.evidence_id)
    ]
    claims = [
        _model_dict(node, f"claim node {node.claim_id}")
        for node in sorted(snapshot.evidence_graph.claim_nodes.values(), key=lambda item: item.claim_id)
    ]
    metadata = [
        _model_dict(item, f"artifact metadata {uri}")
        for uri, item in sorted(snapshot.artifact_metadata.items())
    ]
    sources = []
    for uri, source in sorted(snapshot.artifact_sources.items()):
        sources.append({
            "uri": source.uri,
            "sha256": source.sha256,
            "size_bytes": source.size_bytes,
            "storage_path": source.storage_path,
            "parent_uris": list(source.parent_uris),
            "ancestor_uris": list(source.ancestor_uris),
            "ancestor_hashes": dict(sorted(source.ancestor_hashes.items())),
        })
    contexts = []
    for key, context in sorted(snapshot.audit_contexts.items()):
        contexts.append({
            "key": key,
            "task_id": context.task_id,
            "signature": context.signature,
            "receipt_present": context.receipt_present,
            "status": context.status,
            "overall_passed": context.overall_passed,
            "stop_rule_triggered": context.stop_rule_triggered,
            "contract_fingerprint": context.contract_fingerprint,
            "result_fingerprint": context.result_fingerprint,
            "input_artifacts": list(context.input_artifacts),
            "output_artifacts": list(context.output_artifacts),
            "artifact_hashes": dict(sorted(context.artifact_hashes.items())),
            "auditor_name": context.auditor_name,
            "auditor_version": context.auditor_version,
            "auditor_fingerprint": context.auditor_fingerprint,
            "rejection_reason": context.rejection_reason,
            "report_fingerprint": context.report_fingerprint,
            "report": context.report,
            "receipt": context.receipt,
        })
    payload = {
        "format": SNAPSHOT_FORMAT,
        "schema_version": SNAPSHOT_SCHEMA_VERSION,
        "manifest": _model_dict(snapshot.manifest, "manifest"),
        "config": _jsonable(snapshot.config, path="config"),
        "summary": _jsonable(snapshot.summary, path="summary"),
        "evidence": {
            "nodes": evidence_nodes,
            "claims": claims,
            "edges": _graph_edges(snapshot.evidence_graph),
        },
        "task_results": [
            _model_dict(item, f"task result {item.task_id}")
            for item in sorted(snapshot.task_results, key=lambda value: value.task_id)
        ],
        "audits": {
            "reports": [
                _model_dict(item, f"audit report {item.target_task_id}")
                for item in sorted(snapshot.audit_reports, key=lambda value: (value.target_task_id, value.auditor_name))
            ],
            "contexts": contexts,
        },
        "artifacts": {
            "metadata": metadata,
            "sources": sources,
            "source_hashes": {uri: source.sha256 for uri, source in sorted(snapshot.artifact_sources.items())},
        },
    }
    payload["fingerprints"] = {
        section: _stable_fingerprint(payload[section])
        for section in ("manifest", "config", "summary", "evidence", "task_results", "audits", "artifacts")
    }
    return payload


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path = Path(path)
    if not path.name:
        raise SnapshotWriteError("snapshot path must name a file")
    try:
        atomic_write_json(
            path,
            payload,
            create_parent=True,
            dump_kwargs={
                "indent": 2,
                "sort_keys": True,
                "ensure_ascii": False,
                "allow_nan": False,
            },
            newline="\n",
            trailing_newline=True,
        )
    except (OSError, TypeError, ValueError) as exc:
        raise SnapshotWriteError(f"unable to atomically write snapshot {path}: {exc}") from exc
    return path


def write_study_snapshot(
    path: str | Path,
    *,
    manifest: StudyManifest,
    config: Any,
    summary: Any,
    evidence_graph: EvidenceGraph,
    artifact_registry: Any,
    task_history: Sequence[TaskResult],
    audit_reports: Sequence[ValidationReport],
) -> Path:
    """Validate and atomically write a versioned study snapshot.

    The registry is read for metadata, lineage and durable audit context.  No
    artifact payload is deserialized and no audit receipt is written.
    """

    manifest_model = _coerce_model(manifest, StudyManifest, "manifest")
    if not isinstance(evidence_graph, EvidenceGraph):
        raise SnapshotValidationError("evidence_graph must be an EvidenceGraph")
    tasks = [_coerce_model(item, TaskResult, "task result") for item in task_history]
    reports = [_coerce_model(item, ValidationReport, "audit report") for item in audit_reports]
    evidence_nodes = list(evidence_graph.evidence_nodes.values())
    claims = list(evidence_graph.claim_nodes.values())
    edges = _graph_edges(evidence_graph)
    _validate_graph_shape(evidence_nodes, claims, edges)

    # Build audit contexts before collecting artifact references, because a
    # receipt's input/output hash set is part of the source coverage.
    task_ids = {task.task_id for task in tasks}
    preliminary_artifacts = _registry_artifacts(artifact_registry, manifest_model.study_id)
    receipts = _receipt_values(artifact_registry, task_ids, set(preliminary_artifacts))
    contexts = _audit_contexts(reports, receipts)
    metadata = _collect_artifact_metadata(
        artifact_registry,
        manifest_model,
        tasks,
        evidence_nodes,
        list(contexts.values()),
    )
    sources = _artifact_source_records(artifact_registry, metadata)
    snapshot = StudySnapshot(
        schema_version=SNAPSHOT_SCHEMA_VERSION,
        manifest=manifest_model,
        config=_jsonable(config, path="config"),
        summary=_jsonable(summary, path="summary"),
        evidence_graph=_build_graph_from_models(evidence_nodes, claims, edges),
        task_results=sorted(tasks, key=lambda item: item.task_id),
        audit_reports=sorted(reports, key=lambda item: (item.target_task_id, item.auditor_name)),
        artifact_metadata=dict(sorted(metadata.items())),
        artifact_sources=dict(sorted(sources.items())),
        audit_contexts=contexts,
    )
    _validate_cross_references(
        snapshot.manifest,
        snapshot.task_results,
        snapshot.audit_reports,
        snapshot.evidence_graph,
        snapshot.artifact_metadata,
        snapshot.artifact_sources,
        snapshot.audit_contexts,
    )
    _verify_registry_scope(snapshot, artifact_registry)
    return _atomic_write_json(Path(path), _snapshot_payload(snapshot))


def _require_keys(value: Any, required: set[str], allowed: Optional[set[str]] = None, label: str = "object") -> None:
    if not isinstance(value, Mapping):
        raise SnapshotSchemaError(f"{label} must be an object")
    keys = set(value)
    if not required.issubset(keys):
        raise SnapshotSchemaError(f"{label} is missing keys: {sorted(required - keys)}")
    if allowed is not None and not keys.issubset(allowed):
        raise SnapshotSchemaError(f"{label} has unknown keys: {sorted(keys - allowed)}")


def _load_contexts(raw_values: Any) -> Dict[str, AuditContextSummary]:
    if not isinstance(raw_values, list):
        raise SnapshotSchemaError("audits.contexts must be a list")
    contexts: Dict[str, AuditContextSummary] = {}
    allowed = {
        "key", "task_id", "signature", "receipt_present", "status", "overall_passed",
        "stop_rule_triggered", "contract_fingerprint", "result_fingerprint",
        "input_artifacts", "output_artifacts", "artifact_hashes", "auditor_name",
        "auditor_version", "auditor_fingerprint", "rejection_reason", "report_fingerprint",
        "report", "receipt",
    }
    for index, raw in enumerate(raw_values):
        _require_keys(raw, {"key", "task_id", "receipt_present", "status", "overall_passed", "stop_rule_triggered", "report", "receipt"}, allowed, f"audit context {index}")
        key = raw["key"]
        if not isinstance(key, str) or not key:
            raise SnapshotSchemaError(f"audit context {index} has invalid key")
        if key in contexts:
            raise SnapshotValidationError(f"duplicate audit context key {key!r}")
        if any(type(raw[name]) is not bool for name in ("receipt_present", "overall_passed", "stop_rule_triggered")):
            raise SnapshotSchemaError("audit admission flags must be booleans")
        try:
            contexts[key] = AuditContextSummary(
                key=key,
                task_id=str(raw["task_id"]),
                signature=str(raw.get("signature", "")),
                receipt_present=bool(raw["receipt_present"]),
                status=str(raw["status"]),
                overall_passed=bool(raw["overall_passed"]),
                stop_rule_triggered=bool(raw["stop_rule_triggered"]),
                contract_fingerprint=str(raw.get("contract_fingerprint", "")),
                result_fingerprint=str(raw.get("result_fingerprint", "")),
                input_artifacts=[str(item) for item in raw.get("input_artifacts", []) or []],
                output_artifacts=[str(item) for item in raw.get("output_artifacts", []) or []],
                artifact_hashes={str(key): str(value) for key, value in (raw.get("artifact_hashes", {}) or {}).items()},
                auditor_name=str(raw.get("auditor_name", "")),
                auditor_version=str(raw.get("auditor_version", "")),
                auditor_fingerprint=str(raw.get("auditor_fingerprint", "")),
                rejection_reason=raw.get("rejection_reason"),
                report_fingerprint=str(raw.get("report_fingerprint", "")),
                report=dict(raw["report"]),
                receipt=dict(raw["receipt"]),
            )
        except (TypeError, ValueError) as exc:
            raise SnapshotSchemaError(f"invalid audit context {key!r}") from exc
        if contexts[key].receipt_present and not contexts[key].signature:
            raise SnapshotValidationError(f"durable audit context {key!r} has no signature")
    return dict(sorted(contexts.items()))


def load_study_snapshot(path: str | Path) -> StudySnapshot:
    """Load and strictly validate a snapshot without touching artifact payloads."""

    path = Path(path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            raw = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise SnapshotSchemaError(f"unable to read snapshot {path}: {exc}") from exc
    top_allowed = {
        "format", "schema_version", "manifest", "config", "summary", "evidence",
        "task_results", "audits", "artifacts", "fingerprints",
    }
    _require_keys(raw, top_allowed, top_allowed, "snapshot")
    if raw.get("format") != SNAPSHOT_FORMAT:
        raise SnapshotSchemaError(f"unsupported snapshot format: {raw.get('format')!r}")
    version = raw.get("schema_version")
    if isinstance(version, bool) or not isinstance(version, int) or version != SNAPSHOT_SCHEMA_VERSION:
        raise SnapshotSchemaError(
            f"unsupported snapshot schema version {version!r}; expected {SNAPSHOT_SCHEMA_VERSION}"
        )
    fingerprints = raw.get("fingerprints")
    _require_keys(
        fingerprints,
        {"manifest", "config", "summary", "evidence", "task_results", "audits", "artifacts"},
        {"manifest", "config", "summary", "evidence", "task_results", "audits", "artifacts"},
        "snapshot fingerprints",
    )
    for section in ("manifest", "config", "summary", "evidence", "task_results", "audits", "artifacts"):
        expected_fingerprint = fingerprints[section]
        if not isinstance(expected_fingerprint, str) or expected_fingerprint != _stable_fingerprint(raw[section]):
            raise SnapshotIntegrityError(f"snapshot {section} fingerprint does not match its contents")
    manifest = _coerce_model(raw["manifest"], StudyManifest, "manifest")

    evidence_raw = raw["evidence"]
    _require_keys(evidence_raw, {"nodes", "claims", "edges"}, {"nodes", "claims", "edges"}, "evidence")
    if not all(isinstance(evidence_raw[key], list) for key in ("nodes", "claims", "edges")):
        raise SnapshotSchemaError("evidence nodes, claims and edges must be lists")
    evidence_nodes: List[EvidenceNode] = []
    claims: List[ClaimNode] = []
    for index, item in enumerate(evidence_raw["nodes"]):
        evidence_nodes.append(_coerce_model(item, EvidenceNode, f"evidence node {index}"))
    for index, item in enumerate(evidence_raw["claims"]):
        claims.append(_coerce_model(item, ClaimNode, f"claim node {index}"))
    edges: List[Dict[str, str]] = []
    for index, item in enumerate(evidence_raw["edges"]):
        _require_keys(item, {"source", "target", "relationship"}, {"source", "target", "relationship"}, f"evidence edge {index}")
        if not all(isinstance(item[key], str) for key in ("source", "target", "relationship")):
            raise SnapshotSchemaError(f"evidence edge {index} has non-string fields")
        edges.append(dict(item))
    _validate_graph_shape(evidence_nodes, claims, edges)
    graph = _build_graph_from_models(evidence_nodes, claims, edges)

    tasks_raw = raw["task_results"]
    if not isinstance(tasks_raw, list):
        raise SnapshotSchemaError("task_results must be a list")
    tasks = [_coerce_model(item, TaskResult, f"task result {index}") for index, item in enumerate(tasks_raw)]

    audits_raw = raw["audits"]
    _require_keys(audits_raw, {"reports", "contexts"}, {"reports", "contexts"}, "audits")
    if not isinstance(audits_raw["reports"], list):
        raise SnapshotSchemaError("audits.reports must be a list")
    reports = [_coerce_model(item, ValidationReport, f"audit report {index}") for index, item in enumerate(audits_raw["reports"])]
    contexts = _load_contexts(audits_raw["contexts"])

    artifacts_raw = raw["artifacts"]
    _require_keys(artifacts_raw, {"metadata", "sources", "source_hashes"}, {"metadata", "sources", "source_hashes"}, "artifacts")
    if not isinstance(artifacts_raw["metadata"], list) or not isinstance(artifacts_raw["sources"], list) or not isinstance(artifacts_raw["source_hashes"], Mapping):
        raise SnapshotSchemaError("artifacts.metadata/sources/source_hashes have invalid shapes")
    metadata: Dict[str, ArtifactMetadata] = {}
    for index, item in enumerate(artifacts_raw["metadata"]):
        value = _coerce_model(item, ArtifactMetadata, f"artifact metadata {index}")
        if value.uri in metadata:
            raise SnapshotValidationError(f"duplicate artifact URI {value.uri!r}")
        metadata[value.uri] = value
    source_allowed = {"uri", "sha256", "size_bytes", "storage_path", "parent_uris", "ancestor_uris", "ancestor_hashes"}
    sources: Dict[str, ArtifactSource] = {}
    for index, item in enumerate(artifacts_raw["sources"]):
        _require_keys(item, source_allowed, source_allowed, f"artifact source {index}")
        uri = item.get("uri")
        if not isinstance(uri, str) or not uri:
            raise SnapshotSchemaError(f"artifact source {index} has invalid URI")
        if uri in sources:
            raise SnapshotValidationError(f"duplicate artifact source URI {uri!r}")
        if not isinstance(item.get("ancestor_hashes"), Mapping):
            raise SnapshotSchemaError(f"artifact source {uri} has invalid ancestor_hashes")
        sources[uri] = ArtifactSource(
            uri=uri,
            sha256=str(item["sha256"]),
            size_bytes=int(item["size_bytes"]),
            storage_path=str(item["storage_path"]),
            parent_uris=[str(value) for value in item["parent_uris"]],
            ancestor_uris=[str(value) for value in item["ancestor_uris"]],
            ancestor_hashes={str(key): str(value) for key, value in item["ancestor_hashes"].items()},
        )
    source_hashes = {str(key): str(value) for key, value in artifacts_raw["source_hashes"].items()}
    if source_hashes != {uri: source.sha256 for uri, source in sources.items()}:
        raise SnapshotValidationError("artifact source_hashes index does not match source records")

    snapshot = StudySnapshot(
        schema_version=version,
        manifest=manifest,
        config=_jsonable(raw["config"], path="config"),
        summary=_jsonable(raw["summary"], path="summary"),
        evidence_graph=graph,
        task_results=sorted(tasks, key=lambda item: item.task_id),
        audit_reports=sorted(reports, key=lambda item: (item.target_task_id, item.auditor_name)),
        artifact_metadata=dict(sorted(metadata.items())),
        artifact_sources=dict(sorted(sources.items())),
        audit_contexts=contexts,
        path=path,
    )
    _validate_cross_references(
        snapshot.manifest,
        snapshot.task_results,
        snapshot.audit_reports,
        snapshot.evidence_graph,
        snapshot.artifact_metadata,
        snapshot.artifact_sources,
        snapshot.audit_contexts,
    )
    return snapshot


def render_snapshot_report(path: str | Path, artifact_registry: Any) -> str:
    """Validate a snapshot against live files/receipts and render its report."""

    snapshot = load_study_snapshot(path)
    _verify_registry_scope(snapshot, artifact_registry)
    from eacbp.report.markdown_report import ScientificReportGenerator

    # The generator is deliberately given only the snapshot's task/evidence
    # state.  It may read payloads for table previews, but every payload and
    # metadata record was hash-checked above and the study registry scope is
    # frozen, so a new artifact cannot be mixed into an old report.
    generator = ScientificReportGenerator(
        manifest=snapshot.manifest,
        evidence_graph=snapshot.evidence_graph,
        artifact_registry=artifact_registry,
        task_history=snapshot.task_results,
        audit_reports=snapshot.audit_reports,
    )
    try:
        return generator.generate_markdown()
    except (OSError, KeyError, ValueError) as exc:
        raise SnapshotIntegrityError(f"unable to render validated snapshot report: {exc}") from exc


__all__ = [
    "SNAPSHOT_FORMAT",
    "SNAPSHOT_SCHEMA_VERSION",
    "SnapshotError",
    "SnapshotSchemaError",
    "SnapshotValidationError",
    "SnapshotIntegrityError",
    "SnapshotWriteError",
    "ArtifactSource",
    "AuditContextSummary",
    "StudySnapshot",
    "write_study_snapshot",
    "load_study_snapshot",
    "render_snapshot_report",
]
