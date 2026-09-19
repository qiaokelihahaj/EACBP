"""Persistent artifact registry with integrity-checked payload retrieval."""

from __future__ import annotations

import json
import sys
import threading
from functools import wraps
from copy import deepcopy
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from eacbp._atomic_json import atomic_write_json
from eacbp.artifact.lineage import LineageGraph
from eacbp.artifact.storage import (
    ArtifactAlreadyExistsError,
    ArtifactStorageBackend,
)
from eacbp.artifact.uri import ArtifactURI
from eacbp.schemas.artifact import ArtifactMetadata, ArtifactType
from eacbp.schemas.audit import (
    ArtifactAuditRecord,
    AuditRecordStatus,
    auditor_version as resolve_auditor_version,
    fingerprint_auditor,
    model_to_jsonable,
    stable_fingerprint,
)


class ArtifactRegistryError(RuntimeError):
    """Raised when the durable registry index is invalid or cannot be written."""


class ArtifactAuditError(ArtifactRegistryError):
    """Raised when an audit receipt cannot be created or persisted."""


class ArtifactAuditAccessError(ArtifactAuditError, PermissionError):
    """Raised when an audited query cannot prove admissibility."""


# A descriptive alias for integrations that distinguish a missing receipt from
# a rejected/tampered receipt at their boundary.  All audited access remains
# fail closed regardless of the specific reason.
AuditedArtifactAccessError = ArtifactAuditAccessError


def _synchronized(method):
    @wraps(method)
    def call(self, *args, **kwargs):
        with self._mutex:
            return method(self, *args, **kwargs)
    return call


def _index_json_default(value: Any) -> Any:
    """Serialize metadata values without the lossy ``default=str`` fallback."""

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (pd.Timestamp, datetime, date)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (set, frozenset)):
        return list(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


class ArtifactRegistry:
    """Manage immutable payloads, metadata, and their lineage across restarts."""

    INDEX_FILENAME = ".artifact_registry.json"
    INDEX_VERSION = 1

    def __init__(self, storage_dir: str = ".artifacts"):
        self._mutex = threading.RLock()
        self.storage = ArtifactStorageBackend(base_dir=storage_dir)
        self._index_path = self.storage.base_dir / self.INDEX_FILENAME
        self.lineage = LineageGraph()
        self.registry: Dict[str, ArtifactMetadata] = {}
        self.task_commits = {}
        # Audit records are deliberately separate from task computation
        # receipts.  A computation can be durably published while its audit is
        # still pending or rejected, and raw ``get`` remains available for
        # recovery/diagnostics in either state.
        self.audit_records: Dict[str, ArtifactAuditRecord] = {}
        self._refresh()

    def _refresh(self) -> None:
        # Windows readers can prevent atomic replacement of an open index.
        # Coordinate reads with writers, including initial construction.
        with self.storage.lock():
            self._load_persisted()

    def _rebuild_lineage(self, lineage_data: Optional[Dict[str, Any]] = None) -> None:
        graph = LineageGraph()
        for metadata in self.registry.values():
            graph.add_artifact(metadata)

        # Parent URIs are the canonical source of lineage.  Retain persisted
        # placeholder nodes/edges as well so a registry can rebuild a graph even
        # when an upstream metadata record was produced by an older run.
        if lineage_data:
            for node in lineage_data.get("nodes", []):
                uri = node.get("uri")
                attrs = node.get("attributes", {})
                if isinstance(uri, str):
                    graph.graph.add_node(uri, **attrs)
            for edge in lineage_data.get("edges", []):
                source = edge.get("source")
                target = edge.get("target")
                if isinstance(source, str) and isinstance(target, str):
                    graph.graph.add_edge(
                        source,
                        target,
                        **{k: v for k, v in edge.items() if k not in {"source", "target"}},
                    )
        # Metadata remains authoritative if the persisted graph contains stale
        # node attributes after a partial older write.
        for metadata in self.registry.values():
            graph.metadata_store[metadata.uri] = metadata
        self.lineage = graph

    def _load_persisted(self) -> None:
        if not self._index_path.exists():
            self.registry.clear()
            self.task_commits = {}
            self.audit_records = {}
            self.lineage = LineageGraph()
            return
        try:
            with self._index_path.open("r", encoding="utf-8") as handle:
                payload = json.load(handle)
            if not isinstance(payload, dict):
                raise ValueError("registry index root must be an object")
            records = payload.get("metadata", [])
            if not isinstance(records, list):
                raise ValueError("registry index metadata must be a list")
            rebuilt: Dict[str, ArtifactMetadata] = {}
            for record in records:
                metadata = ArtifactMetadata.model_validate(record)
                canonical_uri = ArtifactURI.parse(metadata.uri).to_string()
                if canonical_uri != metadata.uri:
                    raise ValueError(
                        f"non-canonical artifact URI in registry index: {metadata.uri!r}"
                    )
                rebuilt[metadata.uri] = metadata
            self.registry = rebuilt
            self.task_commits = payload.get("task_commits", {})
            if not isinstance(self.task_commits, dict):
                raise ValueError("registry index task_commits must be an object")
            # ``audit_records`` was added after index version 1.  Missing data
            # is a valid legacy index and intentionally means "unaudited";
            # it must never be inferred as a pass.
            persisted_audits = payload.get("audit_records", {})
            if persisted_audits is None:
                persisted_audits = {}
            if not isinstance(persisted_audits, dict):
                raise ValueError("registry index audit_records must be an object")
            rebuilt_audits: Dict[str, ArtifactAuditRecord] = {}
            for key, record in persisted_audits.items():
                if not isinstance(key, str) or not isinstance(record, Mapping):
                    raise ValueError("registry index audit_records entries must be objects")
                audit = ArtifactAuditRecord.model_validate(record)
                if audit.signature != key:
                    raise ValueError(
                        f"audit record key/signature mismatch for {key!r}"
                    )
                rebuilt_audits[key] = audit
            self.audit_records = rebuilt_audits
            self._rebuild_lineage(payload.get("lineage"))
        except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
            raise ArtifactRegistryError(
                f"Unable to rebuild artifact registry from {self._index_path}: {exc}"
            ) from exc

    def _lineage_payload(self) -> Dict[str, Any]:
        nodes = []
        for uri, attrs in self.lineage.graph.nodes(data=True):
            nodes.append({"uri": uri, "attributes": dict(attrs)})
        edges = []
        for source, target, attrs in self.lineage.graph.edges(data=True):
            edges.append({"source": source, "target": target, **dict(attrs)})
        return {"nodes": nodes, "edges": edges}

    def _index_payload(self) -> Dict[str, Any]:
        metadata_records = [
            metadata.model_dump(mode="json")
            for metadata in sorted(self.registry.values(), key=lambda item: item.uri)
        ]
        return {
            "index_version": self.INDEX_VERSION,
            "metadata": metadata_records,
            "lineage": self._lineage_payload(),
            "task_commits": self.task_commits,
            "audit_records": {
                signature: record.model_dump(mode="json")
                for signature, record in sorted(self.audit_records.items())
            },
        }

    @_synchronized
    def get_task_commit(self, signature):
        self._refresh()
        return deepcopy(self.task_commits.get(signature))

    # ------------------------------------------------------------------
    # Durable audit receipts
    # ------------------------------------------------------------------
    @staticmethod
    def _field(value: Any, name: str, default: Any = None) -> Any:
        if isinstance(value, Mapping):
            return value.get(name, default)
        return getattr(value, name, default)

    @staticmethod
    def _contract_payload(contract: Any) -> Dict[str, Any]:
        if contract is None:
            raise ArtifactAuditError(
                "An explicit task contract/context is required for artifact auditing"
            )
        payload = model_to_jsonable(contract)
        if not isinstance(payload, dict):
            raise ArtifactAuditError("Task contract/context must serialize to an object")
        return deepcopy(payload)

    @classmethod
    def _result_payload(cls, result: Any) -> Dict[str, Any]:
        if result is None:
            raise ArtifactAuditError("An explicit task result is required for artifact auditing")
        payload = model_to_jsonable(result)
        if not isinstance(payload, dict):
            raise ArtifactAuditError("Task result must serialize to an object")
        return deepcopy(payload)

    @classmethod
    def _result_identity(cls, result_payload: Mapping[str, Any]) -> str:
        """Fingerprint result identity without mutable logs or summary metrics."""

        fields = {
            key: result_payload.get(key)
            for key in (
                "task_id",
                "capability",
                "method_used",
                "status",
                "input_artifacts",
                "output_artifacts",
            )
        }
        return stable_fingerprint(fields)

    @staticmethod
    def _canonical_output_uris(values: Any) -> List[str]:
        uris: List[str] = []
        for value in values or []:
            if not isinstance(value, str):
                raise ArtifactAuditError("Task output artifact URI must be a string")
            uri = ArtifactURI.parse(value).to_string()
            if uri not in uris:
                uris.append(uri)
        return uris

    @staticmethod
    def _status_value(value: Any) -> str:
        raw = getattr(value, "value", value)
        return str(raw).lower()

    def _commit_hashes(self, signature: str) -> Dict[str, str]:
        receipt = self.task_commits.get(signature)
        if not isinstance(receipt, Mapping):
            raise ArtifactAuditError(
                f"Task signature {signature!r} has no published computation receipt"
            )
        raw_hashes = receipt.get("output_hashes", {})
        if not isinstance(raw_hashes, Mapping):
            raise ArtifactAuditError(
                f"Task signature {signature!r} has invalid output hash receipt"
            )
        hashes: Dict[str, str] = {}
        for raw_uri, digest in raw_hashes.items():
            if not isinstance(raw_uri, str) or not isinstance(digest, str):
                raise ArtifactAuditError(
                    f"Task signature {signature!r} has an invalid output hash entry"
                )
            hashes[ArtifactURI.parse(raw_uri).to_string()] = digest
        return hashes

    def _snapshot_audit_hashes(
        self,
        signature: str,
        result_payload: Mapping[str, Any],
    ) -> Tuple[List[str], Dict[str, str]]:
        """Snapshot every committed sibling output for one audit context.

        The task receipt is authoritative for the output set.  A result with an
        explicit output list must agree with that set; otherwise a caller could
        audit one output and silently leave a sibling outside the audit.
        """

        commit_hashes = self._commit_hashes(signature)
        committed_outputs = list(commit_hashes)
        result_outputs = self._canonical_output_uris(result_payload.get("output_artifacts", []))
        if result_outputs and set(result_outputs) != set(committed_outputs):
            raise ArtifactAuditError(
                f"Task result outputs do not match published receipt for {signature!r}"
            )
        outputs = result_outputs or committed_outputs
        # Preserve receipt hashes, even if a metadata record is missing or the
        # payload was tampered with.  This allows ``begin_audit`` to durably
        # invalidate an earlier pass before fail-closed validation reports the
        # problem.
        hashes = {uri: commit_hashes.get(uri, "") for uri in outputs}
        return outputs, hashes

    @staticmethod
    def _auditor_name(auditor: Any = None, explicit: Optional[str] = None) -> str:
        if explicit:
            return str(explicit)
        if auditor is not None:
            value = getattr(auditor, "auditor_name", None) or getattr(auditor, "name", None)
            if value:
                return str(value)
        return "unknown_auditor"

    def _new_pending_audit_locked(
        self,
        signature: str,
        contract: Any,
        result: Any,
        *,
        auditor: Any = None,
        auditor_name: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> ArtifactAuditRecord:
        contract_payload = self._contract_payload(contract)
        result_payload = self._result_payload(result)
        outputs, artifact_hashes = self._snapshot_audit_hashes(signature, result_payload)
        task_id = self._field(contract_payload, "task_id") or self._field(result_payload, "task_id")
        if not isinstance(task_id, str) or not task_id:
            raise ArtifactAuditError("Task contract/result must include a task_id")
        result_task_id = self._field(result_payload, "task_id")
        if result_task_id and result_task_id != task_id:
            raise ArtifactAuditError("Task contract and result task identities do not match")
        receipt = self.task_commits.get(signature)
        receipt_result = receipt.get("result", {}) if isinstance(receipt, Mapping) else {}
        receipt_task_id = self._field(receipt_result, "task_id")
        if receipt_task_id and receipt_task_id != task_id:
            raise ArtifactAuditError(
                "Task contract identity does not match the published computation receipt"
            )
        input_artifacts = self._canonical_output_uris(result_payload.get("input_artifacts", []))
        now = datetime.now(timezone.utc)
        return ArtifactAuditRecord(
            signature=signature,
            task_id=task_id,
            contract_fingerprint=stable_fingerprint(contract_payload),
            contract=contract_payload,
            result_fingerprint=self._result_identity(result_payload),
            input_artifacts=input_artifacts,
            output_artifacts=outputs,
            artifact_hashes=artifact_hashes,
            auditor_name=self._auditor_name(auditor, auditor_name),
            auditor_version=resolve_auditor_version(auditor, explicit=auditor_version),
            auditor_fingerprint=fingerprint_auditor(auditor, explicit=auditor_fingerprint),
            status=AuditRecordStatus.PENDING,
            created_at=now,
            updated_at=now,
        )

    def _replace_audit_locked(self, record: ArtifactAuditRecord) -> ArtifactAuditRecord:
        """Persist one audit state transition with in-memory rollback."""

        previous = self.audit_records.get(record.signature)
        try:
            self.audit_records[record.signature] = record
            self._persist_index()
        except BaseException:
            if previous is None:
                self.audit_records.pop(record.signature, None)
            else:
                self.audit_records[record.signature] = previous
            raise
        return record.model_copy(deep=True)

    def _begin_audit_locked(
        self,
        signature: str,
        contract: Any,
        result: Any,
        *,
        auditor: Any = None,
        auditor_name: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> ArtifactAuditRecord:
        pending = self._new_pending_audit_locked(
            signature,
            contract,
            result,
            auditor=auditor,
            auditor_name=auditor_name,
            auditor_version=auditor_version,
            auditor_fingerprint=auditor_fingerprint,
        )
        return self._replace_audit_locked(pending)

    @_synchronized
    def begin_audit(
        self,
        signature: str,
        contract: Any = None,
        result: Any = None,
        *,
        task_contract: Any = None,
        task_result: Any = None,
        auditor: Any = None,
        auditor_name: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> ArtifactAuditRecord:
        """Atomically transition a task signature into ``pending`` audit.

        ``begin_audit`` must be called after ``commit_task``.  It intentionally
        writes before the auditor executes, replacing any prior ``passed``
        receipt so a crash or exception during re-audit cannot reuse stale
        scientific evidence.
        """

        contract = contract if contract is not None else task_contract
        result = result if result is not None else task_result
        signature = str(signature)
        if not signature:
            raise ArtifactAuditError("An audit signature is required")
        with self.storage.lock():
            self._load_persisted()
            return self._begin_audit_locked(
                signature,
                contract,
                result,
                auditor=auditor,
                auditor_name=auditor_name,
                auditor_version=auditor_version,
                auditor_fingerprint=auditor_fingerprint,
            )

    def _verify_audit_outputs_locked(self, record: ArtifactAuditRecord) -> Tuple[bool, str]:
        """Verify all sibling outputs, not only the artifact being queried."""

        try:
            receipt_hashes = self._commit_hashes(record.signature)
        except ArtifactAuditError as exc:
            return False, str(exc)
        if set(receipt_hashes) != set(record.artifact_hashes):
            return False, "audit output set differs from published computation receipt"
        if set(record.output_artifacts) != set(record.artifact_hashes):
            return False, "audit record does not cover every published sibling output"
        for uri in record.artifact_hashes:
            expected = record.artifact_hashes.get(uri, "")
            if not expected or receipt_hashes.get(uri) != expected:
                return False, f"audit hash receipt is incomplete for {uri}"
            metadata = self.registry.get(uri)
            if metadata is None:
                return False, f"audited output metadata is missing for {uri}"
            if metadata.sha256_hash != expected:
                return False, f"audited output metadata hash changed for {uri}"
            try:
                # ``load`` verifies the digest and also rejects a payload that
                # can no longer be deserialized after tampering.
                self.storage.load(
                    metadata.uri,
                    metadata.type,
                    expected_sha256=expected,
                    storage_path=metadata.storage_path,
                )
            except Exception as exc:  # fail closed for every storage error
                return False, f"audited output integrity check failed for {uri}: {exc}"
        return True, ""

    @staticmethod
    def _report_payload(report: Any) -> Dict[str, Any]:
        if report is None:
            return {}
        payload = model_to_jsonable(report)
        return payload if isinstance(payload, dict) else {}

    @staticmethod
    def _report_field(report_payload: Mapping[str, Any], name: str, default: Any = None) -> Any:
        return report_payload.get(name, default)

    @_synchronized
    def record_audit(
        self,
        signature: str,
        contract: Any = None,
        result: Any = None,
        report: Any = None,
        *,
        task_contract: Any = None,
        task_result: Any = None,
        auditor: Any = None,
        auditor_name: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> ArtifactAuditRecord:
        """Durably record a validation report for a published computation.

        If a caller omits ``begin_audit`` this method starts a pending state
        itself.  That convenience path still invalidates a previous pass before
        finalization; orchestration code should call ``begin_audit`` explicitly
        before running the auditor so crashes leave a visible pending receipt.
        A rejected report is persisted as ``rejected`` and never admitted by
        ``get_audited``.
        """

        contract = contract if contract is not None else task_contract
        result = result if result is not None else task_result
        signature = str(signature)
        if not signature:
            raise ArtifactAuditError("An audit signature is required")
        contract_payload = self._contract_payload(contract)
        result_payload = self._result_payload(result)
        contract_fp = stable_fingerprint(contract_payload)
        result_fp = self._result_identity(result_payload)
        report_payload = self._report_payload(report)

        with self.storage.lock():
            self._load_persisted()
            pending = self.audit_records.get(signature)
            if (
                pending is None
                or pending.status != AuditRecordStatus.PENDING
                or pending.contract_fingerprint != contract_fp
                or pending.result_fingerprint != result_fp
            ):
                # Persist the pending transition before inspecting the report,
                # including when replacing an earlier passed record.
                pending = self._begin_audit_locked(
                    signature,
                    contract,
                    result,
                    auditor=auditor,
                    auditor_name=auditor_name,
                    auditor_version=auditor_version,
                    auditor_fingerprint=auditor_fingerprint,
                )

            reasons: List[str] = []
            expected_task_id = self._field(contract_payload, "task_id")
            report_task_id = self._report_field(report_payload, "target_task_id")
            if not report_payload:
                reasons.append("audit report is missing")
            if report_task_id != expected_task_id:
                reasons.append("audit report target task does not match the contract")
            if self._status_value(result_payload.get("status")) != "success":
                reasons.append("only successful computations can become audited evidence")
            if not bool(self._report_field(report_payload, "overall_passed", False)):
                reasons.append("auditor report rejected the computation")
            # A stop rule is a hard veto even if a custom report accidentally
            # leaves overall_passed=True.
            if bool(self._report_field(report_payload, "stop_rule_triggered", False)):
                reasons.append("auditor stop rule was triggered")
            try:
                observed_outputs = self._canonical_output_uris(
                    result_payload.get("output_artifacts", [])
                )
                if observed_outputs and set(observed_outputs) != set(pending.output_artifacts):
                    reasons.append("task result outputs do not match the audited receipt")
            except ArtifactAuditError as exc:
                reasons.append(str(exc))
            for check in self._report_field(report_payload, "checks", []) or []:
                if not isinstance(check, Mapping):
                    continue
                severity = self._status_value(check.get("severity"))
                if check.get("passed") is False and severity in {"error", "stop_rule"}:
                    reasons.append(
                        f"auditor check rejected the computation: {check.get('check_name', 'unknown')}"
                    )

            observed_name = self._report_field(report_payload, "auditor_name")
            final_name = self._auditor_name(auditor, auditor_name) if (auditor or auditor_name) else pending.auditor_name
            if observed_name and final_name not in {"unknown_auditor", observed_name}:
                reasons.append("audit report auditor identity does not match the receipt")
            if not observed_name:
                reasons.append("audit report auditor identity is missing")
            if pending.auditor_name != "unknown_auditor" and observed_name and observed_name != pending.auditor_name:
                reasons.append("audit report auditor identity changed during audit")

            final_version = (
                resolve_auditor_version(auditor, explicit=auditor_version)
                if auditor is not None or auditor_version
                else pending.auditor_version
            )
            final_fingerprint = (
                fingerprint_auditor(auditor, explicit=auditor_fingerprint)
                if auditor is not None or auditor_fingerprint
                else pending.auditor_fingerprint
            )
            if final_version != pending.auditor_version and pending.auditor_version != "unknown":
                reasons.append("auditor version changed during audit")
            if final_fingerprint != pending.auditor_fingerprint and pending.auditor_fingerprint != fingerprint_auditor(None):
                reasons.append("auditor fingerprint changed during audit")

            integrity_ok, integrity_reason = self._verify_audit_outputs_locked(pending)
            if not integrity_ok:
                reasons.append(integrity_reason)

            status = AuditRecordStatus.PASSED if not reasons else AuditRecordStatus.REJECTED
            now = datetime.now(timezone.utc)
            report_passed = bool(self._report_field(report_payload, "overall_passed", False))
            report_stop = bool(self._report_field(report_payload, "stop_rule_triggered", False))
            updated = pending.model_copy(
                update={
                    "auditor_name": observed_name or final_name,
                    "auditor_version": final_version,
                    "auditor_fingerprint": final_fingerprint,
                    "status": status,
                    "overall_passed": report_passed and not report_stop and not reasons,
                    "stop_rule_triggered": report_stop,
                    "report": deepcopy(report_payload),
                    "rejection_reason": "; ".join(reasons) if reasons else None,
                    "updated_at": now,
                }
            )
            return self._replace_audit_locked(updated)

    # ``finalize_audit`` is a semantic alias for orchestration adapters.
    finalize_audit = record_audit

    @_synchronized
    def get_audit_record(
        self,
        signature: str,
        contract: Any = None,
        *,
        task_contract: Any = None,
        context: Any = None,
    ) -> Optional[ArtifactAuditRecord]:
        """Read an audit receipt, optionally requiring a matching context."""

        contract = contract if contract is not None else task_contract
        contract = contract if contract is not None else context
        signature = str(signature)
        self._refresh()
        record = self.audit_records.get(signature)
        if record is None:
            return None
        if contract is not None:
            expected = stable_fingerprint(self._contract_payload(contract))
            if expected != record.contract_fingerprint:
                raise ArtifactAuditAccessError(
                    "Audit receipt context does not match the requested task contract"
                )
        return record.model_copy(deep=True)

    @_synchronized
    def list_audit_records(
        self,
        status: Optional[AuditRecordStatus | str] = None,
    ) -> List[ArtifactAuditRecord]:
        self._refresh()
        expected = self._status_value(status) if status is not None else None
        records = [
            record
            for record in self.audit_records.values()
            if expected is None or record.status.value == expected
        ]
        return [record.model_copy(deep=True) for record in sorted(records, key=lambda item: item.signature)]

    def _require_audited_locked(
        self,
        signature: str,
        contract: Any,
        *,
        expected_auditor_version: Optional[str] = None,
        expected_auditor_fingerprint: Optional[str] = None,
    ) -> ArtifactAuditRecord:
        if contract is None:
            raise ArtifactAuditAccessError(
                "get_audited requires an explicit task contract/context"
            )
        record = self.audit_records.get(signature)
        if record is None:
            raise ArtifactAuditAccessError(
                f"No durable audit receipt exists for task signature {signature!r}"
            )
        if record.status != AuditRecordStatus.PASSED:
            raise ArtifactAuditAccessError(
                f"Task signature {signature!r} is not admitted: {record.status.value}"
            )
        if not record.overall_passed or record.stop_rule_triggered:
            raise ArtifactAuditAccessError(
                "Durable audit admission flags do not certify this artifact"
            )
        contract_payload = self._contract_payload(contract)
        if stable_fingerprint(contract_payload) != record.contract_fingerprint:
            raise ArtifactAuditAccessError(
                "Audited artifact context does not match the requested task contract"
            )
        if (
            expected_auditor_version is not None
            and str(expected_auditor_version) != record.auditor_version
        ):
            raise ArtifactAuditAccessError(
                "Audited artifact was produced by a different auditor version"
            )
        if (
            expected_auditor_fingerprint is not None
            and str(expected_auditor_fingerprint) != record.auditor_fingerprint
        ):
            raise ArtifactAuditAccessError(
                "Audited artifact was produced by a different auditor implementation"
            )
        if self._field(contract_payload, "task_id") != record.task_id:
            raise ArtifactAuditAccessError("Audited artifact task identity does not match the contract")
        report_payload = record.report
        if self._report_field(report_payload, "target_task_id") != record.task_id:
            raise ArtifactAuditAccessError("Durable audit report target does not match the task")
        if not bool(self._report_field(report_payload, "overall_passed", False)):
            raise ArtifactAuditAccessError("Durable audit report did not pass")
        if bool(self._report_field(report_payload, "stop_rule_triggered", False)):
            raise ArtifactAuditAccessError("Durable audit stop rule prevents evidence access")
        integrity_ok, reason = self._verify_audit_outputs_locked(record)
        if not integrity_ok:
            raise ArtifactAuditAccessError(reason)
        return record

    @_synchronized
    def get_audited(
        self,
        uri_str: str,
        signature: Optional[str] = None,
        contract: Any = None,
        *,
        task_signature: Optional[str] = None,
        task_contract: Any = None,
        context: Any = None,
        expected_auditor_version: Optional[str] = None,
        expected_auditor_fingerprint: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> Tuple[ArtifactMetadata, Any]:
        """Return an output only after a matching durable audit passes.

        The explicit signature and contract/context are required by design;
        this method never searches for an arbitrary historical pass.  All
        sibling outputs in the task receipt are hash-checked before the
        requested payload is returned.
        """

        signature = signature if signature is not None else task_signature
        contract = contract if contract is not None else task_contract
        contract = contract if contract is not None else context
        if signature is None:
            raise ArtifactAuditAccessError("get_audited requires an explicit task signature")
        signature = str(signature)
        expected_auditor_version = (
            expected_auditor_version
            if expected_auditor_version is not None
            else auditor_version
        )
        expected_auditor_fingerprint = (
            expected_auditor_fingerprint
            if expected_auditor_fingerprint is not None
            else auditor_fingerprint
        )
        with self.storage.lock():
            self._load_persisted()
            record = self._require_audited_locked(
                signature,
                contract,
                expected_auditor_version=expected_auditor_version,
                expected_auditor_fingerprint=expected_auditor_fingerprint,
            )
            canonical_uri = ArtifactURI.parse(uri_str).to_string()
            if canonical_uri not in record.output_artifacts:
                raise ArtifactAuditAccessError(
                    f"Artifact {canonical_uri!r} is not an output of audited task {signature!r}"
                )
            metadata = self.registry.get(canonical_uri)
            if metadata is None:
                raise ArtifactAuditAccessError(f"Audited artifact metadata is missing for {canonical_uri}")
            try:
                payload = self.storage.load(
                    metadata.uri,
                    metadata.type,
                    expected_sha256=record.artifact_hashes[canonical_uri],
                    storage_path=metadata.storage_path,
                )
            except Exception as exc:
                raise ArtifactAuditAccessError(
                    f"Audited artifact integrity check failed for {canonical_uri}: {exc}"
                ) from exc
            return metadata.model_copy(deep=True), payload

    # A query spelling is provided for integrations that use query terminology;
    # both paths enforce the same explicit signature/context requirements.
    query_audited = get_audited

    def get_audited_payload(self, *args, **kwargs) -> Any:
        return self.get_audited(*args, **kwargs)[1]

    load_audited_payload = get_audited_payload

    @_synchronized
    def get_audited_many(
        self,
        signature: Optional[str] = None,
        contract: Any = None,
        *,
        task_signature: Optional[str] = None,
        task_contract: Any = None,
        context: Any = None,
        expected_auditor_version: Optional[str] = None,
        expected_auditor_fingerprint: Optional[str] = None,
        auditor_version: Optional[str] = None,
        auditor_fingerprint: Optional[str] = None,
    ) -> List[Tuple[ArtifactMetadata, Any]]:
        """Return every output of an audited task after one full verification."""

        signature = signature if signature is not None else task_signature
        contract = contract if contract is not None else task_contract
        contract = contract if contract is not None else context
        if signature is None:
            raise ArtifactAuditAccessError("get_audited_many requires an explicit task signature")
        signature = str(signature)
        expected_auditor_version = (
            expected_auditor_version
            if expected_auditor_version is not None
            else auditor_version
        )
        expected_auditor_fingerprint = (
            expected_auditor_fingerprint
            if expected_auditor_fingerprint is not None
            else auditor_fingerprint
        )
        with self.storage.lock():
            self._load_persisted()
            record = self._require_audited_locked(
                signature,
                contract,
                expected_auditor_version=expected_auditor_version,
                expected_auditor_fingerprint=expected_auditor_fingerprint,
            )
            values: List[Tuple[ArtifactMetadata, Any]] = []
            for uri in record.output_artifacts:
                metadata = self.registry[uri]
                payload = self.storage.load(
                    metadata.uri,
                    metadata.type,
                    expected_sha256=record.artifact_hashes[uri],
                    storage_path=metadata.storage_path,
                )
                values.append((metadata.model_copy(deep=True), payload))
            return values

    query_audited_many = get_audited_many

    @_synchronized
    def commit_task(self, staged, signature, result):
        """Publish all task metadata and its recovery receipt in one index swap."""
        records = staged.list_artifacts()
        # Read-only tasks may return existing inputs as their outputs.
        existing_outputs = [self.get(uri)[0] for uri in
                            set(result.output_artifacts) - {m.uri for m in records}]
        for meta in records:
            staged.get(meta.uri)  # Verify every payload before publication.
            self.storage._validated_path(Path(meta.storage_path))
        with self.storage.lock():
            self._load_persisted()
            if signature in self.task_commits:
                raise ArtifactRegistryError("Task transaction is already committed")
            collisions = set(self.registry) & {m.uri for m in records}
            if collisions:
                raise ArtifactAlreadyExistsError(f"Transaction output collision: {sorted(collisions)}")
            try:
                for meta in records:
                    self.registry[meta.uri] = meta
                self._rebuild_lineage()
                self.task_commits[signature] = {
                    "signature": signature, "phase": "computed", "result": result.model_dump(mode="json"),
                    "output_hashes": {m.uri: m.sha256_hash for m in records + existing_outputs},
                }
                self._persist_index()
            except BaseException:
                self._load_persisted()
                raise

    def _persist_index(self) -> None:
        """Atomically replace the registry index after a successful payload write."""

        try:
            atomic_write_json(
                self._index_path,
                self._index_payload(),
                dump_kwargs={
                    "indent": 2,
                    "ensure_ascii": False,
                    "default": _index_json_default,
                },
            )
        except (OSError, TypeError, ValueError) as exc:
            raise ArtifactRegistryError(
                f"Unable to persist artifact registry index {self._index_path}: {exc}"
            ) from exc

    def _inherit_provenance(
        self,
        parent_uris: List[str],
        summary_metrics: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Carry simulation/origin markers through every derived artifact."""

        parent_metadata = [self.registry.get(uri) for uri in parent_uris]
        parent_metadata = [metadata for metadata in parent_metadata if metadata is not None]
        simulated_values = [
            metadata.summary_metrics["is_simulated"]
            for metadata in parent_metadata
            if "is_simulated" in metadata.summary_metrics
        ]
        if simulated_values:
            # A derived artifact is simulated if any contributing parent is
            # simulated; an explicit false cannot erase that provenance.
            summary_metrics["is_simulated"] = bool(summary_metrics.get("is_simulated", False)) or any(
                bool(value) for value in simulated_values
            )

        origins = [
            metadata.summary_metrics["data_origin"]
            for metadata in parent_metadata
            if "data_origin" in metadata.summary_metrics
        ]
        if origins and "data_origin" not in summary_metrics:
            unique_origins = []
            for origin in origins:
                if origin not in unique_origins:
                    unique_origins.append(origin)
            summary_metrics["data_origin"] = (
                unique_origins[0] if len(unique_origins) == 1 else unique_origins
            )
        return summary_metrics

    @staticmethod
    def _payload_provenance(payload: Any) -> Dict[str, Any]:
        """Extract simulation markers from a root SCData/AnnData payload."""

        if isinstance(payload, dict):
            uns = payload.get("uns", {})
        else:
            uns = getattr(payload, "uns", {})
        if uns is None or not hasattr(uns, "get"):
            return {}
        markers: Dict[str, Any] = {}
        for key in ("is_simulated", "data_origin"):
            if key in uns:
                markers[key] = deepcopy(uns[key])
        return markers

    @_synchronized
    def register(
        self,
        uri_str: str,
        payload: Any,
        artifact_type: ArtifactType,
        study_id: str,
        created_by_task: str,
        operation: str,
        parent_uris: Optional[List[str]] = None,
        parameters: Optional[Dict[str, Any]] = None,
        software_versions: Optional[Dict[str, str]] = None,
        random_seed: int = 42,
        summary_metrics: Optional[Dict[str, Any]] = None,
        container_image: Optional[str] = None,
    ) -> ArtifactMetadata:
        parsed_uri = ArtifactURI.parse(uri_str)
        canonical_uri = parsed_uri.to_string()
        if study_id != parsed_uri.study_id:
            raise ValueError(
                f"study_id {study_id!r} does not match artifact URI study {parsed_uri.study_id!r}"
            )
        canonical_parents = [ArtifactURI.parse(parent).to_string() for parent in (parent_uris or [])]
        if canonical_uri in canonical_parents:
            raise ValueError("An artifact cannot list itself as a parent")
        parameters = deepcopy(parameters or {})
        software_versions = deepcopy(
            software_versions
            or {
                "python": sys.version.split()[0],
                "eacbp": "0.1.0",
            }
        )
        initial_summary = deepcopy(summary_metrics or {})
        payload_markers = self._payload_provenance(payload)
        for key, value in payload_markers.items():
            initial_summary.setdefault(key, value)
        inherited_summary = self._inherit_provenance(canonical_parents, initial_summary)

        with self.storage.lock():
            # A registry object may have been created before another process or
            # instance registered an artifact.  Re-read the atomic index while
            # holding the same lock used by writers before checking uniqueness.
            self._load_persisted()
            inherited_summary = self._inherit_provenance(
                canonical_parents,
                deepcopy(inherited_summary),
            )
            if canonical_uri in self.registry:
                raise ArtifactAlreadyExistsError(
                    f"Artifact '{canonical_uri}' is already registered"
                )
            for parent_uri in canonical_parents:
                if parent_uri in self.lineage.graph and parent_uri in self.lineage.get_descendants(canonical_uri):
                    raise ValueError(
                        f"Adding '{canonical_uri}' under parent '{parent_uri}' would create a lineage cycle"
                    )

            storage_path: Optional[Path] = None
            try:
                storage_path, sha256_hash, size_bytes = self.storage.save(
                    uri_str=canonical_uri,
                    payload=payload,
                    artifact_type=artifact_type,
                )
                metadata = ArtifactMetadata(
                    artifact_id=canonical_uri,
                    uri=canonical_uri,
                    type=artifact_type,
                    study_id=study_id,
                    parent_uris=canonical_parents,
                    created_by_task=created_by_task,
                    operation=operation,
                    parameters=parameters,
                    software_versions=software_versions,
                    container_image=container_image,
                    random_seed=random_seed,
                    sha256_hash=sha256_hash,
                    storage_path=str(storage_path),
                    size_bytes=size_bytes,
                    summary_metrics=inherited_summary,
                )
                self.registry[canonical_uri] = metadata
                self.lineage.add_artifact(metadata)
                try:
                    self._persist_index()
                except Exception:
                    self.registry.pop(canonical_uri, None)
                    self._rebuild_lineage()
                    if storage_path is not None:
                        self.storage.remove(storage_path)
                    raise
                return metadata.model_copy(deep=True)
            except Exception:
                # PayloadSerializer already cleans its temporary file.  If the
                # payload was published but metadata construction/index commit
                # failed, remove only that exact target, never the study tree.
                if storage_path is not None and canonical_uri not in self.registry:
                    try:
                        self.storage.remove(storage_path)
                    except OSError:
                        pass
                raise

    @_synchronized
    def get_metadata(self, uri_str: str) -> ArtifactMetadata:
        canonical_uri = ArtifactURI.parse(uri_str).to_string()
        self._refresh()
        if canonical_uri not in self.registry:
            raise KeyError(f"Artifact URI '{canonical_uri}' is not registered in the system.")
        return self.registry[canonical_uri].model_copy(deep=True)

    def load_payload(self, uri_str: str) -> Any:
        meta = self.get_metadata(uri_str)
        return self.storage.load(
            meta.uri,
            meta.type,
            expected_sha256=meta.sha256_hash,
            storage_path=meta.storage_path,
        )

    def get(self, uri_str: str) -> Tuple[ArtifactMetadata, Any]:
        meta = self.get_metadata(uri_str)
        payload = self.storage.load(
            meta.uri,
            meta.type,
            expected_sha256=meta.sha256_hash,
            storage_path=meta.storage_path,
        )
        return meta, payload

    @_synchronized
    def exists(self, uri_str: str, artifact_type: Optional[ArtifactType] = None) -> bool:
        canonical_uri = ArtifactURI.parse(uri_str).to_string()
        self._refresh()
        metadata = self.registry.get(canonical_uri)
        if metadata is None:
            return False
        if artifact_type is not None and metadata.type != artifact_type:
            return False
        try:
            return self.storage._validated_path(Path(metadata.storage_path)).exists()
        except (OSError, ValueError):
            return False

    @_synchronized
    def list_artifacts(
        self,
        study_id: Optional[str] = None,
        artifact_type: Optional[ArtifactType] = None,
    ) -> List[ArtifactMetadata]:
        self._refresh()
        results = list(self.registry.values())
        if study_id:
            results = [artifact for artifact in results if artifact.study_id == study_id]
        if artifact_type:
            results = [artifact for artifact in results if artifact.type == artifact_type]
        return [artifact.model_copy(deep=True) for artifact in results]
