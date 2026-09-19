"""Durable artifact-audit receipts and fail-closed access tests."""

import json
from unittest.mock import patch

import pytest

from eacbp.artifact.registry import ArtifactAuditAccessError, ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _published(tmp_path):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    transaction = TaskArtifactTransaction(registry)
    transaction.register(
        "json://study/output_a/v1",
        {"value": 1},
        ArtifactType.JSON,
        "study",
        "task",
        "write",
    )
    transaction.register(
        "json://study/output_b/v1",
        {"value": 2},
        ArtifactType.JSON,
        "study",
        "task",
        "write",
        parent_uris=["json://study/output_a/v1"],
    )
    result = TaskResult(
        task_id="task",
        capability="capability",
        method_used="method",
        status=TaskStatus.SUCCESS,
        output_artifacts=[
            "json://study/output_a/v1",
            "json://study/output_b/v1",
        ],
    )
    registry.commit_task(transaction, "sig-v1", result)
    contract = TaskContract(
        task_id="task",
        capability="capability",
        expected_outputs=result.output_artifacts,
    )
    return registry, contract, result


def _pass(registry, contract, result):
    registry.begin_audit(
        "sig-v1",
        contract,
        result,
        auditor_name="independent",
        auditor_version="2026.09",
        auditor_fingerprint="auditor-sha",
    )
    report = ValidationReport(
        auditor_name="independent",
        target_task_id="task",
        overall_passed=True,
    )
    return registry.record_audit(
        "sig-v1",
        contract,
        result,
        report,
        auditor_name="independent",
        auditor_version="2026.09",
        auditor_fingerprint="auditor-sha",
    )


def test_audit_is_durable_and_raw_get_remains_unchanged(tmp_path):
    registry, contract, result = _published(tmp_path)

    # Computation publication and raw diagnostics remain available before the
    # scientific gate is completed.
    assert registry.get("json://study/output_a/v1")[1] == {"value": 1}
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited("json://study/output_a/v1", "sig-v1", contract)

    receipt = _pass(registry, contract, result)
    assert receipt.status.value == "passed"
    assert receipt.artifact_hashes["json://study/output_b/v1"].startswith("sha256:")

    restored = ArtifactRegistry(str(tmp_path / "artifacts"))
    metadata, payload = restored.get_audited(
        "json://study/output_a/v1", task_signature="sig-v1", task_contract=contract
    )
    assert metadata.uri.endswith("output_a/v1")
    assert payload == {"value": 1}


def test_reaudit_invalidates_prior_pass_and_stop_rule_cannot_pass(tmp_path):
    registry, contract, result = _published(tmp_path)
    _pass(registry, contract, result)

    # A new attempt invalidates the old pass before the auditor runs.
    pending = registry.begin_audit("sig-v1", contract, result)
    assert pending.status.value == "pending"
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited("json://study/output_a/v1", "sig-v1", contract)


def test_error_check_rejects_even_inconsistent_overall_flag(tmp_path):
    registry, contract, result = _published(tmp_path)
    registry.begin_audit("sig-v1", contract, result)
    report = ValidationReport(
        auditor_name="unknown_auditor",
        target_task_id="task",
        overall_passed=True,
        checks=[
            {
                "check_name": "fatal",
                "passed": False,
                "severity": "error",
                "message": "bad",
            }
        ],
    )
    receipt = registry.record_audit("sig-v1", contract, result, report)
    assert receipt.status.value == "rejected"

    stopped = ValidationReport(
        auditor_name="unknown_auditor",
        target_task_id="task",
        overall_passed=True,
        stop_rule_triggered=True,
    )
    rejected = registry.record_audit("sig-v1", contract, result, stopped)
    assert rejected.status.value == "rejected"
    assert "stop rule" in (rejected.rejection_reason or "")
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited("json://study/output_a/v1", "sig-v1", contract)


def test_audited_access_requires_context_and_rejects_tampered_sibling(tmp_path):
    registry, contract, result = _published(tmp_path)
    _pass(registry, contract, result)

    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited("json://study/output_a/v1", "sig-v1")
    wrong_context = contract.model_copy(update={"method": "different"})
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited("json://study/output_a/v1", "sig-v1", wrong_context)

    sibling = registry.get_metadata("json://study/output_b/v1")
    with open(sibling.storage_path, "w", encoding="utf-8") as handle:
        json.dump({"value": "tampered"}, handle)
    with pytest.raises(ArtifactAuditAccessError):
        # Querying A verifies B too; a partial task cannot be admitted.
        registry.get_audited("json://study/output_a/v1", "sig-v1", contract)


def test_audited_access_can_pin_auditor_identity(tmp_path):
    registry, contract, result = _published(tmp_path)
    _pass(registry, contract, result)
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited(
            "json://study/output_a/v1",
            "sig-v1",
            contract,
            expected_auditor_version="old-version",
        )
    with pytest.raises(ArtifactAuditAccessError):
        registry.get_audited(
            "json://study/output_a/v1",
            "sig-v1",
            contract,
            expected_auditor_fingerprint="old-auditor",
        )
    assert registry.get_audited(
        "json://study/output_a/v1",
        "sig-v1",
        contract,
        auditor_version="2026.09",
        auditor_fingerprint="auditor-sha",
    )[1] == {"value": 1}


def test_legacy_index_without_audit_records_loads_as_unaudited(tmp_path):
    registry, contract, result = _published(tmp_path)
    index_path = registry.storage.base_dir / registry.INDEX_FILENAME
    payload = json.loads(index_path.read_text(encoding="utf-8"))
    payload.pop("audit_records", None)
    index_path.write_text(json.dumps(payload), encoding="utf-8")

    restored = ArtifactRegistry(str(tmp_path / "artifacts"))
    assert restored.list_audit_records() == []
    with pytest.raises(ArtifactAuditAccessError):
        restored.get_audited("json://study/output_a/v1", "sig-v1", contract)


def test_pending_transition_rolls_back_if_index_persist_fails(tmp_path):
    registry, contract, result = _published(tmp_path)
    with patch.object(registry, "_persist_index", side_effect=OSError("disk full")):
        with pytest.raises(OSError):
            registry.begin_audit("sig-v1", contract, result)
    assert registry.list_audit_records() == []
    assert ArtifactRegistry(str(tmp_path / "artifacts")).list_audit_records() == []
