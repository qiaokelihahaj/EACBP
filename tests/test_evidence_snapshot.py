"""Focused tests for the versioned evidence/report snapshot boundary."""

import json
from pathlib import Path

import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.snapshot import (
    SnapshotIntegrityError,
    SnapshotSchemaError,
    SnapshotValidationError,
    SnapshotWriteError,
    load_study_snapshot,
    render_snapshot_report,
    write_study_snapshot,
)
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import (
    ClaimNode,
    ClaimType,
    ConfidenceScore,
    EvidenceNode,
    EvidenceType,
    LanguageTier,
)
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _fixture(tmp_path: Path, *, audit_passed: bool = True, task_status=TaskStatus.SUCCESS):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    registry.register(
        "json://snap/raw/v1",
        {"raw": True},
        ArtifactType.JSON,
        "snap",
        "task_import",
        "raw_import",
    )
    transaction = TaskArtifactTransaction(registry)
    transaction.register(
        "json://snap/result/v1",
        {"value": 1},
        ArtifactType.JSON,
        "snap",
        "task_result",
        "compute",
        parent_uris=["json://snap/raw/v1"],
    )
    manifest = StudyManifest(
        study_id="snap",
        biological_design=BiologicalDesign(species="human", tissue="brain"),
        data={"raw_artifact_uri": "json://snap/raw/v1"},
    )
    result = TaskResult(
        task_id="task_result",
        status=task_status,
        capability="compute",
        method_used="fixture",
        output_artifacts=["json://snap/result/v1"],
    )
    report = ValidationReport(
        auditor_name="fixture_auditor",
        target_task_id="task_result",
        overall_passed=audit_passed,
        stop_rule_triggered=not audit_passed,
    )
    computed = result.model_copy(update={'status': TaskStatus.SUCCESS})
    registry.commit_task(transaction, 'fixture_signature', computed)
    transaction.close()
    contract = TaskContract(task_id='task_result', capability='compute', method='fixture',
                            expected_outputs=result.output_artifacts)
    registry.record_audit('fixture_signature', contract, computed, report, auditor_name='fixture_auditor')
    graph = EvidenceGraph()
    graph.add_evidence(
        EvidenceNode(
            evidence_id="E1",
            type=EvidenceType.QC_METRICS,
            summary="Observed quality metrics",
            source_task_id="task_result",
            source_artifact_uris=["json://snap/result/v1"],
            audit_passed=audit_passed,
        )
    )
    graph.add_claim(
        ClaimNode(
            claim_id="C1",
            statement="Observed quality metrics",
            language_tier=LanguageTier.LEVEL_1_OBSERVATION,
            claim_type=ClaimType.DESCRIPTIVE,
            support_evidence_ids=["E1"],
            confidence=ConfidenceScore(overall=0.5),
        )
    )
    return registry, manifest, result, report, graph


def _write(tmp_path: Path, **kwargs):
    registry, manifest, result, report, graph = _fixture(tmp_path)
    target = tmp_path / "study.snapshot.json"
    write_study_snapshot(
        target,
        manifest=manifest,
        config={"mode": "fixture"},
        summary={"status": "success"},
        evidence_graph=graph,
        artifact_registry=registry,
        task_history=[result],
        audit_reports=[report],
        **kwargs,
    )
    return target, registry, manifest, result, report, graph


def test_snapshot_roundtrip_rebuilds_graph_and_report_without_execution(tmp_path):
    target, registry, _, _, _, _ = _write(tmp_path)
    snapshot = load_study_snapshot(target)
    assert set(snapshot.evidence_graph.evidence_nodes) == {"E1"}
    assert set(snapshot.evidence_graph.claim_nodes) == {"C1"}
    assert snapshot.artifact_sources["json://snap/result/v1"].ancestor_hashes
    first = render_snapshot_report(target, registry)
    second = render_snapshot_report(target, registry)
    assert first == second
    assert "Scientific Study Report" in first
    assert "Observed quality metrics" in first


@pytest.mark.parametrize("numeric_value,covariates,invalid_field", [
    (0.2, [], None),
    (float("nan"), [], "lfc_se"),
    (float("inf"), [], "lfc_se"),
    (0.2, ["age"], "covariates"),
])
@pytest.mark.parametrize("storage_kind", ["legacy_csv", "table_v1"])
def test_audited_csv_empty_descriptors_do_not_hide_invalid_statistics(tmp_path, monkeypatch, numeric_value, covariates, invalid_field, storage_kind):
    import pandas as pd
    from eacbp.artifact.storage import PayloadSerializer, _sha256_file
    from eacbp.evidence.extraction import EvidenceRowBindingError, extract_evidence

    if storage_kind == "legacy_csv":
        # Construct historical-format bytes explicitly. Production only writes
        # Table v1 now, while this regression still exercises the old decoder.
        original_format = PayloadSerializer.preferred_format
        original_serialize = PayloadSerializer.serialize
        def legacy_format(cls, payload, artifact_type):
            return "native" if artifact_type == ArtifactType.TABLE else original_format(payload, artifact_type)
        def legacy_serialize(cls, payload, artifact_type, target_path, storage_format=None):
            if artifact_type == ArtifactType.TABLE:
                payload.to_csv(target_path, index=True)
                return _sha256_file(target_path)
            return original_serialize(payload, artifact_type, target_path, storage_format)
        monkeypatch.setattr(PayloadSerializer, "preferred_format", classmethod(legacy_format))
        monkeypatch.setattr(PayloadSerializer, "serialize", classmethod(legacy_serialize))

    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri, output_uri = "json://snap/raw/v1", "table://snap/deg/v1"
    registry.register(input_uri, {"raw": True}, ArtifactType.JSON, "snap", "ingest", "ingest")
    transaction = TaskArtifactTransaction(registry)
    transaction.register(output_uri, pd.DataFrame([{
        "gene": "fixture_gene", "log2_fold_change": 1.0, "fdr_q_value": 0.001,
        "lfc_se": numeric_value, "covariates": "", "status": "estimated", "status_reason": "",
    }]), ArtifactType.TABLE, "snap", "deg", "fixture", parent_uris=[input_uri],
        parameters={"condition_a": "stim", "condition_b": "ctrl"})
    contract = TaskContract(task_id="deg", capability="deg", method="fixture",
                            input_artifacts=[input_uri], expected_outputs=[output_uri])
    result = TaskResult(task_id="deg", capability="deg", method_used="fixture", status=TaskStatus.SUCCESS,
                        input_artifacts=[input_uri], output_artifacts=[output_uri],
                        metrics={"is_pseudobulk": True, "covariates": covariates})
    registry.commit_task(transaction, "descriptor_fixture", result)
    transaction.close()
    report = ValidationReport(auditor_name="fixture", target_task_id="deg", overall_passed=True)
    registry.record_audit("descriptor_fixture", contract, result, report, auditor_name="fixture")
    decoded = registry.load_payload(output_uri)
    if storage_kind == "legacy_csv":
        assert pd.isna(decoded.loc[0, "covariates"])
        assert pd.isna(decoded.loc[0, "status_reason"])
    else:
        assert decoded.loc[0, "covariates"] == "" and decoded.loc[0, "status_reason"] == ""
    if invalid_field:
        # New row-bound evidence rejects these already at construction; old
        # selector-free snapshots must still reject the same non-finite data.
        reason = "empty.*covariates" if storage_kind == "table_v1" and invalid_field == "covariates" else "non-finite.*" + invalid_field
        with pytest.raises(EvidenceRowBindingError, match=reason):
            extract_evidence(contract, result, report, registry)
        nodes = [EvidenceNode(evidence_id="invalid_row", type=EvidenceType.PSEUDOBULK_DEG,
                             source_task_id="deg", source_artifact_uris=[output_uri], audit_passed=True,
                             summary="Invalid legacy row", metrics={invalid_field: float("nan")})]
    else:
        nodes = extract_evidence(contract, result, report, registry)
    assert len(nodes) == 1
    if not invalid_field:
        assert nodes[0].metrics["status_reason"] == ""
        assert nodes[0].metrics["covariates"] == []
    graph = EvidenceGraph()
    graph.add_evidence(nodes[0])
    manifest = StudyManifest(study_id="snap", biological_design=BiologicalDesign(species="human", tissue="test"),
                             data={"raw_artifact_uri": input_uri})
    target = tmp_path / "snapshot.json"
    arguments = dict(manifest=manifest, config={}, summary={"status": "success"}, evidence_graph=graph,
                     artifact_registry=registry, task_history=[result], audit_reports=[report])
    if invalid_field:
        with pytest.raises(SnapshotValidationError, match="non-finite.*" + invalid_field):
            write_study_snapshot(target, **arguments)
    else:
        write_study_snapshot(target, **arguments)
        restored = load_study_snapshot(target)
        evidence = next(iter(restored.evidence_graph.evidence_nodes.values()))
        assert evidence.metrics["covariates"] == [] and evidence.metrics["status_reason"] == ""
        assert evidence.metrics["lfc_se"] == 0.2


@pytest.mark.parametrize("target_node", ["evidence", "claim"])
def test_snapshot_rejects_inference_identity_not_owned_by_source(tmp_path, target_node):
    registry, manifest, result, report, graph = _fixture(tmp_path)
    if target_node == "evidence":
        graph.evidence_nodes["E1"].inference_contract_id = "foreign"
    else:
        graph.claim_nodes["C1"].inference_contract_ids = ["foreign"]
    with pytest.raises(SnapshotValidationError, match="inference contract"):
        write_study_snapshot(tmp_path / "invalid.snapshot.json", manifest=manifest,
            config={}, summary={"status": "success"}, evidence_graph=graph,
            artifact_registry=registry, task_history=[result], audit_reports=[report])


def test_snapshot_rejects_payload_tampering_and_missing_files(tmp_path):
    target, registry, _, _, _, _ = _write(tmp_path)
    metadata = registry.get_metadata("json://snap/result/v1")
    Path(metadata.storage_path).write_text('{"value": 9}', encoding="utf-8")
    with pytest.raises(SnapshotIntegrityError, match="changed"):
        render_snapshot_report(target, registry)

    # Restore a fresh fixture and remove the source payload altogether.
    target, registry, _, _, _, _ = _write(tmp_path / "missing")
    metadata = registry.get_metadata("json://snap/raw/v1")
    Path(metadata.storage_path).unlink()
    with pytest.raises(SnapshotIntegrityError, match="missing"):
        render_snapshot_report(target, registry)


def test_snapshot_rejects_schema_version_and_dangling_evidence(tmp_path):
    target, registry, _, _, _, _ = _write(tmp_path)
    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["schema_version"] = 999
    target.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SnapshotSchemaError, match="schema version"):
        load_study_snapshot(target)

    target, registry, _, _, _, _ = _write(tmp_path / "dangling")
    payload = json.loads(target.read_text(encoding="utf-8"))
    payload["evidence"]["claims"][0]["support_evidence_ids"] = ["missing"]
    from eacbp.evidence.snapshot import _stable_fingerprint
    payload['fingerprints']['evidence'] = _stable_fingerprint(payload['evidence'])
    target.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SnapshotValidationError, match="missing evidence"):
        load_study_snapshot(target)


def test_failed_audit_or_scientific_failure_cannot_be_admitted(tmp_path):
    registry, manifest, result, report, graph = _fixture(tmp_path / "failed", audit_passed=False)
    with pytest.raises(SnapshotValidationError, match="passing audit|not admitted"):
        write_study_snapshot(
            tmp_path / "failed.snapshot.json",
            manifest=manifest,
            config={},
            summary={},
            evidence_graph=graph,
            artifact_registry=registry,
            task_history=[result],
            audit_reports=[report],
        )

    registry, manifest, result, report, graph = _fixture(tmp_path / "scientific", task_status=TaskStatus.SCIENTIFIC_FAILURE)
    with pytest.raises(SnapshotValidationError, match="failed task"):
        write_study_snapshot(
            tmp_path / "scientific.snapshot.json",
            manifest=manifest,
            config={},
            summary={},
            evidence_graph=graph,
            artifact_registry=registry,
            task_history=[result],
            audit_reports=[report],
        )


def test_raw_import_without_audit_is_only_an_ancestor(tmp_path):
    target, registry, _, _, _, _ = _write(tmp_path)
    snapshot = load_study_snapshot(target)
    assert "json://snap/raw/v1" in snapshot.artifact_sources["json://snap/result/v1"].ancestor_uris
    assert all(context.task_id != "task_import" for context in snapshot.audit_contexts.values())


def test_atomic_write_failure_keeps_previous_snapshot(tmp_path, monkeypatch):
    target, registry, manifest, result, report, graph = _write(tmp_path)
    original = target.read_bytes()

    import eacbp.evidence.snapshot as snapshot_module

    def fail_replace(*args, **kwargs):
        raise OSError("simulated replace failure")

    monkeypatch.setattr(snapshot_module.os, "replace", fail_replace)
    with pytest.raises(SnapshotWriteError):
        write_study_snapshot(
            target,
            manifest=manifest,
            config={"changed": True},
            summary={},
            evidence_graph=graph,
            artifact_registry=registry,
            task_history=[result],
            audit_reports=[report],
        )
    assert target.read_bytes() == original


def test_snapshot_rejects_revoked_audit_and_missing_commit(tmp_path):
    target, registry, _, result, report, _ = _write(tmp_path / 'revoke')
    contract = registry.get_audit_record('fixture_signature').contract
    registry.begin_audit('fixture_signature', contract, result, auditor_name='fixture_auditor')
    with pytest.raises(SnapshotIntegrityError, match='audit receipt changed'):
        render_snapshot_report(target, registry)
    target, registry, _, _, _, _ = _write(tmp_path / 'missing_commit')
    index = Path(registry.storage.base_dir) / registry.INDEX_FILENAME
    raw = json.loads(index.read_text())
    raw['task_commits'] = {}
    index.write_text(json.dumps(raw))
    with pytest.raises(SnapshotIntegrityError, match='computation receipt'):
        render_snapshot_report(target, registry)


def test_snapshot_detects_content_edits_and_new_artifacts(tmp_path):
    target, registry, _, _, _, _ = _write(tmp_path)
    raw = json.loads(target.read_text())
    raw['evidence']['claims'][0]['statement'] = 'Changed conclusion'
    target.write_text(json.dumps(raw))
    with pytest.raises(SnapshotIntegrityError, match='fingerprint'):
        load_study_snapshot(target)
    target, registry, _, _, _, _ = _write(tmp_path / 'extra')
    registry.register('json://snap/later/v1', {}, ArtifactType.JSON, 'snap', 'later', 'import')
    with pytest.raises(SnapshotIntegrityError, match='scope changed'):
        render_snapshot_report(target, registry)


def test_snapshot_rejects_report_only_or_wrong_task_source(tmp_path):
    registry, manifest, result, report, graph = _fixture(tmp_path)
    node = graph.evidence_nodes['E1']
    node.source_artifact_uris = ['json://snap/raw/v1']
    with pytest.raises(SnapshotValidationError, match='not a task output'):
        write_study_snapshot(tmp_path / 'wrong.json', manifest=manifest, config={}, summary={},
                             evidence_graph=graph, artifact_registry=registry,
                             task_history=[result], audit_reports=[report])
    node.source_artifact_uris = result.output_artifacts
    index = Path(registry.storage.base_dir) / registry.INDEX_FILENAME
    raw = json.loads(index.read_text())
    raw['audit_records'] = {}
    index.write_text(json.dumps(raw))
    with pytest.raises(SnapshotValidationError, match='passing audit'):
        write_study_snapshot(tmp_path / 'unadmitted.json', manifest=manifest, config={}, summary={},
                             evidence_graph=graph, artifact_registry=registry,
                             task_history=[result], audit_reports=[report])
