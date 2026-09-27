"""Control-plane reads stay bounded; relocation preserves scientific gates."""
import json
from pathlib import Path
import shutil
import zipfile

import pytest

from eacbp.artifact.portability import (
    ArtifactBundleError, export_artifact_bundle, import_artifact_bundle,
)
from eacbp.artifact.registry import ArtifactRegistry, ArtifactAuditAccessError
from eacbp.artifact.storage import ArtifactIntegrityError, PayloadSerializer
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.snapshot import write_study_snapshot, render_snapshot_report, load_study_snapshot
from eacbp.orchestrator.resume import ResumeManager
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def published(root):
    registry = ArtifactRegistry(str(root / "artifacts"))
    raw = registry.register("json://s/raw/v1", {"raw": 1}, ArtifactType.JSON, "s", "raw", "import")
    with TaskArtifactTransaction(registry) as transaction:
        for name in ("a", "b"):
            transaction.register(f"json://s/{name}/v1", {"value": name}, ArtifactType.JSON,
                                 "s", "task", "compute", parent_uris=[raw.uri])
        result = TaskResult(task_id="task", capability="compute", method_used="fixture",
                            status=TaskStatus.SUCCESS,
                            input_artifacts=[raw.uri],
                            output_artifacts=["json://s/a/v1", "json://s/b/v1"])
        registry.commit_task(transaction, "signature", result)
    contract = TaskContract(task_id="task", capability="compute", input_artifacts=[raw.uri],
                            expected_outputs=result.output_artifacts)
    report = ValidationReport(auditor_name="fixture", target_task_id="task", overall_passed=True)
    registry.record_audit("signature", contract, result, report)
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="brain"),
                             data={"raw_artifact_uri": raw.uri})
    graph = EvidenceGraph()
    graph.add_evidence(EvidenceNode(evidence_id="E1", type=EvidenceType.QC_METRICS,
                                    summary="Quality", source_task_id="task",
                                    source_artifact_uris=result.output_artifacts, audit_passed=True))
    snapshot = root / "snapshot.json"
    write_study_snapshot(snapshot, manifest=manifest, config={}, summary={}, evidence_graph=graph,
                         artifact_registry=registry, task_history=[result], audit_reports=[report])
    return registry, contract, result, snapshot


def test_verification_and_resume_never_deserialize(tmp_path, monkeypatch):
    registry, contract, result, _ = published(tmp_path)
    def forbidden(*args, **kwargs):
        pytest.fail("control-plane verification deserialized a payload")
    monkeypatch.setattr(PayloadSerializer, "deserialize", forbidden)
    refreshes = []
    original = registry._refresh
    def refresh():
        refreshes.append(True)
        original()
    monkeypatch.setattr(registry, "_refresh", refresh)
    assert set(registry.verify_many(result.output_artifacts)) == set(result.output_artifacts)
    assert len(refreshes) == 1
    manager = ResumeManager(registry)
    assert manager.input_hashes(contract)
    journal = type("Journal", (), {"entries": {}})()
    assert manager.lookup(journal, contract, "signature").reused
    Path(registry.get_metadata(result.output_artifacts[0]).storage_path).write_text("tampered")
    with pytest.raises(ArtifactIntegrityError):
        manager.lookup(journal, contract, "signature")


def test_audited_queries_deserialize_each_sibling_once(tmp_path, monkeypatch):
    registry, contract, result, _ = published(tmp_path)
    loaded = []
    original = PayloadSerializer.deserialize
    def deserialize(path, artifact_type):
        loaded.append(path)
        return original(path, artifact_type)
    monkeypatch.setattr(PayloadSerializer, "deserialize", deserialize)
    assert len(registry.get_audited_many("signature", contract)) == 2
    assert len(loaded) == 2
    loaded.clear()
    assert registry.get_audited(result.output_artifacts[0], "signature", contract)[1] == {"value": "a"}
    assert len(loaded) == 2


def test_directory_copy_and_bundle_preserve_report_and_receipts(tmp_path):
    source = tmp_path / "source"
    registry, contract, result, snapshot = published(source)
    before = render_snapshot_report(snapshot, registry)
    index = json.loads((source / "artifacts" / registry.INDEX_FILENAME).read_text())
    assert index["index_version"] == 2
    assert all(not Path(item["storage_path"]).is_absolute() for item in index["metadata"])
    moved = tmp_path / "moved"
    shutil.copytree(source, moved)
    relocated = ArtifactRegistry(str(moved / "artifacts"))
    assert render_snapshot_report(moved / "snapshot.json", relocated) == before
    bundle = export_artifact_bundle(source / "artifacts", tmp_path / "study.zip",
                                    snapshot_paths={"snapshot.json": snapshot})
    imported = import_artifact_bundle(bundle, tmp_path / "imported")
    restored = ArtifactRegistry(str(imported / "artifacts"))
    assert render_snapshot_report(imported / "snapshot.json", restored) == before
    assert restored.get_task_commit("signature") == registry.get_task_commit("signature")
    assert restored.get_audit_record("signature") == registry.get_audit_record("signature")
    assert restored.get_audited(result.output_artifacts[0], "signature", contract)[1] == {"value": "a"}
    restored.begin_audit("signature", contract, result)
    with pytest.raises(ArtifactAuditAccessError):
        restored.get_audited(result.output_artifacts[0], "signature", contract)
    with pytest.raises(FileExistsError):
        import_artifact_bundle(bundle, imported)


@pytest.mark.parametrize("mode", ["tamper", "traversal", "absolute", "duplicate"])
def test_invalid_bundle_never_publishes_destination(tmp_path, mode):
    source = tmp_path / "source"
    registry, _, _, snapshot = published(source)
    bundle = export_artifact_bundle(registry.storage.base_dir, tmp_path / "good.zip",
                                    snapshot_paths={"snapshot.json": snapshot})
    malicious = tmp_path / "bad.zip"
    with zipfile.ZipFile(bundle) as original, zipfile.ZipFile(malicious, "w") as changed:
        for key in original.namelist():
            content = original.read(key)
            if mode == "tamper" and key == "snapshot.json":
                content += b" "
            changed.writestr(key, content)
        if mode == "traversal":
            changed.writestr("../escape", "bad")
        elif mode == "absolute":
            changed.writestr("C:/escape", "bad")
        elif mode == "duplicate":
            changed.writestr("SNAPSHOT.JSON", "bad")
    destination = tmp_path / "destination"
    with pytest.raises(ArtifactBundleError):
        import_artifact_bundle(malicious, destination)
    assert not destination.exists()
    assert not (tmp_path / "escape").exists()


def test_snapshot_v2_lineage_is_derived_and_checked(tmp_path):
    registry, _, _, snapshot_path = published(tmp_path)
    snapshot = load_study_snapshot(snapshot_path)
    assert snapshot.evidence_graph.evidence_nodes["E1"].data_origin_uris == ["json://s/raw/v1"]
    # Even a caller that recomputes a section fingerprint cannot forge roots.
    from eacbp.evidence.snapshot import _stable_fingerprint, SnapshotValidationError
    payload = json.loads(snapshot_path.read_text())
    payload["evidence"]["nodes"][0]["data_origin_uris"] = ["json://s/a/v1"]
    payload["fingerprints"]["evidence"] = _stable_fingerprint(payload["evidence"])
    snapshot_path.write_text(json.dumps(payload))
    with pytest.raises(SnapshotValidationError, match="provenance"):
        load_study_snapshot(snapshot_path)


def test_legacy_absolute_registry_and_snapshot_export_to_portable_v2(tmp_path):
    from eacbp.evidence.snapshot import _stable_fingerprint
    source = tmp_path / "source"
    registry, _, _, snapshot = published(source)
    metadata = {item.uri: item for item in registry.list_artifacts()}
    index_path = registry.storage.base_dir / registry.INDEX_FILENAME
    legacy_index = json.loads(index_path.read_text())
    legacy_index["index_version"] = 1
    for item in legacy_index["metadata"]:
        item["storage_path"] = metadata[item["uri"]].storage_path
    index_path.write_text(json.dumps(legacy_index))
    legacy_snapshot = json.loads(snapshot.read_text())
    legacy_snapshot["schema_version"] = 1
    for section in ("metadata", "sources"):
        for item in legacy_snapshot["artifacts"][section]:
            item["storage_path"] = metadata[item["uri"]].storage_path
    legacy_snapshot["fingerprints"]["artifacts"] = _stable_fingerprint(legacy_snapshot["artifacts"])
    snapshot.write_text(json.dumps(legacy_snapshot))
    bundle = export_artifact_bundle(registry.storage.base_dir, tmp_path / "legacy.zip",
                                    snapshot_paths={"snapshot.json": snapshot})
    destination = import_artifact_bundle(bundle, tmp_path / "relocated")
    restored = ArtifactRegistry(str(destination / "artifacts"))
    assert load_study_snapshot(destination / "snapshot.json").schema_version == 2
    assert render_snapshot_report(destination / "snapshot.json", restored)
    # Export did not silently migrate the original archive in place.
    assert json.loads(index_path.read_text())["index_version"] == 1


def test_export_rejects_active_study_and_tampered_source(tmp_path):
    from eacbp.orchestrator.checkpoint import StudyJournal
    registry, _, result, snapshot = published(tmp_path / "source")
    bundle = tmp_path / "study.zip"
    with StudyJournal(registry.storage.base_dir, "s"):
        with pytest.raises(ArtifactBundleError, match="running"):
            export_artifact_bundle(registry.storage.base_dir, bundle,
                                   snapshot_paths={"snapshot.json": snapshot})
    assert not bundle.exists()
    Path(registry.get_metadata(result.output_artifacts[0]).storage_path).write_text("tampered")
    with pytest.raises(ArtifactIntegrityError):
        export_artifact_bundle(registry.storage.base_dir, bundle,
                               snapshot_paths={"snapshot.json": snapshot})
    assert not bundle.exists()


def test_artifact_only_bundle_is_generic_but_not_an_application_study(tmp_path):
    registry, _, _, _ = published(tmp_path / "source")
    bundle = export_artifact_bundle(registry.storage.base_dir, tmp_path / "artifacts.zip")
    assert import_artifact_bundle(bundle, tmp_path / "generic").is_dir()
    with pytest.raises(ArtifactBundleError, match="requires"):
        import_artifact_bundle(bundle, tmp_path / "study", require_study=True)
    assert not (tmp_path / "study").exists()


@pytest.mark.parametrize("invalid", [None, "incomplete", "version", "study_id", "manifest", "no_config"])
def test_study_import_requires_a_matching_completed_envelope(tmp_path, invalid):
    registry, _, _, snapshot = published(tmp_path / "source")
    model = load_study_snapshot(snapshot).manifest
    envelope = {"schema_version": 1, "import_completed": True, "study_id": model.study_id,
                "manifest": model.model_dump(mode="json"), "config": {}}
    if invalid == "incomplete":
        envelope["import_completed"] = False
    elif invalid == "version":
        envelope["schema_version"] = True
    elif invalid == "study_id":
        envelope["study_id"] = "unrelated"
    elif invalid == "manifest":
        envelope["manifest"]["biological_design"]["tissue"] = "heart"
    elif invalid == "no_config":
        envelope.pop("config")
    config_path = tmp_path / "source" / "run_config.json"
    config_path.write_text(json.dumps(envelope))
    bundle = export_artifact_bundle(registry.storage.base_dir, tmp_path / "study.zip",
                                    snapshot_paths={"snapshot.json": snapshot},
                                    run_files={"run_config.json": config_path})
    destination = tmp_path / "imported"
    if invalid is None:
        assert import_artifact_bundle(bundle, destination, require_study=True) == destination
        assert export_artifact_bundle(
            registry.storage.base_dir, tmp_path / "strict-study.zip", require_study=True,
            snapshot_paths={"snapshot.json": snapshot}, run_files={"run_config.json": config_path},
        ).is_file()
    else:
        with pytest.raises(ArtifactBundleError):
            import_artifact_bundle(bundle, destination, require_study=True)
        assert not destination.exists()
        strict_bundle = tmp_path / "strict-study.zip"
        with pytest.raises(ArtifactBundleError):
            export_artifact_bundle(
                registry.storage.base_dir, strict_bundle, require_study=True,
                snapshot_paths={"snapshot.json": snapshot}, run_files={"run_config.json": config_path},
            )
        assert not strict_bundle.exists()
