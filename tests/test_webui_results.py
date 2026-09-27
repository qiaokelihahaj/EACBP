"""Browser projections retain the saved evidence and storage boundaries."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor.base import ValidationReport
from eacbp.capabilities.sc_data import SCData
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.snapshot import write_study_snapshot
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.webui.jobs import write_json
from eacbp.webui.results import result_views
from eacbp.webui.service import WebService


def fixture_run(tmp_path, *, admitted=True, method="kmeans_marker_embedding_v1", simulated=False,
                raw_padding=0, deg_rows=3):
    run = tmp_path / "runs" / "fixture"
    registry = ArtifactRegistry(str(run / "artifacts"))
    registry.register("json://fixture/raw/v1", {"raw": True, "padding": "x" * raw_padding}, ArtifactType.JSON,
                      "fixture", "import", "import", summary_metrics={"is_simulated": simulated})
    data = SCData(np.ones((5, 2)), pd.DataFrame({"cluster": ["0", "0", "1", "1", "1"],
                  "cell_type": ["A", "A", "B", "B", "B"]}), pd.DataFrame(index=["A", "B"]),
                  obsm={"X_embedding_2d": np.arange(10).reshape(5, 2),
                        "X_umap": np.arange(10).reshape(5, 2)})
    if deg_rows == 3:
        table = pd.DataFrame({"gene": ["gene_a", "gene_b", "not_tested"],
                              "log2_fold_change": [1., -.5, np.nan],
                              "fdr_q_value": [.04, .2, np.nan], "p_value": [.01, .1, np.nan]})
    else:
        genes = [f"gene_{index:04d}" for index in range(deg_rows - 3)] + [
            "target_after_200", "001", "NA"]
        table = pd.DataFrame({"gene": genes,
                              "log2_fold_change": np.linspace(-2, 2, deg_rows),
                              "fdr_q_value": np.linspace(.001, .99, deg_rows),
                              "p_value": np.linspace(.0001, .999, deg_rows)})
    tasks, reports = [], []
    for capability, payload, kind in (("qc", data, ArtifactType.ANNDATA),
                                      ("clustering", data, ArtifactType.ANNDATA),
                                      ("deg", table, ArtifactType.TABLE)):
        task_id = "task_" + capability
        uri = f"{'table' if kind == ArtifactType.TABLE else 'adata'}://fixture/{capability}/v1"
        transaction = TaskArtifactTransaction(registry)
        transaction.register(uri, payload, kind, "fixture", task_id, capability,
                             parent_uris=["json://fixture/raw/v1"])
        task = TaskResult(task_id=task_id, status=TaskStatus.SUCCESS,
                          capability=capability, method_used=method if capability == "clustering" else "fixture",
                          output_artifacts=[uri], metrics={"initial_cells": 8,
                            "condition_a": "treated", "condition_b": "control",
                            "statistical_unit": "donor_pseudobulk"})
        contract = TaskContract(task_id=task_id, capability=capability, expected_outputs=[uri])
        registry.commit_task(transaction, "signature_" + capability, task)
        transaction.close()
        report = ValidationReport(auditor_name="fixture", target_task_id=task_id, overall_passed=admitted)
        registry.record_audit("signature_" + capability, contract, task, report, auditor_name="fixture")
        if not admitted:
            task.status = TaskStatus.SCIENTIFIC_FAILURE
        tasks.append(task)
        reports.append(report)
    manifest = StudyManifest(study_id="fixture", biological_design={"species": "human", "tissue": "fixture"},
                             data={"raw_artifact_uri": "json://fixture/raw/v1"})
    write_study_snapshot(run / "snapshot.json", manifest=manifest, config={"mode": "real"},
                         summary={}, evidence_graph=EvidenceGraph(), artifact_registry=registry,
                         task_history=tasks, audit_reports=reports)
    return run, registry


def test_results_show_admitted_qc_embedding_and_deg_without_mutation(tmp_path):
    run, registry = fixture_run(tmp_path, simulated=True)
    before = {path: path.read_bytes() for path in run.rglob("*") if path.is_file()}
    response = result_views(run)
    assert response["state"] == "ready", response
    assert response["qc"]["initial_cells"] == 8
    assert response["qc"]["retained_cells"] == 5
    assert response["qc"]["filtered_cells"] == 3
    assert response["embedding"]["label"] == "二维展示投影（非 UMAP）"
    assert response["embedding"]["groups"] == [{"label": "0", "count": 2}, {"label": "1", "count": 3}]
    assert response["deg"][0]["rows"][2]["fdr_q_value"] is None
    assert response["deg"][0]["condition_a"] == "treated"
    assert any("模拟" in item for item in response["warnings"])
    assert "未重新执行" in response["notice"]
    json.dumps(response, allow_nan=False)
    assert all(path.read_bytes() == content for path, content in before.items())


def test_only_real_umap_method_uses_umap_label_and_sample_is_bounded(tmp_path, monkeypatch):
    run, _ = fixture_run(tmp_path, method="scanpy_leiden_umap_v1")
    monkeypatch.setattr("eacbp.webui.results.MAX_POINTS", 2)
    response = result_views(run)
    assert response["embedding"]["label"] == "UMAP"
    assert response["embedding"]["sampled"] is True
    assert response["embedding"]["shown_cells"] == 2
    assert response["embedding"]["total_cells"] == 5


@pytest.mark.parametrize("tamper", ["payload", "ancestor", "snapshot", "receipt"])
def test_results_fail_closed_after_tampering(tmp_path, tamper):
    run, registry = fixture_run(tmp_path)
    if tamper in {"payload", "ancestor"}:
        uri = "table://fixture/deg/v1" if tamper == "payload" else "json://fixture/raw/v1"
        Path(registry.get_metadata(uri).storage_path).write_bytes(b"corrupted")
    elif tamper == "snapshot":
        path = run / "snapshot.json"
        data = json.loads(path.read_text())
        data["config"]["mode"] = "demo"
        path.write_text(json.dumps(data))
    else:
        path = run / "artifacts" / ArtifactRegistry.INDEX_FILENAME
        data = json.loads(path.read_text())
        data["audit_records"]["signature_qc"]["status"] = "rejected"
        path.write_text(json.dumps(data))
    response = result_views(run)
    assert response["state"] == "unavailable"
    assert response["qc"] is None and response["embedding"] is None and response["deg"] == []


def test_failed_audits_are_never_charts(tmp_path):
    run, _ = fixture_run(tmp_path, admitted=False)
    response = result_views(run)
    assert response["state"] == "unavailable"
    assert response["qc"] is None and not response["deg"]
    assert response["warnings"]


def test_snapshot_pass_without_live_durable_receipt_is_not_admitted(tmp_path):
    run, _ = fixture_run(tmp_path)
    path = run / "artifacts" / ArtifactRegistry.INDEX_FILENAME
    data = json.loads(path.read_text())
    data.pop("audit_records")
    path.write_text(json.dumps(data))
    response = result_views(run)
    assert response["state"] == "unavailable"
    assert response["qc"] is None and response["embedding"] is None and response["deg"] == []
    assert "audit receipt is missing" in response["notice"]


def test_expanded_h5ad_limit_is_checked_before_payload_load(tmp_path):
    h5py = pytest.importorskip("h5py")
    from eacbp.webui.results import _check_expanded_size
    path = tmp_path / "compressed.h5ad"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", shape=(10_000, 10_000), dtype="f4", chunks=True, compression="gzip")
    assert path.stat().st_size < 32_000_000
    with pytest.raises(ValueError, match="大小限制"):
        _check_expanded_size(path)


def test_h5ad_external_links_are_not_followed(tmp_path):
    h5py = pytest.importorskip("h5py")
    from eacbp.webui.results import _check_expanded_size
    path = tmp_path / "external.h5ad"
    with h5py.File(path, "w") as handle:
        handle["obs"] = h5py.ExternalLink(str(tmp_path / "outside.h5ad"), "/obs")
    with pytest.raises(ValueError, match="外部"):
        _check_expanded_size(path)


def test_size_limits_and_paths_and_active_jobs(tmp_path, monkeypatch):
    service = WebService(tmp_path, tmp_path / "runs")
    with pytest.raises(ValueError):
        service.results("../outside")
    run, registry = fixture_run(tmp_path)
    monkeypatch.setattr("eacbp.webui.results.MAX_FILE_BYTES", 5)
    response = service.results("fixture")
    assert response["state"] == "unavailable" and any("大小限制" in item for item in response["warnings"])
    write_json(service.manager.jobs_dir / "queued" / "job.json", {
        "id": "queued", "status": "queued", "created_at": "", "created_epoch": 9999999999,
        "run_dir": str(run)})
    assert "后台操作" in service.results("fixture")["notice"]


def test_large_ancestor_is_hash_checked_without_hiding_small_results(tmp_path, monkeypatch):
    run, registry = fixture_run(tmp_path, raw_padding=100_000)
    monkeypatch.setattr("eacbp.webui.results.MAX_FILE_BYTES", 64_000)
    raw = Path(registry.get_metadata("json://fixture/raw/v1").storage_path)
    assert raw.stat().st_size > 64_000
    response = result_views(run)
    assert response["state"] == "ready", response
    assert response["deg"] and response["qc"] and response["embedding"]
    # Skipping ancestor deserialization must not skip its integrity gate.
    raw.write_bytes(b"corrupt")
    assert result_views(run)["state"] == "unavailable"


def test_oversized_chart_outputs_are_skipped_but_small_deg_remains(tmp_path, monkeypatch):
    run, registry = fixture_run(tmp_path)
    deg_size = Path(registry.get_metadata("table://fixture/deg/v1").storage_path).stat().st_size
    monkeypatch.setattr("eacbp.webui.results.MAX_FILE_BYTES", deg_size + 1)
    response = result_views(run)
    assert response["state"] == "ready", response
    assert response["qc"] is None and response["embedding"] is None
    assert response["deg"][0]["total_genes"] == 3
    assert any("task_qc" in item and "大小限制" in item for item in response["warnings"])
    assert any("task_clustering" in item and "大小限制" in item for item in response["warnings"])


def test_total_streaming_hash_budget_still_fails_closed(tmp_path, monkeypatch):
    run, _ = fixture_run(tmp_path)
    monkeypatch.setattr("eacbp.webui.results.MAX_HASH_BYTES", 1)
    response = result_views(run)
    assert response["state"] == "unavailable" and not response["deg"]
    assert "读取预算" in response["notice"]


def test_registry_cannot_reference_file_outside_run(tmp_path):
    run, registry = fixture_run(tmp_path)
    path = run / "artifacts" / ArtifactRegistry.INDEX_FILENAME
    data = json.loads(path.read_text())
    data["index_version"] = 1
    data["metadata"][0]["storage_path"] = str(tmp_path / "outside.csv")
    path.write_text(json.dumps(data))
    response = result_views(run)
    assert response["state"] == "unavailable"
    assert response["qc"] is None


def test_audited_deg_table_searches_and_pages_beyond_preview_limit(tmp_path):
    fixture_run(tmp_path, deg_rows=260)
    service = WebService(tmp_path, tmp_path / "runs")

    first = service.deg_table("fixture", task_id="task_deg", page=1, size=50)
    assert first["total_items"] == 260
    assert first["total_pages"] == 6
    assert len(first["rows"]) == 50
    assert first["rows"][0]["gene"] == "gene_0000"

    last = service.deg_table("fixture", task_id="task_deg", page=6, size=50)
    assert len(last["rows"]) == 10
    assert last["rows"][-3]["gene"] == "target_after_200"
    assert [row["gene"] for row in last["rows"][-2:]] == ["001", "NA"]

    searched = service.deg_table("fixture", task_id="task_deg", page=1, size=20,
                                 query="target_after_200")
    assert searched["total_items"] == 1
    assert searched["rows"][0]["gene"] == "target_after_200"
    assert "001" in [row["gene"] for row in service.deg_table(
        "fixture", task_id="task_deg", query="001")["rows"]]
    assert service.deg_table("fixture", task_id="task_deg", query="NA")["rows"][0]["gene"] == "NA"

    filtered = service.deg_table("fixture", task_id="task_deg", page=1, size=100,
                                 fdr_max=.01, significant_only=True)
    assert filtered["total_items"] == 3
    assert all(row["fdr_q_value"] <= .01 for row in filtered["rows"])


@pytest.mark.parametrize("arguments", [
    {"page": 0}, {"page": -1}, {"page": 1.5}, {"size": 0}, {"size": 101},
    {"fdr_max": -0.01}, {"fdr_max": 1.01}, {"fdr_max": float("nan")},
    {"significant_only": True},
])
def test_deg_page_rejects_invalid_pagination_and_filters(tmp_path, arguments):
    fixture_run(tmp_path, deg_rows=220)
    service = WebService(tmp_path, tmp_path / "runs")
    with pytest.raises(ValueError):
        service.deg_table("fixture", task_id="task_deg", **arguments)


def test_full_deg_csv_export_uses_current_audit_and_rejects_tampering(tmp_path):
    run, registry = fixture_run(tmp_path, deg_rows=260)
    service = WebService(tmp_path, tmp_path / "runs")
    source, filename = service.deg_table_csv("fixture", task_id="task_deg")
    assert filename.endswith(".csv")
    assert "target_after_200" in source.read_text(encoding="utf-8")

    metadata = registry.get_metadata("table://fixture/deg/v1")
    Path(metadata.storage_path).write_text("gene,log2_fold_change,fdr_q_value\nforged,9,.001\n")
    with pytest.raises((ValueError, RuntimeError)):
        service.deg_table_csv("fixture", task_id="task_deg")


def test_large_preview_skipped_deg_remains_available_via_verified_catalog(tmp_path, monkeypatch):
    fixture_run(tmp_path, deg_rows=260)
    service = WebService(tmp_path, tmp_path / "runs")
    monkeypatch.setattr("eacbp.webui.results.MAX_FILE_BYTES", 1)

    preview = service.results("fixture")
    assert preview["deg"] == []
    catalog = service.deg_tables("fixture")
    assert catalog["tables"][0]["task_id"] == "task_deg"
    assert service.deg_table("fixture", task_id="task_deg")["total_items"] == 260
