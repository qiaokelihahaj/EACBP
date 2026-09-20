"""Web boundary, process lifecycle and original scientific-gate regression tests."""
import http.client
import json
from pathlib import Path
import threading
import time

import numpy as np
import pandas as pd
import pytest

from eacbp.webui.jobs import JobManager, lease, read_json, write_json
from eacbp.webui.server import WebServer
from eacbp.webui.service import WebService


@pytest.fixture
def service(tmp_path):
    return WebService(tmp_path, tmp_path / "runs")


def input_data(tmp_path):
    ad = pytest.importorskip("anndata")
    rng = np.random.default_rng(24)
    counts = rng.poisson(4, (80, 40)).astype(np.float32)
    data = ad.AnnData(counts, obs=pd.DataFrame(index=[f"cell{i}" for i in range(80)]),
                      var=pd.DataFrame(index=[f"g{i}" for i in range(40)]))
    path = tmp_path / "controlled_fixture.h5ad"
    data.write_h5ad(path)
    return path


def form(data):
    return {"data": str(data), "study_id": "ui_test", "species": "human", "tissue": "synthetic_test",
            "title": "Synthetic software validation only", "method_profile": "baseline", "min_genes": 1}


def wait_job(manager, job_id, timeout=90):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        job = next(j for j in manager.jobs() if j["id"] == job_id)
        if job["status"] not in {"running", "queued"}:
            return job
        time.sleep(.15)
    pytest.fail("Worker did not finish within the test timeout")


def test_form_plan_and_metadata_keep_input_read_only(service, tmp_path):
    data = input_data(tmp_path)
    before = data.read_bytes()
    info = service.dataset(data)
    assert info["n_cells"] == 80 and info["n_genes"] == 40
    result = service.preview({**form(data), "condition_a": "treated", "condition_b": "control",
                              "donor_col": "patient", "advanced_analysis": True,
                              "config_overrides": {"capability_parameters": {"deg": {"paired": True}}}})
    assert result["plan"]["phase"] == "before_dataset_audit"
    assert result["config"]["capability_parameters"]["deg"]["condition_a"] == "treated"
    assert result["config"]["capability_parameters"]["deg"]["paired"] is True
    assert result["config"]["capability_parameters"]["dataset_audit"]["donor_col"] == "patient"
    assert not list(service.manager.runs_dir.glob("*/run_config.json"))
    assert data.read_bytes() == before


@pytest.mark.parametrize("changes", [
    {"study_id": "../../escape"}, {"condition_a": "A"}, {"condition_a": "A", "condition_b": "A"},
    {"advanced_analysis": "false"}, {"paired": True}, {"config_overrides": {"mode": "demo"}},
    {"config_overrides": {"resume": True}}, {"max_mito_pct": 101},
    {"config_overrides": {"capability_parameters": {"deg": {"paired": True}}}},
    {"paired": True, "advanced_analysis": True,
     "config_overrides": {"method_overrides": {"deg": "donor_pseudobulk_welch_v1"}}},
])
def test_form_rejects_ambiguous_or_unsafe_controls(service, tmp_path, changes):
    data = tmp_path / "input.h5ad"
    data.touch()
    with pytest.raises((ValueError, TypeError)):
        service.prepare({**form(data), **changes})


def test_file_and_run_boundaries(service, tmp_path):
    with pytest.raises(ValueError):
        service.browse("..")
    with pytest.raises(ValueError):
        service.data_path(tmp_path.parent / "outside.h5ad")
    with pytest.raises(ValueError):
        service.manager.run_path("../escape")
    with pytest.raises(ValueError):
        service.manager.run_path(".webui/jobs")
    with pytest.raises(ValueError):
        service.manager.run_path("")


def test_optional_null_overrides_use_typed_config_defaults(service, tmp_path):
    data = tmp_path / "input.h5ad"
    data.touch()
    _, _, _, config = service.prepare({**form(data), "config_overrides": {
        "method_overrides": None, "capability_parameters": None,
    }})
    assert "method_overrides" not in config and "capability_parameters" not in config
    with pytest.raises(ValueError):
        service.prepare({**form(data), "config_overrides": {"method_overrides": False}})


def test_stale_worker_and_corrupt_records_block_or_recover(service):
    directory = service.manager.jobs_dir / "stale"
    write_json(directory / "job.json", {"id": "stale", "status": "running", "created_at": "",
                                         "created_epoch": 0, "run_dir": ""})
    with lease(directory / "worker.lock"):
        assert service.manager.jobs()[0]["status"] == "running"
    assert service.manager.jobs()[0]["status"] == "interrupted"
    (directory / "job.json").write_text("bad JSON")
    assert service.manager.jobs()[0]["status"] == "unknown"
    with pytest.raises(BlockingIOError):
        service.manager.submit("run", run_dir=service.manager.runs_dir / "new")


def test_spawn_failure_is_persisted(service, monkeypatch):
    def fail(*args, **kwargs):
        raise OSError("cannot spawn")
    monkeypatch.setattr("eacbp.webui.jobs.subprocess.Popen", fail)
    with pytest.raises(OSError, match="cannot spawn"):
        service.manager.submit("run", run_dir=service.manager.runs_dir / "new")
    job = service.manager.jobs()[0]
    assert job["status"] == "failed" and job["error"] == "cannot spawn"
    assert not (service.manager.runs_dir / "new").exists()


@pytest.fixture
def http_server(service):
    server = WebServer(service, port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield server
    server.shutdown()
    server.server_close()
    thread.join(timeout=5)


def request(server, path, method="GET", body=None, headers=None):
    connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=15)
    connection.request(method, path, body=body, headers=headers or {})
    response = connection.getresponse()
    result = (response.status, dict(response.getheaders()), response.read())
    connection.close()
    return result


def test_http_auth_origin_and_packaged_assets(http_server):
    server = http_server
    status, headers, body = request(server, "/")
    assert status == 200 and server.token.encode() in body
    assert "frame-ancestors 'none'" in headers["Content-Security-Policy"]
    assert request(server, "/app.js")[0] == 200
    assert request(server, "/style.css")[0] == 200
    assert request(server, "/api/settings")[0] == 403
    auth = {"X-EACBP-Token": server.token}
    assert request(server, "/api/settings", headers=auth)[0] == 200
    assert request(server, "/api/settings", headers={**auth, "Origin": "https://untrusted.invalid"})[0] == 403
    assert request(server, "/", headers={"Host": "untrusted.invalid"})[0] == 403
    assert request(server, "/api/settings", headers={**auth, "Sec-Fetch-Site": "cross-site"})[0] == 403
    assert request(server, "/api/files?path=..", headers=auth)[0] == 400
    assert request(server, "/api/download?id=x&name=../../pyproject.toml", headers=auth)[0] == 400
    assert request(server, "/api/preview", "POST", "[]", {**auth, "Content-Type": "application/json"})[0] == 400


def test_worker_real_run_restart_resume_and_report_gate(service, tmp_path):
    """Real h5ad + real subprocess, retained across a fresh manager instance."""
    source = input_data(tmp_path)
    prepared = form(source)
    result = service.start(prepared)
    first = wait_job(service.manager, result["job"]["id"])
    assert first["status"] == "success", first
    run_dir = Path(first["run_dir"])
    report = (run_dir / "report.md").read_bytes()
    restarted = JobManager(tmp_path, service.manager.runs_dir)
    assert restarted.runs()[0]["status"] == "success"
    source.unlink()  # Existing resume must not re-import the external source.
    resumed = restarted.submit("resume", run_dir=run_dir)
    assert wait_job(restarted, resumed["id"])["status"] == "success"
    detail = restarted.detail(result["run_id"])
    assert any(event["kind"] == "task_resumed" for event in detail["events"])
    assert not any(event["kind"] == "attempt_started" for event in detail["events"])
    rebuilt = restarted.submit("report", run_dir=run_dir)
    assert wait_job(restarted, rebuilt["id"])["status"] == "success"
    assert (run_dir / "report.md").read_bytes() == report
    from eacbp.artifact.registry import ArtifactRegistry
    registry = ArtifactRegistry(str(run_dir / "artifacts"))
    raw = registry.get_metadata("adata://ui_test/raw/v1")
    Path(raw.storage_path).write_bytes(b"corrupted payload")
    rejected = restarted.submit("report", run_dir=run_dir)
    failure = wait_job(restarted, rejected["id"])
    assert failure["status"] == "failed"
    assert (run_dir / "report.md").read_bytes() == report


def test_active_queued_job_prevents_double_submission(service):
    write_json(service.manager.jobs_dir / "queued" / "job.json", {
        "id": "queued", "status": "queued", "created_at": "", "created_epoch": time.time(), "run_dir": ""})
    with pytest.raises(BlockingIOError):
        service.manager.submit("run", run_dir=service.manager.runs_dir / "new")


def test_cli_webui_is_lazy_and_forwards_options(monkeypatch):
    from eacbp import cli
    captured = {}
    monkeypatch.setattr("eacbp.webui.server.serve", lambda **kwargs: captured.update(kwargs) or 0)
    assert cli.main(["webui", "--port", "9876", "--no-browser"]) == 0
    assert captured["port"] == 9876 and captured["open_browser"] is False
