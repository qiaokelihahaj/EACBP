"""Directory changes preserve active work and invalidate stale browser contexts."""
from concurrent.futures import ThreadPoolExecutor
import http.client
import json
from pathlib import Path
import threading
import time

import pytest

from eacbp.webui.jobs import JobManager, write_json
from eacbp.webui.preferences import startup_directories
from eacbp.webui.server import WebServer
from eacbp.webui.service import WebService


@pytest.fixture
def service(tmp_path):
    original = tmp_path / "original_data"
    original.mkdir()
    return WebService(original, tmp_path / "original_results", settings_file=tmp_path / "config" / "webui.json")


def directories(tmp_path):
    workspace = tmp_path / "separate data"
    workspace.mkdir(exist_ok=True)
    return {"workspace": str(workspace), "runs_dir": str(tmp_path / "separate results")}


def queued_job(manager):
    write_json(manager.jobs_dir / "pending" / "job.json", {
        "id": "pending", "status": "queued", "created_at": "2026-09-22T00:00:00Z",
        "created_epoch": time.time(), "run_dir": str(manager.runs_dir / "pending"),
    })


def test_save_both_roots_and_reload_without_moving_history(service, tmp_path):
    previous = service.manager
    old_data = previous.workspace / "original.h5ad"
    old_data.write_bytes(b"untouched input")
    history = previous.runs_dir / "old" / "run_config.json"
    write_json(history, {"study_id": "old", "status": "success"})
    before = history.read_bytes()
    paths = directories(tmp_path)
    new_data = Path(paths["workspace"]) / "selected.h5ad"
    new_data.write_bytes(b"new input")
    result = service.update_directories(paths)
    assert result == {**paths, "saved": True}
    assert service.manager.workspace == Path(paths["workspace"])
    assert service.manager.runs_dir == Path(paths["runs_dir"])
    assert service.data_path("selected.h5ad") == new_data
    with pytest.raises(ValueError, match="outside"):
        service.data_path(old_data)
    assert history.read_bytes() == before and old_data.read_bytes() == b"untouched input"
    assert service.manager.runs() == []
    workspace, runs, settings = startup_directories(settings_file=service.settings_file)
    assert workspace == service.manager.workspace and runs == service.manager.runs_dir
    assert settings == service.settings_file
    assert service.settings()["package_dir"] != str(workspace)
    service.update_directories({"workspace": str(previous.workspace), "runs_dir": str(previous.runs_dir)})
    assert service.manager.runs()[0]["study_id"] == "old"


def test_startup_defaults_and_cli_precedence(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    settings = tmp_path / "preferences.json"
    workspace, runs, _ = startup_directories(settings_file=settings)
    assert workspace == tmp_path and runs == tmp_path / "outputs" / "runs"
    assert not settings.exists()
    saved = directories(tmp_path)
    write_json(settings, {"schema_version": 1, **saved})
    workspace, runs, _ = startup_directories(workspace=tmp_path, settings_file=settings)
    assert workspace == tmp_path and runs == Path(saved["runs_dir"])
    workspace, runs, _ = startup_directories(runs_dir=tmp_path / "explicit", settings_file=settings)
    assert workspace == Path(saved["workspace"]) and runs == tmp_path / "explicit"


@pytest.mark.parametrize("payload", ["not json", "[]", '{"schema_version":2}', '{"schema_version":true}'])
def test_invalid_preferences_are_not_silently_ignored(tmp_path, payload):
    path = tmp_path / "preferences.json"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(ValueError, match="目录设置"):
        startup_directories(settings_file=path)


@pytest.mark.parametrize("change", [
    {"workspace": "relative"}, {"runs_dir": "relative"}, {"workspace": ""},
    {"runs_dir": None}, {"unknown": "value"},
])
def test_invalid_fields_leave_roots_unchanged(service, tmp_path, change):
    manager = service.manager
    with pytest.raises(ValueError):
        service.update_directories({**directories(tmp_path), **change})
    assert service.manager is manager
    assert not service.settings_file.exists()


def test_missing_input_or_file_output_leaves_roots_unchanged(service, tmp_path):
    manager = service.manager
    paths = directories(tmp_path)
    with pytest.raises(ValueError, match="输入数据目录不存在"):
        service.update_directories({**paths, "workspace": str(tmp_path / "missing")})
    output = Path(paths["runs_dir"])
    output.write_text("do not replace", encoding="utf-8")
    with pytest.raises(OSError):
        service.update_directories(paths)
    assert service.manager is manager and not service.settings_file.exists()
    assert output.read_text(encoding="utf-8") == "do not replace"


@pytest.mark.parametrize("target", [False, True])
@pytest.mark.parametrize("corrupt", [False, True])
def test_active_or_unknown_jobs_block_switch_at_either_root(service, tmp_path, target, corrupt):
    paths = directories(tmp_path)
    manager = service.manager
    occupied = JobManager(Path(paths["workspace"]), Path(paths["runs_dir"])) if target else manager
    queued_job(occupied)
    if corrupt:
        (occupied.jobs_dir / "pending" / "job.json").write_text("broken", encoding="utf-8")
    with pytest.raises(BlockingIOError, match="状态待检查"):
        service.update_directories(paths)
    assert service.manager is manager and not service.settings_file.exists()


def test_failed_persistence_does_not_apply_new_roots(service, tmp_path, monkeypatch):
    manager = service.manager
    write_json(service.settings_file, {"schema_version": 1, "workspace": str(manager.workspace), "runs_dir": str(manager.runs_dir)})
    before = service.settings_file.read_bytes()

    def fail(*args, **kwargs):
        raise PermissionError("preference destination is read-only")

    monkeypatch.setattr("eacbp.webui.service.write_json", fail)
    with pytest.raises(PermissionError):
        service.update_directories(directories(tmp_path))
    assert service.manager is manager and service.settings_file.read_bytes() == before


@pytest.fixture
def http_server(service):
    server = WebServer(service, port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield server
    server.shutdown()
    server.server_close()
    thread.join(timeout=5)


def request(server, path, *, token=None, payload=None, origin=None):
    headers = {"X-EACBP-Token": token} if token else {}
    if origin:
        headers["Origin"] = origin
    if payload is not None:
        headers["Content-Type"] = "application/json"
    connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=15)
    connection.request("POST" if payload is not None else "GET", path,
                       body=json.dumps(payload) if payload is not None else None, headers=headers)
    response = connection.getresponse()
    result = response.status, response.read()
    connection.close()
    return result


def test_settings_save_invalidates_old_tabs_and_preserves_auth_boundary(http_server, tmp_path):
    server = http_server
    paths = directories(tmp_path)
    old_token = server.token
    assert request(server, "/api/settings", payload=paths)[0] == 403
    assert request(server, "/api/settings", token=old_token, payload=paths, origin="https://example.invalid")[0] == 403
    assert request(server, "/api/settings", token=old_token, payload={**paths, "workspace": "relative"})[0] == 400
    assert server.token == old_token
    assert request(server, "/api/settings", token=old_token, payload=paths)[0] == 200
    assert server.token != old_token
    assert request(server, "/api/files", token=old_token)[0] == 403
    assert request(server, "/api/run", token=old_token, payload={})[0] == 403
    assert server.token.encode() in request(server, "/")[1]
    status, data = request(server, "/api/settings", token=server.token)
    assert status == 200 and json.loads(data)["workspace"] == paths["workspace"]


def test_two_old_tabs_cannot_apply_conflicting_roots(http_server, tmp_path):
    server = http_server
    token = server.token
    a = directories(tmp_path)
    b = {**a, "runs_dir": str(tmp_path / "alternative_results")}
    with ThreadPoolExecutor(max_workers=2) as pool:
        calls = [pool.submit(request, server, "/api/settings", token=token, payload=payload) for payload in [a, b]]
        assert sorted(call.result()[0] for call in calls) == [200, 403]
    saved = json.loads(server.service.settings_file.read_text(encoding="utf-8"))
    assert saved["runs_dir"] == str(server.service.manager.runs_dir)


def test_job_submission_and_directory_save_are_serialized(http_server, tmp_path, monkeypatch):
    server = http_server
    started, release = threading.Event(), threading.Event()
    paths = directories(tmp_path)
    previous = server.service.manager

    def start(payload):
        started.set()
        assert release.wait(timeout=5)
        queued_job(previous)
        return {"run_id": "pending"}

    monkeypatch.setattr(server.service, "_start", start)
    token = server.token
    with ThreadPoolExecutor(max_workers=2) as pool:
        job = pool.submit(request, server, "/api/run", token=token, payload={})
        assert started.wait(timeout=5)
        settings = pool.submit(request, server, "/api/settings", token=token, payload=paths)
        release.set()
        assert job.result()[0] == 202
        assert settings.result()[0] == 409
    assert server.service.manager is previous and server.token == token


def test_new_worker_uses_selected_data_and_output_roots(service, tmp_path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    paths = directories(tmp_path)
    source = Path(paths["workspace"]) / "SYNTHETIC_ONLY.h5ad"
    data = ad.AnnData(np.random.default_rng(24).poisson(4, (80, 40)).astype(np.float32),
                      obs=pd.DataFrame(index=[f"c{i}" for i in range(80)]),
                      var=pd.DataFrame(index=[f"g{i}" for i in range(40)]))
    data.uns["is_simulated"] = True
    data.write_h5ad(source)
    before = source.read_bytes()
    service.update_directories(paths)
    response = service.start({"data": "SYNTHETIC_ONLY.h5ad", "study_id": "directory_test",
                              "species": "human", "tissue": "synthetic_test",
                              "method_profile": "baseline", "min_genes": 1})
    deadline = time.monotonic() + 60
    while time.monotonic() < deadline:
        job = next(j for j in service.manager.jobs() if j["id"] == response["job"]["id"])
        if job["status"] not in {"running", "queued"}:
            break
        time.sleep(.1)
    else:
        pytest.fail("directory validation worker timed out")
    assert job["status"] == "success", job
    assert Path(job["run_dir"]).is_relative_to(Path(paths["runs_dir"]))
    saved_request = json.loads((service.manager.jobs_dir / job["id"] / "request.json").read_text(encoding="utf-8"))
    assert saved_request["workspace"] == paths["workspace"]
    assert source.read_bytes() == before
