"""Prepared-experiment recovery boundaries, durable jobs and bounded display.

These tests exercise dispatch and persistence without running scientific models.
"""
from copy import deepcopy
import http.client
from html.parser import HTMLParser
import io
import json
from pathlib import Path
import threading
import time

import pytest

from eacbp.application import recovery
from eacbp.application.execution_journal import ExecutionJournal
from eacbp.webui import worker
from eacbp.webui.jobs import (JobManager, RECOVERY_OPERATION, lease, read_json,
                              recovery_directory, recovery_progress, write_json)
from eacbp.webui.server import WebServer
from eacbp.webui.service import WebService


DIGEST = "a" * 64


@pytest.fixture
def service(tmp_path):
    return WebService(tmp_path, tmp_path / "runs")


@pytest.fixture
def prepared(service):
    path = service.manager.runs_dir / "prepared-experiment"
    path.mkdir()
    (path / "bindings.json").write_text("{}", encoding="utf-8")
    (path / "report.md").write_bytes(b"Original historical report\n")
    return path


def payload(prepared):
    return {"prepared": str(prepared), "expected_bindings_sha256": DIGEST}


def ready(prepared):
    return {"schema_version": 1, "status": "ready", "can_resume": True,
            "output_dir": str(prepared), "expected_bindings_sha256": DIGEST,
            "identity": {"verified": True}, "blockers": [], "journal": None,
            "branches": [{"branch_id": "primary", "native_exists": True,
                          "state": "pending_re_admission", "reuse_authorized": False,
                          "requires_native_readmission": True}],
            "biological_stability_claim_supported": False}


def tree_bytes(path):
    return {str(p.relative_to(path)): p.read_bytes() for p in path.rglob("*") if p.is_file()}


@pytest.fixture
def no_spawn(monkeypatch):
    class Process:
        def wait(self):
            return 0
    calls = []
    def popen(*args, **kwargs):
        calls.append((args, kwargs))
        return Process()
    monkeypatch.setattr("eacbp.webui.jobs.subprocess.Popen", popen)
    return calls


def queued(service, prepared, no_spawn):
    job = service.manager.submit(RECOVERY_OPERATION, run_dir=prepared,
                                 expected_bindings_sha256=DIGEST, recovery_inspection=ready(prepared))
    return job, service.manager.jobs_dir / job["id"]


@pytest.mark.parametrize("input_file", [False, True])
def test_inspection_is_read_only_and_does_not_submit(service, prepared, monkeypatch, input_file):
    seen = []
    def inspect(path, *, expected_bindings_sha256):
        seen.append((path, expected_bindings_sha256))
        return ready(path)
    monkeypatch.setattr(recovery, "inspect_recovery", inspect)
    monkeypatch.setattr(service.manager, "submit", lambda *a, **k: pytest.fail("Inspection submitted work"))
    before = tree_bytes(service.manager.runs_dir)
    target = prepared / "bindings.json" if input_file else prepared
    result = service.inspect_recovery(payload(target))
    assert seen == [(prepared, DIGEST)]
    assert result["branches"][0]["reuse_authorized"] is False
    assert result["biological_stability_claim_supported"] is False
    assert tree_bytes(service.manager.runs_dir) == before


def test_actual_shared_inspection_missing_lock_is_blocked_without_creating_it(service, prepared):
    before = tree_bytes(service.manager.runs_dir)
    report = service.inspect_recovery(payload(prepared))
    assert report["status"] == "blocked" and report["can_resume"] is False
    assert report["blockers"][0]["code"] == "executor_lock"
    assert tree_bytes(service.manager.runs_dir) == before


@pytest.mark.parametrize("sha", [None, "", "a" * 63, "A" * 64, "g" * 64, "a" * 65,
                                 "a" * 64 + "\n", " " + "a" * 64, 123])
def test_expected_sha_must_be_explicit_lowercase_hex(service, prepared, sha, monkeypatch):
    monkeypatch.setattr(recovery, "inspect_recovery", lambda *a, **k: pytest.fail("Invalid SHA reached inspection"))
    value = payload(prepared)
    if sha is None:
        value.pop("expected_bindings_sha256")
    else:
        value["expected_bindings_sha256"] = sha
    with pytest.raises(ValueError):
        service.inspect_recovery(value)
    assert not list(service.manager.jobs_dir.glob("*/job.json"))


@pytest.mark.parametrize("target", ["root", "webui", "outside", "missing", "otherfile"])
def test_only_existing_prepared_directory_inside_current_root(service, prepared, tmp_path, target):
    paths = {"root": service.manager.runs_dir, "webui": service.manager.jobs_dir,
             "outside": tmp_path, "missing": service.manager.runs_dir / "absent",
             "otherfile": prepared / "report.md"}
    with pytest.raises((ValueError, FileNotFoundError)):
        service.inspect_recovery(payload(paths[target]))


def test_blocked_inspection_never_enqueues(service, prepared, monkeypatch):
    report = {**ready(prepared), "status": "blocked", "can_resume": False,
              "blockers": [{"code": "source_identity", "reason": "source changed"}]}
    monkeypatch.setattr(recovery, "inspect_recovery", lambda *a, **k: report)
    before = tree_bytes(service.manager.runs_dir)
    assert service.resume_recovery(payload(prepared)) == report
    assert tree_bytes(service.manager.runs_dir) == before


def test_ready_enqueue_has_explicit_pins_and_is_not_a_fake_study(service, prepared, monkeypatch, no_spawn):
    inspection = {**ready(prepared), "journal": {"attempts": [{"provenance": "do not persist in poll snapshot"}]}}
    monkeypatch.setattr(recovery, "inspect_recovery", lambda *a, **k: inspection)
    response = service.resume_recovery(payload(prepared / "bindings.json"))
    job = response["job"]
    request = read_json(service.manager.jobs_dir / job["id"] / "request.json")
    assert job["operation"] == RECOVERY_OPERATION and job["status"] == "queued"
    assert request["expected_bindings_sha256"] == DIGEST
    assert request["approved_runs_root"] == str(service.manager.runs_dir)
    assert request["run_dir"] == str(prepared)
    assert "journal" not in request["recovery_inspection"]
    assert len(no_spawn) == 1
    assert service.manager.runs() == []
    with pytest.raises(BlockingIOError):
        service.manager.submit("run", run_dir=service.manager.runs_dir / "ordinary")
    assert not (prepared / "run_config.json").exists()


@pytest.mark.parametrize("status,expected_exit,expected_job", [
    ("completed", 0, "success"), ("partial", 1, "failed"),
    ("failed", 1, "failed"), ("blocked", 1, "failed"),
])
def test_worker_dispatch_preserves_actual_result_and_failure_reason(service, prepared, no_spawn, monkeypatch,
                                                                   status, expected_exit, expected_job):
    job, directory = queued(service, prepared, no_spawn)
    result = {"status": status, "biological_stability_claim_supported": False,
              "failures": [{"branch_id": "resolution_0_5", "reason": "independent audit refused"}],
              "run_provenance": {"large": {"artifacts": ["full evidence stays on disk"]}}}
    seen = []
    def resume(path, *, expected_bindings_sha256):
        seen.append((path, expected_bindings_sha256))
        return deepcopy(result)
    monkeypatch.setattr(recovery, "resume_recovery", resume)
    monkeypatch.setattr("eacbp.application.study_service.resume_study",
                        lambda **kwargs: pytest.fail("Recovery used single-study resume"))
    report_before = (prepared / "report.md").read_bytes()
    assert worker.run_job(directory) == expected_exit
    persisted = read_json(directory / "job.json")
    assert seen == [(prepared, DIGEST)]
    assert persisted["status"] == expected_job and persisted["recovery_status"] == status
    assert persisted["result"] == result
    assert persisted["biological_stability_claim_supported"] is False
    if status != "completed":
        assert "independent audit refused" in persisted["error"]
    displayed = service.manager.recovery_detail(job["id"])
    assert "run_provenance" not in displayed["job"]["result"]
    assert displayed["job"]["result"]["failures"] == result["failures"]
    assert (prepared / "report.md").read_bytes() == report_before


def test_worker_partial_reads_saved_failure_details_without_changing_comparison(service, prepared, no_spawn, monkeypatch):
    _, directory = queued(service, prepared, no_spawn)
    comparison = prepared / "comparisons" / "one" / "comparison.json"
    write_json(comparison, {"status": "partial", "failures": [{"reason": "missing admitted output"}],
                            "run_provenance": {"not_for_poll": "inventory"}})
    before = comparison.read_bytes()
    monkeypatch.setattr(recovery, "resume_recovery", lambda *a, **k: {
        "status": "partial", "comparison": str(comparison), "biological_stability_claim_supported": False})
    assert worker.run_job(directory) == 1
    persisted = read_json(directory / "job.json")
    assert "missing admitted output" in persisted["error"]
    assert persisted["recovery_failure_details"] == {"failures": [{"reason": "missing admitted output"}]}
    assert comparison.read_bytes() == before


@pytest.mark.parametrize("change", ["root", "path", "sha", "missing_path", "run", "resume", "report"])
def test_worker_revalidates_queue_root_path_and_sha_before_dispatch(service, prepared, no_spawn, monkeypatch, change):
    _, directory = queued(service, prepared, no_spawn)
    request = read_json(directory / "request.json")
    if change == "root":
        request["approved_runs_root"] = str(service.manager.workspace)
    elif change == "path":
        request["run_dir"] = str(service.manager.workspace)
    elif change == "sha":
        request["expected_bindings_sha256"] = "b" * 64
    elif change == "missing_path":
        prepared.rename(prepared.with_name("moved"))
    else:
        request["operation"] = change
        request["webui_context"] = {"must_not_write_after_rejection": True}
    write_json(directory / "request.json", request)
    monkeypatch.setattr(recovery, "resume_recovery", lambda *a, **k: pytest.fail("Changed request executed"))
    for api in ("run_study", "resume_study", "report_study"):
        monkeypatch.setattr("eacbp.application.study_service." + api,
                            lambda **k: pytest.fail("Changed request executed ordinary API"))
    assert worker.run_job(directory) == 1
    persisted = read_json(directory / "job.json")
    assert persisted["status"] == "failed" and persisted["error"]
    assert "result" not in persisted
    if prepared.exists():
        assert not (prepared / "webui_context.json").exists()


def test_actual_shared_worker_rejects_changed_bindings_before_executor(service, prepared, no_spawn, monkeypatch):
    _, directory = queued(service, prepared, no_spawn)
    (prepared / ".execution-journal.lock").write_bytes(b"0")
    def changed(*a, **k):
        raise ValueError("Pinned bindings changed after enqueue")
    monkeypatch.setattr(recovery.execution, "_load_preparation", changed)
    monkeypatch.setattr(recovery.execution, "execute_clustering_sensitivity",
                        lambda *a, **k: pytest.fail("Changed bindings executed"))
    assert worker.run_job(directory) == 1
    persisted = read_json(directory / "job.json")
    assert persisted["result"]["status"] == "blocked"
    assert "Pinned bindings changed" in persisted["error"]


def test_real_detached_worker_persists_blocked_dispatch_without_scientific_execution(service, prepared):
    before = tree_bytes(prepared)
    job = service.manager.submit(RECOVERY_OPERATION, run_dir=prepared,
                                 expected_bindings_sha256=DIGEST)
    deadline = time.monotonic() + 30
    directory = service.manager.jobs_dir / job["id"]
    while time.monotonic() < deadline:
        persisted = read_json(directory / "job.json")
        if persisted["status"] not in {"queued", "running"}:
            break
        time.sleep(.1)
    assert persisted["status"] == "failed"
    assert persisted["pid"] > 0 and persisted["started_at"] and persisted["finished_at"]
    assert persisted["result"]["status"] == "blocked"
    assert persisted["result"]["blockers"][0]["code"] == "executor_lock"
    assert "Finished: blocked" in (directory / "worker.log").read_text(encoding="utf-8")
    assert tree_bytes(prepared) == before
    restarted = JobManager(service.manager.workspace, service.manager.runs_dir)
    assert restarted.recovery_detail(job["id"])["job"]["error"] == persisted["error"]


def journal_fixture(prepared):
    identity = {key: DIGEST for key in ("bindings_sha256", "plan_file_sha256", "plan_sha256", "source_sha256")}
    identity["execution_environment"] = {"source_sha256": DIGEST}
    with ExecutionJournal(prepared, identity=identity) as journal:
        failed = journal.start_attempt("postprocessing")
        journal.finish_attempt(failed, status="postprocessing_failed", details={
            "reason": "report write fault", "report": "comparison/report.md",
            "provenance": {"artifacts": ["do not display"]}})
        journal.start_attempt("native_admission", "primary")


def test_progress_keeps_unfinished_and_failure_reasons_without_repair_or_fullinspect(prepared, monkeypatch):
    journal_fixture(prepared)
    head = prepared / "execution_journal" / "head.json"
    saved = read_json(head)
    saved["sequence"] = 0
    write_json(head, saved)
    before = tree_bytes(prepared)
    monkeypatch.setattr(recovery, "inspect_recovery", lambda *a, **k: pytest.fail("Polling fully inspected"))
    progress = recovery_progress(prepared)
    assert progress["integrity_scope"] == "diagnostic_only_not_native_admission"
    assert [item["status"] for item in progress["attempts"]] == ["postprocessing_failed", "started"]
    assert progress["attempts"][0]["finish"]["details"]["reason"] == "report write fault"
    assert "provenance" not in progress["attempts"][0]["finish"]["details"]
    assert progress["attempts"][1]["finish"] is None
    assert "elapsed_seconds" not in progress["attempts"][1]["start"]
    assert tree_bytes(prepared) == before


def test_progress_head_read_failure_still_displays_events(prepared):
    journal_fixture(prepared)
    (prepared / "execution_journal" / "head.json").write_text("not JSON", encoding="utf-8")
    before = tree_bytes(prepared)
    progress = recovery_progress(prepared)
    assert progress["warnings"] and len(progress["attempts"]) == 2
    assert tree_bytes(prepared) == before


@pytest.mark.parametrize("mode", ["single_large", "total_budget"])
def test_progress_limits_disk_budget_and_excludes_inventory(prepared, mode):
    events = prepared / "execution_journal" / "events"
    events.mkdir(parents=True)
    write_json(events.parent / "head.json", {"sequence": 20})
    count, inventory_size = (1, 1_100_000) if mode == "single_large" else (20, 600_000)
    for sequence in range(1, count + 1):
        write_json(events / f"{sequence:08d}.json", {"kind": "attempt_started", "stage": "native_execution",
                   "attempt_id": str(sequence), "status": "started", "details": {
                       "provenance": {"files": {"large": "x" * inventory_size}}}})
    progress = recovery_progress(prepared)
    assert progress["truncated"] is True and progress["warnings"]
    assert len(json.dumps(progress)) < 20_000
    if mode == "total_budget":
        assert 0 < len(progress["attempts"]) < count
        assert any("8 MB" in warning for warning in progress["warnings"])


def test_root_switch_does_not_expose_old_recovery_job(service, prepared, no_spawn, tmp_path):
    job, directory = queued(service, prepared, no_spawn)
    saved = read_json(directory / "job.json")
    saved.update(status="failed", error="durable original failure")
    write_json(directory / "job.json", saved)
    before = (prepared / "report.md").read_bytes()
    new_workspace = tmp_path / "new-workspace"
    new_workspace.mkdir()
    service.update_directories({"workspace": str(new_workspace), "runs_dir": str(tmp_path / "new-runs")})
    with pytest.raises(ValueError):
        service.inspect_recovery(payload(prepared))
    with pytest.raises(FileNotFoundError):
        service.manager.recovery_detail(job["id"])
    assert service.manager.jobs() == []
    assert (prepared / "report.md").read_bytes() == before
    assert read_json(directory / "job.json")["error"] == "durable original failure"


def test_progress_read_is_bounded_if_event_grows_after_stat(prepared, monkeypatch):
    events = prepared / "execution_journal" / "events"
    events.mkdir(parents=True)
    write_json(events.parent / "head.json", {"sequence": 1})
    path = events / "00000001.json"
    write_json(path, {"kind": "attempt_started", "attempt_id": "one", "stage": "preflight"})
    old_size = path.stat().st_size
    before = tree_bytes(prepared)
    requested = []
    class Growing(io.BytesIO):
        def read(self, size=-1):
            requested.append(size)
            return super().read(size)
    original_open = Path.open
    def open_path(self, *args, **kwargs):
        if self == path and args == ("rb",):
            return Growing(b"x" * 2_000_000)
        return original_open(self, *args, **kwargs)
    monkeypatch.setattr(Path, "open", open_path)
    progress = recovery_progress(prepared)
    assert requested == [old_size + 1]
    assert progress["truncated"] and any("读取时改变" in item for item in progress["warnings"])
    monkeypatch.setattr(Path, "open", original_open)
    assert tree_bytes(prepared) == before


@pytest.fixture
def http_server(service):
    server = WebServer(service, port=0)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield server
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def request(server, path, *, token=None, payload=None, origin=None):
    connection = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=15)
    headers = {}
    if token:
        headers["X-EACBP-Token"] = token
    if origin:
        headers["Origin"] = origin
    body = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        body = json.dumps(payload)
    connection.request("GET" if payload is None else "POST", path, body=body, headers=headers)
    response = connection.getresponse()
    status, body = response.status, response.read()
    connection.close()
    return status, body


@pytest.mark.parametrize("route", ["/api/recovery/inspect", "/api/recovery/resume", "/api/recovery/job?id=" + "a" * 32])
def test_recovery_routes_keep_token_and_origin_boundary(http_server, prepared, route):
    data = None if "/job?" in route else payload(prepared)
    assert request(http_server, route, payload=data)[0] == 403
    assert request(http_server, route, token=http_server.token, payload=data, origin="https://outside.example")[0] == 403


def test_http_inspect_enqueue_job_and_no_duplicate_job_id(http_server, prepared, monkeypatch, no_spawn):
    monkeypatch.setattr(recovery, "inspect_recovery", lambda path, **k: ready(path))
    server = http_server
    before = tree_bytes(prepared)
    status, body = request(server, "/api/recovery/inspect", token=server.token, payload=payload(prepared))
    assert status == 200 and json.loads(body)["can_resume"] is True
    assert tree_bytes(prepared) == before and server.service.manager.jobs() == []
    status, body = request(server, "/api/recovery/resume", token=server.token, payload=payload(prepared))
    assert status == 202
    job_id = json.loads(body)["job"]["id"]
    status, body = request(server, "/api/recovery/job?id=" + job_id, token=server.token)
    assert status == 200 and json.loads(body)["job"]["status"] == "queued"
    assert request(server, "/api/recovery/job?id=" + job_id + "&id=" + job_id, token=server.token)[0] == 400
    assert request(server, "/api/recovery/job?id=../escape", token=server.token)[0] == 400
    assert request(server, "/api/recovery/resume", token=server.token, payload=payload(prepared))[0] == 409


def test_http_blocked_result_does_not_look_enqueued(http_server, prepared):
    status, body = request(http_server, "/api/recovery/resume", token=http_server.token, payload=payload(prepared))
    report = json.loads(body)
    assert status == 200 and report["status"] == "blocked" and "job" not in report
    assert http_server.service.manager.jobs() == []


def test_packaged_recovery_form_requires_sha_and_explains_diagnostic_scope():
    static = Path(__file__).parents[1] / "eacbp" / "webui" / "static"
    class Elements(HTMLParser):
        def __init__(self):
            super().__init__()
            self.ids = {}
        def handle_starttag(self, tag, attrs):
            fields = dict(attrs)
            if "id" in fields:
                assert fields["id"] not in self.ids
                self.ids[fields["id"]] = fields
    parser = Elements()
    html = (static / "index.html").read_text(encoding="utf-8")
    parser.feed(html)
    sha = parser.ids["recovery-bindings-sha"]
    assert "required" in sha and sha["pattern"] == "[0-9a-f]{64}" and "value" not in sha
    assert "disabled" in parser.ids["resume-recovery"]
    assert "原生分支" in html and "外部科学结论需要独立验证" in html
    script = (static / "app.js").read_text(encoding="utf-8")
    assert "未记录结束" in script and "/api/recovery/job?id=" in script
    assert "state.recoveryKey!==JSON.stringify(recoveryPayload())" in script
