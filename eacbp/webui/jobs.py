"""Persistent, process-isolated jobs for the local WebUI.

The browser server never executes scientific work in an HTTP request. Workers
use the package CLI API, including its original run locks and audit gates.
OS file leases distinguish a surviving worker from a stale running record.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import threading
import time
from uuid import uuid4

from eacbp._atomic_json import atomic_write_json


RECOVERY_OPERATION = "recover_clustering"


def study_preview_binding(*, data, source_sha256, manifest, config, plan):
    """Bind the confirmed source and effective plan without trusting browser hashes."""
    from eacbp.schemas.audit import stable_fingerprint

    if not isinstance(source_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", source_sha256) is None:
        raise ValueError("分析输入的 SHA-256 绑定缺失或无效，请重新检查方案。")
    return {"schema_version": 1, "source_sha256": source_sha256,
            "fingerprint": stable_fingerprint({"data": str(Path(data).resolve()),
                "source_sha256": source_sha256, "manifest": manifest, "config": config, "plan": plan})}


def recovery_summary(value):
    """Small display projection; full executor evidence remains in job.json."""
    keys = ("schema_version", "status", "can_resume", "output_dir", "expected_bindings_sha256",
            "blockers", "error", "reason", "failures", "native_failures", "postprocessing_failures",
            "native_reused_branch_ids", "native_started_branch_ids", "comparison_path", "report_path",
            "comparison", "report", "primary_native_summary", "completion_scope",
            "biological_stability_claim_supported")
    return {key: value[key] for key in keys if key in value} if isinstance(value, dict) else None


def recovery_display_event(event):
    keys = ("attempt_id", "stage", "branch_id", "status", "kind", "sequence", "session_id",
            "timestamp_utc", "identity_sha256", "event_sha256", "elapsed_seconds")
    details_keys = ("reason", "error", "error_type", "admitted", "reused", "path", "sha256",
                    "report", "report_sha256", "status", "json_path", "report_path", "comparison", "run_dir")
    result = {key: event[key] for key in keys if key in event}
    details = event.get("details", {})
    if not isinstance(details, dict):
        raise ValueError("阶段事件 details 必须是 JSON 对象")
    result["details"] = {key: details[key] for key in details_keys if key in details}
    return result


def recovery_directory(runs_root: Path, prepared: str | Path) -> Path:
    """Resolve an existing prepared target within the currently approved root."""
    if not isinstance(prepared, (str, Path)) or not str(prepared).strip():
        raise ValueError("请选择已准备的实验目录或 bindings.json")
    root = runs_root.resolve()
    path = contained(root, prepared)
    if path.is_file():
        if path.name != "bindings.json":
            raise ValueError("仅支持已准备的目录或 bindings.json")
        path = path.parent
    if path == root or path.is_relative_to(root / ".webui"):
        raise ValueError("请选择结果目录内的独立已准备实验")
    if not path.is_dir():
        raise FileNotFoundError("已准备的实验目录不存在")
    return path


def recovery_sha(value: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError("expected_bindings_sha256 必须显式提供 64 位小写十六进制 SHA-256")
    return value


def recovery_progress(prepared: Path):
    """Bounded, read-only journal display; never an execution/admission gate.

    The full inspection is performed explicitly and by the worker. Polling only
    reads published diagnostic metadata, without repairing journal heads or
    loading scientific payloads. A concurrently published event can appear on
    the next poll.
    """
    journal = contained(prepared, prepared / "execution_journal")
    result = {"path": str(journal), "integrity_scope": "diagnostic_only_not_native_admission",
              "attempts": [], "warnings": [], "truncated": False}
    if not journal.is_dir():
        result["warnings"].append("尚无持久阶段记录；检查或恢复时会核验实际状态。")
        return result
    budget = 8_000_000
    def bounded_json(path):
        nonlocal budget
        path = contained(prepared, path)
        size = path.stat().st_size
        if size > 1_000_000:
            raise ValueError("阶段记录超过显示大小限制")
        if size > budget:
            result["truncated"] = True
            raise ValueError("阶段轮询已达到 8 MB 总读盘预算")
        budget -= size
        with path.open("rb") as stream:
            data = stream.read(size + 1)
        if len(data) > size:
            raise ValueError("阶段记录在读取时改变")
        value = json.loads(data.decode("utf-8"))
        if not isinstance(value, dict):
            raise ValueError("阶段记录必须是 JSON 对象")
        return value
    try:
        result["head"] = bounded_json(journal / "head.json")
    except (OSError, ValueError, TypeError) as exc:
        result["warnings"].append(f"阶段 head 暂无法显示：{type(exc).__name__}: {exc}")
    try:
        # Sequence filenames are published once by the shared executor.
        paths = sorted(p for p in (journal / "events").glob("*.json")
                       if re.fullmatch(r"[0-9]{8}\.json", p.name))
        result["event_count"] = len(paths)
        result["truncated"] = len(paths) > 300
        attempts = {}
        # Read newest records first so a bounded poll retains current failures.
        events = []
        for path in reversed(paths[-300:]):
            try:
                events.append(recovery_display_event(bounded_json(path)))
            except (OSError, ValueError, TypeError) as exc:
                result["truncated"] = True
                result["warnings"].append(f"阶段事件显示已截断：{type(exc).__name__}: {exc}")
                break
        for event in reversed(events):
            attempt_id = event.get("attempt_id")
            if not isinstance(attempt_id, str) or event.get("kind") not in {"attempt_started", "attempt_finished"}:
                raise ValueError("阶段事件格式无法显示")
            item = attempts.setdefault(attempt_id, {"attempt_id": attempt_id,
                "stage": event.get("stage"), "branch_id": event.get("branch_id"),
                "status": "started", "start": None, "finish": None})
            if event["kind"] == "attempt_started":
                item["start"] = event
            else:
                item.update(status=event.get("status"), finish=event)
        result["attempts"] = list(attempts.values())
        if len(paths) > 300:
            result["warnings"].append("仅显示最近 300 条事件；完整记录保存在实验目录。")
    except (OSError, ValueError, TypeError) as exc:
        result["warnings"].append(f"持久阶段记录暂无法完整显示：{type(exc).__name__}: {exc}")
    return result


def read_json(path: Path, default=None):
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value):
    atomic_write_json(path, value, create_parent=True,
                      dump_kwargs={"ensure_ascii": False, "indent": 2, "allow_nan": False})


def now():
    return datetime.now(timezone.utc).isoformat()


def contained(root: Path, value: str | Path) -> Path:
    candidate = Path(value).expanduser()
    path = (candidate if candidate.is_absolute() else root / candidate).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("路径必须位于配置的目录内 / path is outside the configured root")
    return path


@contextmanager
def lease(path: Path):
    """Nonblocking cross-process lease, released by the OS on process death."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        if handle.tell() == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise BlockingIOError("已有操作正在运行，请等待完成") from exc
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def tail(path: Path, limit=96_000) -> str:
    if not path.is_file():
        return ""
    with path.open("rb") as stream:
        size = stream.seek(0, os.SEEK_END)
        stream.seek(max(0, size - limit))
        data = stream.read(limit)
    if size > limit:
        data = data.partition(b"\n")[2]
    return data.decode("utf-8", errors="replace")


class JobManager:
    def __init__(self, workspace: Path, runs_dir: Path):
        self.workspace = workspace.resolve()
        self.runs_dir = runs_dir.resolve()
        self.jobs_dir = self.runs_dir / ".webui" / "jobs"
        self.jobs_dir.mkdir(parents=True, exist_ok=True)
        self._guard = threading.Lock()

    def jobs(self):
        result = []
        for path in self.jobs_dir.glob("*/job.json"):
            if not path.resolve().is_relative_to(self.jobs_dir.resolve()):
                continue
            try:
                job = read_json(path)
                if job["status"] in {"queued", "running"}:
                    try:
                        with lease(path.parent / "worker.lock"):
                            # Re-read under the lease: a worker may just have finished.
                            job = read_json(path)
                            if job["status"] == "running" or (
                                job["status"] == "queued" and time.time() - job["created_epoch"] > 120
                            ):
                                job.update(status="interrupted", finished_at=now(),
                                           error="后台进程已退出；保留了结果，请检查日志后恢复。")
                                write_json(path, job)
                    except BlockingIOError:
                        pass
                if job.get("operation") == RECOVERY_OPERATION and "result" in job:
                    job = {**job, "result": recovery_summary(job["result"]),
                           "result_scope": "display_summary_full_result_in_job_json"}
                result.append(job)
            except (OSError, ValueError, KeyError, TypeError):
                # A corrupt job record must not permit another computation.
                result.append({"id": path.parent.name, "status": "unknown",
                               "error": "无法读取后台作业记录，请检查 " + str(path),
                               "created_at": "", "run_dir": ""})
        return sorted(result, key=lambda item: item["created_at"], reverse=True)

    def submit(self, operation: str, *, run_dir: Path, manifest=None, config=None, data=None,
               webui_context=None, expected_bindings_sha256=None, recovery_inspection=None,
               preview_binding=None):
        if operation not in {"run", "resume", "report", RECOVERY_OPERATION}:
            raise ValueError("Unknown operation")
        run_dir = contained(self.runs_dir, run_dir)
        if operation == RECOVERY_OPERATION:
            run_dir = recovery_directory(self.runs_dir, run_dir)
            recovery_sha(expected_bindings_sha256)
        if run_dir == self.runs_dir:
            raise ValueError("请选择独立运行目录")
        with self._guard, lease(self.jobs_dir.parent / "submit.lock"):
            if any(job["status"] in {"queued", "running", "unknown"} for job in self.jobs()):
                raise BlockingIOError("已有后台操作正在运行。首版每次只运行一个作业。")
            if operation == "run" and run_dir.exists():
                raise ValueError("运行目录已存在，请创建新运行")
            if operation not in {"run", RECOVERY_OPERATION}:
                saved = read_json(run_dir / "run_config.json", {})
                if saved.get("import_completed") is not True:
                    raise ValueError("该运行未完成数据导入，不能恢复或重建报告")
            job_id = uuid4().hex
            directory = self.jobs_dir / job_id
            directory.mkdir()
            job = {"id": job_id, "operation": operation, "status": "queued",
                   "run_dir": str(run_dir), "created_at": now(), "created_epoch": time.time()}
            if operation == "run":
                job["preview_binding"] = preview_binding
            if operation == RECOVERY_OPERATION:
                job["expected_bindings_sha256"] = expected_bindings_sha256
            write_json(directory / "request.json", {
                "operation": operation, "run_dir": str(run_dir), "workspace": str(self.workspace),
                "manifest": manifest, "config": config, "data": str(data) if data else None,
                "webui_context": webui_context,
                "preview_binding": preview_binding,
                "approved_runs_root": str(self.runs_dir),
                "expected_bindings_sha256": expected_bindings_sha256,
                "recovery_inspection": recovery_summary(recovery_inspection),
            })
            write_json(directory / "job.json", job)
            try:
                env = dict(os.environ, PYTHONUNBUFFERED="1", PYTHONIOENCODING="utf-8")
                # Supports source checkouts as well as installed wheels, even when
                # the selected data workspace is not the source directory.
                package_parent = str(Path(__file__).resolve().parents[2])
                env["PYTHONPATH"] = package_parent + os.pathsep + env.get("PYTHONPATH", "")
                kwargs = {"creationflags": subprocess.CREATE_NO_WINDOW} if os.name == "nt" else {"start_new_session": True}
                with (directory / "worker.log").open("ab") as log:
                    process = subprocess.Popen(
                        [sys.executable, "-m", "eacbp.webui.worker", str(directory)],
                        cwd=self.workspace, env=env, stdin=subprocess.DEVNULL,
                        stdout=log, stderr=subprocess.STDOUT, **kwargs,
                    )
                # Reap children without tying their lifetime to the web server.
                threading.Thread(target=process.wait, daemon=True).start()
            except Exception as exc:
                job.update(status="failed", error=str(exc), finished_at=now())
                write_json(directory / "job.json", job)
                raise
            return job

    def new_run_dir(self, study_id: str):
        name = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + study_id + "-" + uuid4().hex[:8]
        return self.runs_dir / name

    def run_path(self, run_id: str):
        path = contained(self.runs_dir, run_id)
        if path == self.runs_dir or path.is_relative_to(self.jobs_dir.parent):
            raise ValueError("Invalid run directory")
        return path

    def runs(self):
        jobs = self.jobs()
        paths = {Path(j["run_dir"]) for j in jobs if j.get("run_dir")
                 and j.get("operation") != RECOVERY_OPERATION}
        # Historical CLI and compatibility-script layouts; never recurse into
        # large artifact trees to discover runs.
        for pattern in ("*/run_config.json", "*/*/run_config.json", "*/*/*/run_config.json"):
            paths.update(p.parent for p in self.runs_dir.glob(pattern))
        rows = []
        for path in paths:
            try:
                path = self.run_path(str(path))
                saved = read_json(path / "run_config.json", {})
                latest = next((j for j in jobs if j.get("run_dir") == str(path)), None)
                status = saved.get("status", "queued")
                if latest and latest["status"] in {"queued", "running", "interrupted", "unknown"}:
                    status = latest["status"]
                elif latest and latest["operation"] != "report" and latest["status"] == "failed":
                    status = "failed"
                rows.append({"id": path.relative_to(self.runs_dir).as_posix(), "path": str(path),
                             "study_id": saved.get("study_id", path.name),
                             "title": saved.get("manifest", {}).get("title", ""), "status": status,
                             "updated_at": latest["created_at"] if latest else datetime.fromtimestamp(
                                 (path / "run_config.json").stat().st_mtime, timezone.utc).isoformat(),
                             "job": latest, "can_resume": saved.get("import_completed") is True,
                             "can_report": saved.get("import_completed") is True and (path / "snapshot.json").is_file(),
                             "has_report": (path / "report.md").is_file()})
            except (OSError, ValueError, TypeError):
                rows.append({"id": str(path), "study_id": path.name, "status": "unreadable",
                             "updated_at": "", "error": "无法读取运行配置"})
        return sorted(rows, key=lambda row: row["updated_at"], reverse=True)

    def recovery_detail(self, job_id: str):
        if not isinstance(job_id, str) or re.fullmatch(r"[0-9a-f]{32}", job_id) is None:
            raise ValueError("Invalid recovery job id")
        directory = contained(self.jobs_dir, self.jobs_dir / job_id)
        job = next((item for item in self.jobs() if item["id"] == job_id), None)
        if job is None:
            raise FileNotFoundError("恢复作业记录不存在")
        if job.get("operation") != RECOVERY_OPERATION:
            raise ValueError("该作业不是实验恢复作业")
        request = read_json(directory / "request.json", {})
        try:
            prepared = recovery_directory(self.runs_dir, job["run_dir"])
            progress = recovery_progress(prepared)
        except (OSError, ValueError, TypeError) as exc:
            progress = {"attempts": [], "warnings": [str(exc)],
                        "integrity_scope": "diagnostic_only_not_native_admission"}
        return {"job": job, "inspection": request.get("recovery_inspection"),
                "progress": progress, "log": tail(directory / "worker.log"),
                "biological_stability_claim_supported": False}

    def detail(self, run_id: str):
        path = self.run_path(run_id)
        row = next((row for row in self.runs() if row["id"] == run_id), None)
        if row is None:
            raise FileNotFoundError("运行记录不存在")
        events = []
        latest = row.get("job")
        event_paths = list((path / "artifacts" / "_runs" / "events").glob("*/*.jsonl"))
        if latest and latest.get("operation") in {"run", "resume"}:
            # Until a new invocation has a log, do not present the previous
            # invocation's completed tasks as progress of a just-queued resume.
            event_paths = [item for item in event_paths if item.stat().st_mtime >= latest["created_epoch"]]
        if event_paths:
            event_path = max(event_paths, key=lambda item: item.stat().st_mtime_ns)
            for line in tail(event_path, 256_000).splitlines():
                try:
                    events.append(json.loads(line))
                except ValueError:
                    continue  # A writer can be midway through the last event.
        log = tail(self.jobs_dir / latest["id"] / "worker.log") if latest else ""
        return {**row, "summary": read_json(path / "summary.json", {}),
                "saved": read_json(path / "run_config.json", {}),
                "events": events[-300:], "log": log,
                "report": tail(path / "report.md", 2_000_000),
                "report_state": "saved_copy_not_revalidated"}
