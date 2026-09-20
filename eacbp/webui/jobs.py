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
import subprocess
import sys
import threading
import time
from uuid import uuid4

from eacbp._atomic_json import atomic_write_json


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
                result.append(job)
            except (OSError, ValueError, KeyError, TypeError):
                # A corrupt job record must not permit another computation.
                result.append({"id": path.parent.name, "status": "unknown",
                               "error": "无法读取后台作业记录，请检查 " + str(path),
                               "created_at": "", "run_dir": ""})
        return sorted(result, key=lambda item: item["created_at"], reverse=True)

    def submit(self, operation: str, *, run_dir: Path, manifest=None, config=None, data=None):
        if operation not in {"run", "resume", "report"}:
            raise ValueError("Unknown operation")
        run_dir = contained(self.runs_dir, run_dir)
        if run_dir == self.runs_dir:
            raise ValueError("请选择独立运行目录")
        with self._guard, lease(self.jobs_dir.parent / "submit.lock"):
            if any(job["status"] in {"queued", "running", "unknown"} for job in self.jobs()):
                raise BlockingIOError("已有后台操作正在运行。首版每次只运行一个作业。")
            if operation == "run" and run_dir.exists():
                raise ValueError("运行目录已存在，请创建新运行")
            if operation != "run":
                saved = read_json(run_dir / "run_config.json", {})
                if saved.get("import_completed") is not True:
                    raise ValueError("该运行未完成数据导入，不能恢复或重建报告")
            job_id = uuid4().hex
            directory = self.jobs_dir / job_id
            directory.mkdir()
            job = {"id": job_id, "operation": operation, "status": "queued",
                   "run_dir": str(run_dir), "created_at": now(), "created_epoch": time.time()}
            write_json(directory / "request.json", {
                "operation": operation, "run_dir": str(run_dir), "workspace": str(self.workspace),
                "manifest": manifest, "config": config, "data": str(data) if data else None,
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
        paths = {Path(j["run_dir"]) for j in jobs if j.get("run_dir")}
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
