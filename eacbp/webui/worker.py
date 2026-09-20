"""Detached worker: do not bypass the CLI lifecycle or its scientific gates."""
from __future__ import annotations

import os
from pathlib import Path
import sys
import traceback

from .jobs import lease, now, read_json, write_json


def run_job(directory: Path) -> int:
    with lease(directory / "worker.lock"):
        job = read_json(directory / "job.json")
        # A late worker must not execute a job already declared interrupted.
        if job["status"] != "queued":
            return 1
        job.update(status="running", started_at=now(), pid=os.getpid())
        write_json(directory / "job.json", job)
        try:
            with lease(directory.parent.parent / "compute.lock"):
                from eacbp import cli
                from eacbp.schemas.study import StudyManifest
                request = read_json(directory / "request.json")
                print(f"EACBP {request['operation']} → {request['run_dir']}", flush=True)
                if request["operation"] == "run":
                    result = cli.run_study(
                        manifest=StudyManifest.model_validate(request["manifest"]),
                        config=request["config"], data=request["data"], run_dir=request["run_dir"],
                    )
                elif request["operation"] == "resume":
                    result = cli.resume_study(run_dir=request["run_dir"])
                elif request["operation"] == "report":
                    result = cli.report_study(run_dir=request["run_dir"])
                else:
                    raise ValueError("Unknown operation")
                job.update(status="success" if result.get("status") == "success" else "failed",
                           result=result)
                if job["status"] == "failed":
                    job["error"] = "分析有失败或阻断步骤，请查看任务事件和结果摘要。"
                print(f"Finished: {result.get('status')}", flush=True)
        except BaseException as exc:
            job.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            traceback.print_exc()
        finally:
            job["finished_at"] = now()
            write_json(directory / "job.json", job)
        return 0 if job["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(run_job(Path(sys.argv[1]).resolve()))
