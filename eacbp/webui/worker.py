"""Detached worker: do not bypass the CLI lifecycle or its scientific gates."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import traceback

from .jobs import (RECOVERY_OPERATION, contained, lease, now, read_json, recovery_directory,
                   recovery_sha, study_preview_binding, write_json)


def recovery_error(result):
    """Keep the executor's durable reasons rather than flattening partial work."""
    reasons = []
    for item in result.get("blockers", []):
        if isinstance(item, dict):
            reasons.append(str(item.get("reason") or item.get("code") or item))
    for key in ("error", "reason", "failures", "native_failures", "postprocessing_failures"):
        if result.get(key):
            reasons.append(str(result[key]))
    return "；".join(reasons) or "实验恢复未完成，请查看持久阶段记录与结果中的失败原因。"


def recovery_failure_details(prepared, result):
    """Read bounded saved diagnostic reasons; never grant scientific admission."""
    if result.get("status") not in {"partial", "failed"} or not result.get("comparison"):
        return {}
    try:
        path = contained(prepared, result["comparison"])
        if path.stat().st_size > 8_000_000:
            raise ValueError("比较诊断超过 8 MB 显示限制")
        with path.open("rb") as stream:
            data = stream.read(8_000_001)
        if len(data) > 8_000_000:
            raise ValueError("比较诊断超过 8 MB 显示限制")
        comparison = json.loads(data.decode("utf-8"))
        if not isinstance(comparison, dict):
            raise ValueError("比较诊断必须是 JSON 对象")
        return {key: comparison[key] for key in ("reason", "failures", "primary_native_summary")
                if key in comparison}
    except (OSError, ValueError, TypeError) as exc:
        return {"reason": f"比较失败详情暂无法读取：{type(exc).__name__}: {exc}"}


def run_job(directory: Path) -> int:
    with lease(directory / "worker.lock"):
        job = read_json(directory / "job.json")
        # A late worker must not execute a job already declared interrupted.
        if job["status"] != "queued":
            return 1
        job.update(status="running", started_at=now(), pid=os.getpid())
        write_json(directory / "job.json", job)
        request = None
        request_validated = False
        try:
            with lease(directory.parent.parent / "compute.lock"):
                from eacbp.application import study_service
                from eacbp.schemas.study import StudyManifest
                request = read_json(directory / "request.json")
                if request.get("operation") != job.get("operation"):
                    raise ValueError("后台作业的操作绑定已改变")
                if request.get("run_dir") != job.get("run_dir"):
                    raise ValueError("后台作业的运行目录绑定已改变")
                request_validated = True
                print(f"EACBP {request['operation']} → {request['run_dir']}", flush=True)
                if request["operation"] == "run":
                    from eacbp.orchestrator.planning import preview_study_plan

                    binding = request.get("preview_binding")
                    if not isinstance(binding, dict) or binding != job.get("preview_binding"):
                        raise ValueError("后台作业的方案确认绑定缺失或已改变，请重新检查方案。")
                    manifest = StudyManifest.model_validate(request["manifest"])
                    current = study_preview_binding(data=request["data"],
                        source_sha256=binding.get("source_sha256"),
                        manifest=manifest.model_dump(mode="json"), config=request["config"],
                        plan=preview_study_plan(manifest, request["config"]))
                    if binding != current:
                        raise ValueError("后台作业的输入路径、配置或方案已改变，请重新检查方案。")
                    result = study_service.run_study(
                        manifest=manifest,
                        config=request["config"], data=request["data"], run_dir=request["run_dir"],
                        expected_source_sha256=binding["source_sha256"],
                    )
                elif request["operation"] == "resume":
                    result = study_service.resume_study(run_dir=request["run_dir"])
                elif request["operation"] == "report":
                    result = study_service.report_study(run_dir=request["run_dir"])
                elif request["operation"] == RECOVERY_OPERATION:
                    from eacbp.application.recovery import resume_recovery

                    # Bind to this job queue's approved root, not a root chosen
                    # by a modified request. Resolve again after queueing.
                    runs_root = directory.resolve().parents[2]
                    if Path(request["approved_runs_root"]).resolve() != runs_root:
                        raise ValueError("恢复作业的批准结果目录已改变")
                    prepared = recovery_directory(runs_root, request["run_dir"])
                    if job["run_dir"] != request["run_dir"] or (
                        job.get("expected_bindings_sha256") != request.get("expected_bindings_sha256")
                    ):
                        raise ValueError("恢复作业的路径或 SHA 绑定已改变")
                    digest = recovery_sha(request["expected_bindings_sha256"])
                    result = resume_recovery(prepared, expected_bindings_sha256=digest)
                else:
                    raise ValueError("Unknown operation")
                recovery = request["operation"] == RECOVERY_OPERATION
                expected_status = "completed" if recovery else "success"
                job.update(status="success" if result.get("status") == expected_status else "failed",
                           result=result)
                if recovery:
                    job.update(recovery_status=result.get("status"),
                               biological_stability_claim_supported=False)
                    job["recovery_failure_details"] = recovery_failure_details(prepared, result)
                if job["status"] == "failed":
                    job["error"] = recovery_error({**result, **job["recovery_failure_details"]}) if recovery else "分析有失败或阻断步骤，请查看任务事件和结果摘要。"
                print(f"Finished: {result.get('status')}", flush=True)
        except BaseException as exc:
            job.update(status="failed", error=f"{type(exc).__name__}: {exc}")
            traceback.print_exc()
        finally:
            if request_validated and job.get("operation") == "run" and request.get("webui_context"):
                try:
                    run_dir = Path(request["run_dir"])
                    if run_dir.is_dir():
                        write_json(run_dir / "webui_context.json", request["webui_context"])
                except Exception as exc:
                    print(f"Unable to save optional WebUI context: {type(exc).__name__}: {exc}", flush=True)
            job["finished_at"] = now()
            write_json(directory / "job.json", job)
        return 0 if job["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(run_job(Path(sys.argv[1]).resolve()))
