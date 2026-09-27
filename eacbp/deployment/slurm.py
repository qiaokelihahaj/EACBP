"""Thin Slurm controls for one EACBP CLI invocation per scheduled job.

This module only wraps Slurm commands. Study execution and recovery remain the
responsibility of ``eacbp run``, ``eacbp resume``, and ``eacbp report`` inside
the container.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sysconfig
import sys
from typing import Mapping, Sequence


_RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_JOB_ID_RE = re.compile(r"^[0-9]+(?:_[0-9]+)?$")
_MEMORY_RE = re.compile(r"^[1-9][0-9]*[KMGTP]?$", re.IGNORECASE)
_TIME_RE = re.compile(r"^(?:[0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$|^[0-9]+:[0-9]{2}$")
_TOKEN_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


class DeploymentError(ValueError):
    """A deployment argument or Slurm response is invalid."""


def default_batch_script() -> Path:
    """Locate the one canonical batch template in source and wheel installs."""
    relative = Path("share") / "eacbp" / "slurm" / "run_eacbp_container.sbatch"
    installed = Path(sysconfig.get_path("data")) / relative
    if installed.is_file():
        return installed
    source_tree = Path(__file__).resolve().parents[2] / "slurm" / "run_eacbp_container.sbatch"
    return source_tree if source_tree.is_file() else installed


def _required_path(env: Mapping[str, str], name: str, *, directory: bool = False) -> Path:
    value = env.get(name, "").strip()
    if not value:
        raise DeploymentError(f"Set {name} before submitting the Slurm job")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise DeploymentError(f"{name} must be an absolute host path")
    try:
        resolved = path.resolve(strict=True)
    except OSError as exc:
        raise DeploymentError(f"{name} does not exist: {path}") from exc
    expected = resolved.is_dir() if directory else resolved.is_file()
    if not expected:
        kind = "directory" if directory else "file"
        raise DeploymentError(f"{name} must name an existing {kind}: {resolved}")
    return resolved


def _file_under(root: Path, relative: str, name: str) -> Path:
    rel = Path(relative)
    if not relative.strip() or rel.is_absolute() or ".." in rel.parts:
        raise DeploymentError(f"{name} must be a relative path contained by its input directory")
    try:
        resolved = (root / rel).resolve(strict=True)
    except OSError as exc:
        raise DeploymentError(f"{name} does not exist under {root}: {relative}") from exc
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise DeploymentError(f"{name} must resolve to a file inside {root}")
    return resolved


def validate_job_environment(env: Mapping[str, str] | None = None) -> dict[str, str]:
    """Validate and normalize the host paths needed by the batch template."""
    values = dict(os.environ if env is None else env)
    action = values.get("EACBP_ACTION", "").strip().lower()
    if action not in {"run", "resume", "report"}:
        raise DeploymentError("EACBP_ACTION must be run, resume, or report")
    run_id = values.get("EACBP_RUN_ID", "").strip()
    if not _RUN_ID_RE.fullmatch(run_id) or run_id in {".", ".."}:
        raise DeploymentError("EACBP_RUN_ID must be one safe path component using letters, digits, '.', '_' or '-'")

    image = _required_path(values, "EACBP_IMAGE")
    if image.suffix.lower() != ".sif":
        raise DeploymentError("EACBP_IMAGE must point to a prepared .sif file")
    output_root = _required_path(values, "EACBP_OUTPUTS", directory=True)
    if not os.access(output_root, os.W_OK | os.X_OK):
        raise DeploymentError(f"EACBP_OUTPUTS is not writable: {output_root}")
    run_dir = output_root / run_id
    if action == "run" and run_dir.exists():
        raise DeploymentError(f"Run directory already exists; choose a new EACBP_RUN_ID: {run_dir}")
    if action in {"resume", "report"} and not run_dir.is_dir():
        raise DeploymentError(f"EACBP_ACTION={action} requires an existing run directory: {run_dir}")

    values.update({
        "EACBP_ACTION": action,
        "EACBP_RUN_ID": run_id,
        "EACBP_IMAGE": str(image),
        "EACBP_OUTPUTS": str(output_root),
    })
    expected_sha = values.get("EACBP_IMAGE_SHA256", "").strip().lower()
    if expected_sha:
        if not re.fullmatch(r"[0-9a-f]{64}", expected_sha):
            raise DeploymentError("EACBP_IMAGE_SHA256 must contain 64 hexadecimal characters")
        values["EACBP_IMAGE_SHA256"] = expected_sha

    if action == "run":
        data_root = _required_path(values, "EACBP_DATA_DIR", directory=True)
        config_root = _required_path(values, "EACBP_CONFIG_DIR", directory=True)
        data_file = _file_under(data_root, values.get("EACBP_DATA", ""), "EACBP_DATA")
        manifest = _file_under(config_root, values.get("EACBP_MANIFEST", ""), "EACBP_MANIFEST")
        values.update({
            "EACBP_DATA_DIR": str(data_root),
            "EACBP_DATA": str(data_file.relative_to(data_root)),
            "EACBP_CONFIG_DIR": str(config_root),
            "EACBP_MANIFEST": str(manifest.relative_to(config_root)),
        })
        config = values.get("EACBP_CONFIG", "").strip()
        if config:
            config_file = _file_under(config_root, config, "EACBP_CONFIG")
            values["EACBP_CONFIG"] = str(config_file.relative_to(config_root))
    else:
        # Runs persist absolute resource paths. A resume/report job may need
        # those same container mount points even though it does not re-import
        # the source h5ad or reload the original manifest/config files.
        for name in ("EACBP_DATA_DIR", "EACBP_CONFIG_DIR", "EACBP_REFS_DIR"):
            if values.get(name, "").strip():
                values[name] = str(_required_path(values, name, directory=True))
    if action == "run" and values.get("EACBP_REFS_DIR", "").strip():
        values["EACBP_REFS_DIR"] = str(_required_path(values, "EACBP_REFS_DIR", directory=True))
    return values


def _slurm_value(name: str, value: str, pattern: re.Pattern[str]) -> str:
    if not pattern.fullmatch(value):
        raise DeploymentError(f"Invalid {name}: {value!r}")
    return value


def build_submit_command(
    *,
    script: Path,
    logs_dir: Path,
    cwd: Path,
    cpus_per_task: int = 4,
    memory: str = "16G",
    time_limit: str = "04:00:00",
    partition: str | None = None,
    account: str | None = None,
) -> list[str]:
    """Build an argv-only ``sbatch`` command; values are never shell-expanded."""
    if cpus_per_task < 1:
        raise DeploymentError("cpus_per_task must be at least 1")
    memory = _slurm_value("memory", memory, _MEMORY_RE)
    time_limit = _slurm_value("time limit", time_limit, _TIME_RE)
    command = [
        "sbatch", "--parsable", "--export=ALL", "--job-name=eacbp",
        f"--chdir={cwd}",
        f"--output={logs_dir}/%x-%j.out",
        f"--error={logs_dir}/%x-%j.err",
        f"--cpus-per-task={cpus_per_task}", f"--mem={memory}", f"--time={time_limit}",
    ]
    for name, value in (("partition", partition), ("account", account)):
        if value:
            command.append(f"--{name}={_slurm_value(name, value, _TOKEN_RE)}")
    command.append(str(script))
    return command


def _parse_job_id(raw: str) -> str:
    job_id = raw.strip().split(";", 1)[0]
    if not _JOB_ID_RE.fullmatch(job_id):
        raise DeploymentError(f"Unexpected job id returned by Slurm: {raw.strip()!r}")
    return job_id


def _require_command(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise DeploymentError(f"Required Slurm command is unavailable: {name}")
    return executable


def _run(
    command: Sequence[str], *, env: Mapping[str, str] | None = None
) -> subprocess.CompletedProcess[str]:
    try:
        return subprocess.run(command, check=False, text=True, capture_output=True, env=env)
    except OSError as exc:
        raise DeploymentError(f"Could not run {command[0]}: {exc}") from exc


def _submit(args: argparse.Namespace) -> int:
    values = validate_job_environment()
    script = args.script.expanduser().resolve(strict=True)
    if not script.is_file():
        raise DeploymentError(f"Slurm batch script is not a file: {script}")
    logs_dir = args.logs_dir.expanduser().resolve()
    cwd = args.cwd.expanduser().resolve(strict=True)
    if not cwd.is_dir():
        raise DeploymentError(f"Submit working directory is not a directory: {cwd}")
    command = build_submit_command(
        script=script, logs_dir=logs_dir, cwd=cwd,
        cpus_per_task=args.cpus_per_task, memory=args.mem, time_limit=args.time,
        partition=args.partition, account=args.account,
    )
    if args.dry_run:
        print(json.dumps({
            "mode": "dry-run",
            "action": values["EACBP_ACTION"],
            "run_id": values["EACBP_RUN_ID"],
            "command": command,
            "stdout_log": str(logs_dir / "eacbp-%j.out"),
            "stderr_log": str(logs_dir / "eacbp-%j.err"),
            "scheduler_contacted": False,
        }, ensure_ascii=False, indent=2))
        return 0

    _require_command("sbatch")
    logs_dir.mkdir(parents=True, exist_ok=True)
    if not os.access(logs_dir, os.W_OK | os.X_OK):
        raise DeploymentError(f"Slurm log directory is not writable: {logs_dir}")
    completed = _run(command, env=values)
    if completed.returncode:
        if completed.stderr:
            print(completed.stderr.rstrip(), file=sys.stderr)
        return completed.returncode
    job_id = _parse_job_id(completed.stdout)
    print(json.dumps({
        "job_id": job_id,
        "action": values["EACBP_ACTION"],
        "run_id": values["EACBP_RUN_ID"],
        "stdout_log": str(logs_dir / f"eacbp-{job_id}.out"),
        "stderr_log": str(logs_dir / f"eacbp-{job_id}.err"),
    }, ensure_ascii=False, indent=2))
    return 0


def _validate_job_id(raw: str) -> str:
    value = raw.strip()
    if not _JOB_ID_RE.fullmatch(value):
        raise DeploymentError("JOB_ID must be a numeric Slurm job id, optionally with an array index")
    return value


def _status(args: argparse.Namespace) -> int:
    job_id = _validate_job_id(args.job_id)
    squeue = _require_command("squeue")
    live = _run([squeue, "--noheader", "--jobs", job_id, "--format=%T"])
    if live.returncode:
        if live.stderr:
            print(live.stderr.rstrip(), file=sys.stderr)
        return live.returncode
    state = live.stdout.strip()
    source = "squeue"
    if not state:
        sacct = _require_command("sacct")
        finished = _run([sacct, "--noheader", "--allocations", "--jobs", job_id, "--format=State", "--parsable2"])
        if finished.returncode:
            if finished.stderr:
                print(finished.stderr.rstrip(), file=sys.stderr)
            return finished.returncode
        state = next((line.split("|", 1)[0].strip() for line in finished.stdout.splitlines() if line.strip()), "UNKNOWN")
        source = "sacct"
    print(json.dumps({"job_id": job_id, "state": state, "source": source}, ensure_ascii=False))
    return 0


def _cancel(args: argparse.Namespace) -> int:
    job_id = _validate_job_id(args.job_id)
    completed = _run([_require_command("scancel"), job_id])
    if completed.returncode:
        if completed.stderr:
            print(completed.stderr.rstrip(), file=sys.stderr)
        return completed.returncode
    print(json.dumps({"job_id": job_id, "cancel_requested": True}, ensure_ascii=False))
    return 0


def _tail(path: Path, line_count: int) -> str:
    if not path.is_file():
        return "(log file is not available yet)"
    with path.open("rb") as stream:
        stream.seek(0, os.SEEK_END)
        position = stream.tell()
        block_size = 8192
        data = b""
        while position > 0 and data.count(b"\n") <= line_count:
            size = min(block_size, position)
            position -= size
            stream.seek(position)
            data = stream.read(size) + data
    lines = data.splitlines(keepends=True)
    return b"".join(lines[-line_count:]).decode("utf-8", errors="replace").rstrip()


def _logs(args: argparse.Namespace) -> int:
    job_id = _validate_job_id(args.job_id)
    if args.lines < 1:
        raise DeploymentError("--lines must be at least 1")
    logs_dir = args.logs_dir.expanduser().resolve(strict=True)
    for stream in ("out", "err"):
        log_path = logs_dir / f"eacbp-{job_id}.{stream}"
        print(f"--- {stream}: {log_path} ---")
        print(_tail(log_path, args.lines))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m eacbp.deployment", description="Thin Slurm controls for EACBP container jobs")
    commands = parser.add_subparsers(dest="command", required=True)
    submit = commands.add_parser("submit", help="submit one EACBP CLI invocation")
    submit.add_argument("--script", type=Path, default=default_batch_script())
    submit.add_argument("--logs-dir", type=Path, required=True)
    submit.add_argument("--cwd", type=Path, default=Path.cwd())
    submit.add_argument("--cpus-per-task", type=int, default=4)
    submit.add_argument("--mem", default="16G")
    submit.add_argument("--time", default="04:00:00")
    submit.add_argument("--partition")
    submit.add_argument("--account")
    submit.add_argument("--dry-run", action="store_true", help="validate inputs and print the planned sbatch command without contacting Slurm")
    submit.set_defaults(handler=_submit)

    for name, handler in (("status", _status), ("cancel", _cancel)):
        command = commands.add_parser(name, help=f"{name} one Slurm job")
        command.add_argument("job_id")
        command.set_defaults(handler=handler)
    logs = commands.add_parser("logs", help="show the last lines of stdout and stderr")
    logs.add_argument("job_id")
    logs.add_argument("--logs-dir", type=Path, required=True)
    logs.add_argument("--lines", type=int, default=100)
    logs.set_defaults(handler=_logs)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.handler(args)
    except (DeploymentError, OSError) as exc:
        print(f"eacbp deployment: {exc}", file=sys.stderr)
        return 2


__all__ = ["DeploymentError", "build_parser", "build_submit_command", "main", "validate_job_environment"]
