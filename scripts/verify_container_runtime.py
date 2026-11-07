"""Exercise the built CPU image through a synthetic CLI lifecycle.

The fixture is deliberately synthetic and small. This verifies packaging,
container permissions, and run portability; it does not validate biology.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Sequence


FIXTURE_BUILDER = r'''
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

rng = np.random.default_rng(42)
donors = np.repeat([f"donor_{i}" for i in range(8)], 3)
conditions = np.repeat(["control"] * 4 + ["treated"] * 4, 3)
batches = np.repeat([f"batch_{i % 2}" for i in range(8)], 3)
counts = rng.poisson(5, size=(len(donors), 32)).astype(np.int32)
gene_names = ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r"] + [
    f"synthetic_gene_{i:03d}" for i in range(27)
]
obs = pd.DataFrame(
    {"condition": conditions, "donor": donors, "batch": batches},
    index=[f"cell_{i:03d}" for i in range(len(donors))],
)
var = pd.DataFrame({"gene_name": gene_names}, index=gene_names)
fixture = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
fixture.uns["is_simulated"] = True
fixture.uns["data_origin"] = "synthetic_container_smoke"
fixture.write_h5ad("/fixture/input.h5ad", compression="gzip")
'''


RUNTIME_PROBE = r'''
import importlib
import json
import os
from pathlib import Path
import tempfile

assert os.geteuid() == 10001, f"expected uid 10001, got {os.geteuid()}"
modules = (
    "anndata", "cellrank", "celltypist", "decoupler", "eacbp", "liana",
    "pydeseq2", "scanpy", "skimage",
)
versions = {}
for name in modules:
    module = importlib.import_module(name)
    versions[name] = getattr(module, "__version__", "imported")

for variable in ("HOME", "XDG_CACHE_HOME", "NUMBA_CACHE_DIR", "MPLCONFIGDIR", "TMPDIR"):
    directory = Path(os.environ[variable])
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=directory) as probe:
        probe.write(b"writable by the runtime uid")

try:
    descriptor = os.open("/data/input.h5ad", os.O_WRONLY | os.O_APPEND)
except OSError:
    pass
else:
    os.close(descriptor)
    raise AssertionError("the read-only input mount accepted a write open")

with tempfile.NamedTemporaryFile(dir="/outputs") as probe:
    probe.write(b"writable output mount")
print(json.dumps({"uid": os.geteuid(), "imports": versions, "cache_paths_writable": True,
                  "input_mount_read_only": True, "output_mount_writable": True}, sort_keys=True))
'''


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _make_writable(path: Path) -> None:
    """Best-effort mode preparation for Linux; Windows bind mounts use ACLs."""
    try:
        path.chmod(0o777)
    except OSError:
        pass


def _docker_mount(source: Path, target: str, *, readonly: bool = False) -> str:
    mode = ",readonly" if readonly else ""
    return f"type=bind,source={source.resolve()},target={target}{mode}"


def _run(
    docker: str,
    image: str,
    name: str,
    mounts: Sequence[tuple[Path, str, bool]],
    command: Sequence[str],
    logs: Path,
    *,
    readonly_rootfs: bool = True,
    scratch_tmpfs: bool = True,
) -> subprocess.CompletedProcess[str]:
    args = [docker, "run", "--rm", "--network", "none"]
    if readonly_rootfs:
        args.append("--read-only")
    if scratch_tmpfs:
        args.extend((
            "--tmpfs",
            "/scratch:rw,nosuid,nodev,uid=10001,gid=10001,mode=0700,size=512m",
        ))
    for source, target, readonly in mounts:
        args.extend(("--mount", _docker_mount(source, target, readonly=readonly)))
    args.extend((image, *command))
    completed = subprocess.run(args, capture_output=True, text=True, check=False)
    (logs / f"{name}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (logs / f"{name}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        rendered = subprocess.list2cmdline(args) if os.name == "nt" else " ".join(args)
        raise RuntimeError(
            f"container stage {name!r} failed with exit code {completed.returncode}\n"
            f"command: {rendered}\nstdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    print(f"PASS {name}", flush=True)
    return completed


def _assert_cli_success(result: subprocess.CompletedProcess[str], stage: str) -> dict:
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"{stage} did not return JSON on stdout: {result.stdout!r}") from exc
    if payload.get("status") != "success":
        raise RuntimeError(f"{stage} did not report success: {payload}")
    return payload


def _prepare_layout(root: Path) -> tuple[Path, Path, Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    data = root / "data"
    config = root / "config"
    outputs = root / "outputs"
    logs = root / "logs"
    for directory in (data, config, outputs, logs):
        directory.mkdir()
        _make_writable(directory)

    manifest = {
        "study_id": "container_runtime_smoke",
        "title": "Synthetic container runtime smoke; not a biological analysis",
        "biological_design": {
            "species": "synthetic",
            "tissue": "synthetic",
            "conditions": ["control", "treated"],
            "target_cell_types": [],
        },
        "experimental_design": {
            "biological_unit": "donor",
            "batches": ["batch_0", "batch_1"],
            "total_samples": 8,
            "donor_replicates_per_condition": {"control": 4, "treated": 4},
        },
    }
    config_payload = {
        "method_profile": "baseline",
        "capability_parameters": {"qc": {"min_genes": 1, "max_mito_pct": 100.0}},
    }
    (config / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (config / "config.json").write_text(json.dumps(config_payload, indent=2), encoding="utf-8")
    return data, config, outputs, logs


def verify(image: str, docker: str, root: Path) -> None:
    data, config, outputs, logs = _prepare_layout(root)
    _run(
        docker, image, "fixture-build", [(data, "/fixture", False)],
        ("/opt/venv/bin/python", "-c", FIXTURE_BUILDER), logs,
    )
    input_path = data / "input.h5ad"
    if not input_path.is_file():
        raise RuntimeError("container did not create the synthetic H5AD fixture")
    input_hash = _sha256(input_path)

    mounts = [(data, "/data", True), (config, "/config", True), (outputs, "/outputs", False)]
    probe = _run(
        docker, image, "runtime-permissions-and-imports-default-scratch", mounts,
        ("/opt/venv/bin/python", "-c", RUNTIME_PROBE), logs,
        readonly_rootfs=False, scratch_tmpfs=False,
    )
    probe_payload = json.loads(probe.stdout)
    if probe_payload.get("uid") != 10001 or not probe_payload.get("input_mount_read_only"):
        raise RuntimeError(f"runtime permission probe returned unexpected result: {probe_payload}")

    _run(
        docker, image, "runtime-permissions-readonly-rootfs", mounts,
        ("/opt/venv/bin/python", "-c", RUNTIME_PROBE), logs,
    )

    run_dir = outputs / "run"
    run = _run(
        docker, image, "run", mounts,
        ("/opt/venv/bin/eacbp", "run", "--manifest", "/config/manifest.json",
         "--config", "/config/config.json", "--data", "/data/input.h5ad",
         "--run-dir", "/outputs/run"), logs,
    )
    run_payload = _assert_cli_success(run, "run")
    if run_payload.get("is_simulated") is not True:
        raise RuntimeError("run summary did not preserve the AnnData is_simulated provenance marker")
    if not (run_dir / "snapshot.json").is_file() or not (run_dir / "report.md").is_file():
        raise RuntimeError("run did not persist the expected snapshot and report")

    resume = _run(
        docker, image, "resume", mounts,
        ("/opt/venv/bin/eacbp", "resume", "--run-dir", "/outputs/run"), logs,
    )
    _assert_cli_success(resume, "resume")
    event_logs = list((run_dir / "artifacts" / "_runs" / "events").rglob("*.jsonl"))
    if not event_logs:
        raise RuntimeError("run/resume did not persist an event log")
    event_kinds = {
        json.loads(line)["kind"]
        for event_log in event_logs
        for line in event_log.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }
    if "task_resumed" not in event_kinds:
        raise RuntimeError("resume completed without recording any task_resumed event")

    report = _run(
        docker, image, "report", mounts,
        ("/opt/venv/bin/eacbp", "report", "--run-dir", "/outputs/run",
         "--output", "/outputs/run/rebuilt.md"), logs,
    )
    _assert_cli_success(report, "report")
    if not (run_dir / "rebuilt.md").is_file():
        raise RuntimeError("report command did not create its requested output")

    bundle_path = outputs / "study.eacbp.zip"
    exported = _run(
        docker, image, "export", mounts,
        ("/opt/venv/bin/eacbp", "export", "--run-dir", "/outputs/run",
         "--output", "/outputs/study.eacbp.zip"), logs,
    )
    _assert_cli_success(exported, "export")
    if not bundle_path.is_file() or bundle_path.stat().st_size == 0:
        raise RuntimeError("export command did not create a non-empty bundle")

    imported = _run(
        docker, image, "import", mounts,
        ("/opt/venv/bin/eacbp", "import", "--bundle", "/outputs/study.eacbp.zip",
         "--run-dir", "/outputs/restored"), logs,
    )
    _assert_cli_success(imported, "import")
    if not (outputs / "restored" / "run_config.json").is_file():
        raise RuntimeError("import command did not restore run configuration")

    restored_report = _run(
        docker, image, "restored-report", mounts,
        ("/opt/venv/bin/eacbp", "report", "--run-dir", "/outputs/restored",
         "--output", "/outputs/restored/rebuilt.md"), logs,
    )
    _assert_cli_success(restored_report, "restored report")
    report_text = (outputs / "restored" / "rebuilt.md").read_text(encoding="utf-8")
    if "container_runtime_smoke" not in report_text:
        raise RuntimeError("restored report is missing the synthetic study identifier")
    if _sha256(input_path) != input_hash:
        raise RuntimeError("the synthetic input file changed during the read-only container run")

    print(f"PASS input SHA-256 unchanged: {input_hash}", flush=True)
    print(f"PASS run id: {run_payload.get('run_id', 'recorded')}", flush=True)
    print(f"Evidence directory: {root.resolve()}", flush=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image", default="eacbp:py312-cpu", help="built image tag or digest")
    parser.add_argument("--docker", default="docker", help="Docker CLI executable")
    parser.add_argument(
        "--work-dir", type=Path,
        help="retain synthetic input, run outputs, and per-stage logs in a new/empty directory",
    )
    args = parser.parse_args(argv)
    docker = shutil.which(args.docker)
    if docker is None:
        print(f"Docker CLI not found: {args.docker}", file=sys.stderr)
        return 2

    if args.work_dir:
        root = args.work_dir.expanduser().resolve()
        if root.exists() and any(root.iterdir()):
            print(f"work directory must be new or empty: {root}", file=sys.stderr)
            return 2
        try:
            verify(args.image, docker, root)
        except Exception as exc:
            print(f"container runtime smoke failed: {exc}", file=sys.stderr)
            return 1
        return 0

    try:
        with tempfile.TemporaryDirectory(prefix="eacbp-container-smoke-") as temp_dir:
            verify(args.image, docker, Path(temp_dir))
    except Exception as exc:
        print(f"container runtime smoke failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
