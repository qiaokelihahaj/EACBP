"""Replay question-only native analyses in an installed, isolated EACBP image.

No reference directory is exposed to the analysis container. The project names
select external teaching inputs/questions, never analysis policies or answers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess


PROJECTS = {
    "pbmc3k": {
        "file": "pbmc3k_raw.h5ad", "run": "pbmc3k-autonomous-v1", "species": "human", "tissue": "PBMC",
        "question": "Explore this human PBMC dataset, discover the main cell populations and the genes characterizing them.",
    },
    "kang2018": {
        "file": "kang_2018.h5ad", "run": "kang2018-autonomous-v1", "species": "human", "tissue": "PBMC",
        "question": "Compare interferon-stimulated and control human PBMCs within each observed cell population using independent biological donors, and assess donor influence.",
    },
    "pancreas_velocity": {
        "file": "endocrinogenesis_day15.h5ad", "run": "pancreas-autonomous-v1", "species": "mouse", "tissue": "embryonic pancreas",
        "question": "Explore the cell populations and their characteristic genes in this embryonic pancreas dataset, and infer the developmental direction using RNA velocity if the data support it.",
    },
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="data, outputs and caches root; choose F: on this host")
    parser.add_argument("--image", required=True, help="installed EACBP image (no source-code mount)")
    parser.add_argument("--project", choices=PROJECTS, required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--run-name", help="fresh output name for a replay; does not select analysis parameters")
    args = parser.parse_args(argv)
    root = args.root.resolve()
    if root.drive.upper() == "C:":
        parser.error("This curriculum keeps downloaded data and outputs outside C:")
    project = PROJECTS[args.project]
    run_name = args.run_name or project["run"]
    if Path(run_name).name != run_name or run_name in {".", ".."} or any(c in run_name for c in "/\\:"):
        parser.error("run-name must be a single directory name")
    source = root / "datasets" / project["file"]
    if not source.is_file():
        parser.error(f"Missing input: {source}")
    run = root / "runs" / run_name
    if run.exists():
        parser.error(f"Run already exists: {run}; preserve it and use another curriculum root")
    for relative in ("runs", "cache/tmp", "cache/home", "cache/cache", "dependencies/site-packages", "dependencies/linux-site-packages"):
        (root / relative).mkdir(parents=True, exist_ok=True)
    command = [
        "docker", "run", "--rm", "--network", "none", "--read-only", "--user", "10001:10001",
        "--mount", f"type=bind,source={source.as_posix()},target=/input/data.h5ad,readonly",
        "--mount", f"type=bind,source={(root/'runs').as_posix()},target=/results",
        "--mount", f"type=bind,source={(root/'cache').as_posix()},target=/scratch",
        "--mount", f"type=bind,source={(root/'dependencies'/'site-packages').as_posix()},target=/deps,readonly",
        "--mount", f"type=bind,source={(root/'dependencies'/'linux-site-packages').as_posix()},target=/linux-deps,readonly",
        "--env", "PYTHONPATH=/linux-deps:/deps", "--entrypoint", "/opt/venv/bin/eacbp", args.image,
        "analyze", "--question", project["question"], "--data", "/input/data.h5ad",
        "--species", project["species"], "--tissue", project["tissue"],
        "--run-dir", f"/results/{run_name}", "--threads", str(args.threads),
    ]
    receipt = root / "runs" / f"{run_name}-execution-command.json"
    if receipt.exists():
        parser.error(f"Existing execution receipt must be preserved: {receipt}")
    receipt.write_text(json.dumps({"project": args.project, "command": command,
                                  "answer_reference_mounted": False}, indent=2), encoding="utf-8")
    with (root / "runs" / f"{run_name}-execution.log").open("x", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    print(json.dumps({"project": args.project, "returncode": completed.returncode,
                      "run_root": str(run), "command_receipt": str(receipt)}, indent=2))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
