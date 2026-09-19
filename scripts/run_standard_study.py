"""Run a standard study through the package CLI lifecycle.

This script remains as a compatibility entry point for existing notebooks and
batch jobs. Importing, artifact registration, orchestration, snapshotting, and
report generation all belong to :mod:`eacbp.cli`; keeping this adapter small
prevents the two entry points from drifting apart.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys
from uuid import uuid4

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from eacbp.cli import run_study
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest


def _load_analysis_config(path: Path, parser: argparse.ArgumentParser) -> dict:
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        parser.error(f"Unable to read analysis-config JSON: {exc}")
    if not isinstance(state, dict):
        parser.error("analysis-config must be a JSON object")
    supported = {
        "capability_parameters",
        "target_parameters",
        "method_overrides",
        "analysis_extensions",
        "advanced_analysis",
    }
    unknown = set(state) - supported
    if unknown:
        parser.error(f"Unsupported analysis-config keys: {sorted(unknown)}")
    return state


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--study-id", default="standard_study")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "outputs" / "runs")
    parser.add_argument("--target-cell-type", action="append", default=[])
    parser.add_argument("--root-cell-id")
    parser.add_argument("--marker-reference", type=Path, help="JSON mapping cell types to marker gene names")
    parser.add_argument("--terminal-states", type=Path, help="JSON mapping CellRank fates to terminal cell IDs")
    parser.add_argument("--condition-a")
    parser.add_argument("--condition-b")
    parser.add_argument("--paga", action="store_true")
    parser.add_argument(
        "--analysis-config",
        type=Path,
        help="JSON analysis configuration: capability_parameters, method_overrides and optional analysis_extensions",
    )
    parser.add_argument(
        "--advanced",
        action="store_true",
        help="Use PyDESeq2 and donor leave-one-out sensitivity; requires advanced-statistics dependencies",
    )
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", args.study_id):
        parser.error("study-id must contain only letters, numbers, underscores or hyphens")
    if not args.data.is_file():
        parser.error("data must be an existing h5ad file")

    state = _load_analysis_config(args.analysis_config, parser) if args.analysis_config else {}
    parameters = state.setdefault("capability_parameters", {})
    if not isinstance(parameters, dict):
        parser.error("capability_parameters must be an object")

    # Command-line values retain precedence over a config file, matching the
    # historical script. A configured root is also a valid source for the
    # terminal-state compatibility option.
    trajectory = parameters.setdefault("trajectory_inference", {})
    if not isinstance(trajectory, dict):
        parser.error("capability_parameters.trajectory_inference must be an object")
    if args.root_cell_id:
        trajectory["root_cell_id"] = args.root_cell_id
    if args.paga:
        trajectory["run_paga"] = True

    if args.marker_reference:
        try:
            marker_reference = json.loads(args.marker_reference.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"Unable to read marker-reference JSON: {exc}")
        # Preserve the old option's replacement semantics: an explicit CLI
        # marker map replaces the complete clustering parameter block.
        parameters["clustering"] = {"marker_reference": marker_reference}

    if args.condition_a is not None or args.condition_b is not None:
        for name in ("deg", "differential_abundance"):
            condition_parameters = parameters.setdefault(name, {})
            if not isinstance(condition_parameters, dict):
                parser.error(f"capability_parameters.{name} must be an object")
            condition_parameters.update(condition_a=args.condition_a, condition_b=args.condition_b)

    state["method_profile"] = "standard"
    if args.advanced:
        state["advanced_analysis"] = True

    if args.terminal_states:
        try:
            terminal_states = json.loads(args.terminal_states.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            parser.error(f"Unable to read terminal-states JSON: {exc}")
        if args.root_cell_id is None and not trajectory.get("root_cell_id"):
            parser.error("terminal-states requires root-cell-id (on the command line or in trajectory_inference)")
        state["cellrank_terminal_states"] = terminal_states

    manifest = StudyManifest(
        study_id=args.study_id,
        biological_design=BiologicalDesign(
            species=args.species,
            tissue=args.tissue,
            target_cell_types=args.target_cell_type,
        ),
        data=DataSpec(),
    )
    run_dir = args.output_dir.resolve() / args.study_id / uuid4().hex
    summary = run_study(manifest=manifest, config=state, data=args.data, run_dir=run_dir)
    print(f"Study status: {summary['status']}\nReport: {summary['report']}")
    return 0 if summary["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
