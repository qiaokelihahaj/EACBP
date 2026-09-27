"""Command-line adapter for the shared study application service."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Optional, Sequence

from eacbp.application.study_service import (
    run_study, resume_study, report_study, inspect_run, cleanup_artifacts, export_study, import_study,
    load_manifest as _load_manifest, load_config as _load_config,
    preview_study as _planning_preview, _resolve_run_dir,
    normalise_config_paths as _normalise_config_paths,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="eacbp", description="Evidence-aware Agentic Computational Biology Platform")
    parser.add_argument("--version", action="version", version="eacbp 0.1.0")
    commands = parser.add_subparsers(dest="command", required=True)

    webui = commands.add_parser("webui", help="open the local browser interface")
    webui.add_argument("--workspace", type=Path, help="input directory (overrides saved WebUI settings)")
    webui.add_argument("--runs-dir", type=Path, help="results directory (overrides saved WebUI settings)")
    webui.add_argument("--settings-file", type=Path, help="directory preferences file (default: user configuration directory)")
    webui.add_argument("--port", type=int, default=8765)
    webui.add_argument("--no-browser", action="store_true", help="do not open a browser automatically")

    plan = commands.add_parser("plan", help="preview a static study plan")
    plan.add_argument("--manifest", required=True, help="manifest JSON file or JSON object")
    plan.add_argument("--config", help="optional run-config JSON file or JSON object")

    run = commands.add_parser("run", help="import h5ad and execute a new isolated run")
    run.add_argument("--manifest", required=True, help="manifest JSON file or JSON object")
    run.add_argument("--config", help="optional run-config JSON file or JSON object")
    run.add_argument("--data", required=True, type=Path, help="source h5ad file")
    run.add_argument("--run-dir", required=True, type=Path, help="new isolated run directory")

    resume = commands.add_parser("resume", help="resume an existing run from its registry")
    resume.add_argument("--run-dir", required=True, type=Path)

    inspect_parser = commands.add_parser("inspect", help="summarize a run event log")
    inspect_parser.add_argument("--run-dir", type=Path)
    inspect_parser.add_argument("--events", "--event-log", dest="events", type=Path)
    inspect_parser.add_argument("target", nargs="?", type=Path, help="run directory or event JSONL path")

    cleanup = commands.add_parser("cleanup", help="preview or apply artifact transaction cleanup")
    cleanup.add_argument("artifact_root", nargs="?", type=Path)
    cleanup.add_argument("--artifact-root", dest="artifact_root_option", type=Path)
    cleanup.add_argument("--older-than-days", type=int, default=30)
    cleanup.add_argument("--apply", action="store_true", help="apply the cleanup; default is preview")

    report = commands.add_parser("report", help="rebuild Markdown from a verified run snapshot")
    report.add_argument("--run-dir", required=True, type=Path)
    report.add_argument("--output", type=Path)
    export = commands.add_parser("export", help="export verified research artifacts and snapshots")
    export.add_argument("--run-dir", required=True, type=Path)
    export.add_argument("--output", required=True, type=Path)
    restore = commands.add_parser("import", help="validate and import a research bundle into a new directory")
    restore.add_argument("--bundle", required=True, type=Path)
    restore.add_argument("--run-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "webui":
            from eacbp.webui.server import serve
            return serve(workspace=args.workspace, runs_dir=args.runs_dir, port=args.port,
                         open_browser=not args.no_browser, settings_file=args.settings_file)
        if args.command == "plan":
            manifest = _load_manifest(args.manifest)
            result = _planning_preview(manifest, _load_config(args.config))
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("valid", True) else 1
        if args.command == "run":
            manifest = _load_manifest(args.manifest)
            config = _load_config(args.config)
            result = run_study(manifest=manifest, config=config, data=args.data, run_dir=args.run_dir)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("status") == "success" else 1
        if args.command == "resume":
            result = resume_study(run_dir=args.run_dir)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("status") == "success" else 1
        if args.command == "inspect":
            target = args.target
            run_dir = args.run_dir
            events = args.events
            if target is not None:
                if target.suffix.casefold() == ".jsonl":
                    events = target
                else:
                    run_dir = target
            result = inspect_run(run_dir=run_dir, events=events)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0
        if args.command == "cleanup":
            root = args.artifact_root_option or args.artifact_root
            if root is None:
                parser.error("cleanup requires artifact_root or --artifact-root")
            result = cleanup_artifacts(artifact_root=root, older_than_days=args.older_than_days, apply=args.apply)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 1 if result.get("blocked") else 0
        if args.command == "report":
            result = report_study(run_dir=args.run_dir, output=args.output)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0
        if args.command == "export":
            result = export_study(run_dir=args.run_dir, output=args.output)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0
        if args.command == "import":
            result = import_study(bundle=args.bundle, run_dir=args.run_dir)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0
        parser.error(f"unknown command: {args.command}")
    except Exception as exc:
        # A failed run has already written summary/error/journal data.  Keep
        # stderr concise for shell users while returning a conventional
        # non-zero code for automation.
        print(f"eacbp: {type(exc).__name__}: {exc}", file=sys.stderr)
        if args.command in {"run", "resume"}:
            run_path = _resolve_run_dir(args.run_dir)
            if (run_path / "summary.json").is_file():
                try:
                    print((run_path / "summary.json").read_text(encoding="utf-8"))
                except OSError:
                    pass
        return 1


__all__ = [
    "build_parser",
    "cleanup_artifacts",
    "inspect_run",
    "main",
    "report_study",
    "resume_study",
    "run_study",
]
