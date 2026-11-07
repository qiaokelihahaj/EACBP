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

    analyze = commands.add_parser("analyze", help="analyze observed input from a scientific question using native deterministic policies")
    analyze.add_argument("--question", required=True, help="scientific question; no expected genes or workflow configuration")
    analyze.add_argument("--data", required=True, type=Path, help="source h5ad file")
    analyze.add_argument("--species", required=True, help="objective species metadata")
    analyze.add_argument("--tissue", required=True, help="objective tissue metadata")
    analyze.add_argument("--run-dir", required=True, type=Path, help="new isolated run directory")
    analyze.add_argument("--threads", type=int, default=1, help="available CPU threads")
    analyze.add_argument("--max-runtime-hours", type=float, default=12.0, help="resource time constraint")

    resume = commands.add_parser("resume", help="resume an existing run from its registry")
    resume.add_argument("--run-dir", required=True, type=Path)

    recovery_status = commands.add_parser("recovery-status", help="inspect a prepared clustering experiment without changing its journal")
    recovery_status.add_argument("--prepared", required=True, type=Path, help="existing prepared experiment directory")
    recovery_status.add_argument("--bindings-sha256", required=True, help="independently checked frozen bindings SHA256")

    resume_experiment = commands.add_parser("resume-experiment", help="re-admit native results and explicitly resume a prepared clustering experiment")
    resume_experiment.add_argument("--prepared", required=True, type=Path, help="existing prepared experiment directory")
    resume_experiment.add_argument("--bindings-sha256", required=True, help="independently checked frozen bindings SHA256")

    finalize = commands.add_parser("finalize-completed-journal", help="rebuild presentation from an already complete audited native journal; no compute or re-audit")
    finalize.add_argument("--run-dir", required=True, type=Path)
    finalize.add_argument("--compute-image", required=True, help="independently checked original image sha256 digest")
    finalize.add_argument("--source-proof", required=True, type=Path, help="original installed-source/environment proof")
    finalize.add_argument("--source-proof-sha256", required=True, help="independently checked original proof file SHA256")

    verify = commands.add_parser("verify-evidence", help="verify a pinned acceptance index and report engineering/scientific verdicts separately")
    verify.add_argument("--index", required=True, type=Path)
    verify.add_argument("--index-sha256", required=True, help="independently checked acceptance index SHA256")
    verify.add_argument("--source-archive", type=Path, help="explicit historical source mapping if the current checkout has evolved")
    verify.add_argument("--source-archive-sha256", help="independently checked source archive receipt SHA256")
    verify.add_argument("--output", type=Path, help="new JSON verification receipt; existing files are preserved")
    verify.add_argument("--report", type=Path, help="new human-readable verification report")

    descriptive = commands.add_parser("review-descriptive", help="reconcile pinned external descriptive tables without statistical inference")
    descriptive.add_argument("--bundle", required=True, type=Path, help="frozen_scoring_tables_v1 JSON manifest")
    descriptive.add_argument("--bundle-sha256", required=True, help="explicit SHA256 of the reviewed manifest")
    descriptive.add_argument("--output", type=Path, help="new JSON review receipt")
    descriptive.add_argument("--report", type=Path, help="new Markdown report including design limitations")

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
        if args.command == "analyze":
            from eacbp.application.autonomous import analyze_study
            result = analyze_study(question=args.question, data=args.data, species=args.species,
                                   tissue=args.tissue, run_dir=args.run_dir, threads=args.threads,
                                   max_runtime_hours=args.max_runtime_hours)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("status") == "success" and result.get("analysis_request_status") == "complete" else 1
        if args.command == "resume":
            result = resume_study(run_dir=args.run_dir)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("status") == "success" else 1
        if args.command in {"recovery-status", "resume-experiment"}:
            from eacbp.application.recovery import inspect_recovery, resume_recovery
            if args.command == "recovery-status":
                result = inspect_recovery(args.prepared, expected_bindings_sha256=args.bindings_sha256)
                succeeded = result.get("status") == "ready" and result.get("can_resume") is True
            else:
                result = resume_recovery(args.prepared, expected_bindings_sha256=args.bindings_sha256)
                succeeded = result.get("status") == "completed"
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if succeeded else 1
        if args.command == "finalize-completed-journal":
            from eacbp.application.finalization import finalize_completed_journal
            result = finalize_completed_journal(run_dir=args.run_dir, compute_image=args.compute_image,
                                                source_proof=args.source_proof, source_proof_sha256=args.source_proof_sha256)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0 if result.get("status") == "success" else 1
        if args.command == "verify-evidence":
            from eacbp.application.verification import verify_acceptance_index, format_verification_summary
            from eacbp._atomic_files import publish_new_text
            outputs = [path.resolve() for path in (args.output, args.report) if path is not None]
            protected = {path.resolve() for path in (args.index, args.source_archive) if path is not None}
            # Protect declared inputs even when a missing file makes verification
            # fail. A failure receipt must never occupy that evidence location.
            for source in (args.index, args.source_archive):
                if source is None or not source.is_file():
                    continue
                try:
                    declarations = json.loads(source.read_text(encoding="utf-8-sig"))
                except (ValueError, UnicodeError):
                    continue  # The verifier reports invalid input JSON.
                if not isinstance(declarations, dict):
                    continue
                references = []
                entries = declarations.get("evidence_files", [])
                if isinstance(entries, list):
                    references.extend(item.get("path") for item in entries if isinstance(item, dict))
                sources = declarations.get("source_files", {})
                if isinstance(sources, dict):
                    references.extend(sources)
                archives = declarations.get("sources", [])
                if isinstance(archives, list):
                    for item in archives:
                        if isinstance(item, dict):
                            references.extend((item.get("original_path"), item.get("archive_path")))
                protected.update(Path(value).resolve() for value in references if isinstance(value, str))
            if len(outputs) != len(set(outputs)) or any(path in protected or path.exists() for path in outputs):
                raise FileExistsError("Verification outputs must be distinct new files, separate from input evidence")
            result = verify_acceptance_index(args.index, expected_index_sha256=args.index_sha256,
                archive_receipt_path=args.source_archive,
                expected_archive_receipt_sha256=args.source_archive_sha256)
            document = json.dumps(result.model_dump(mode="json"), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
            summary = format_verification_summary(result)
            if args.output is not None:
                publish_new_text(args.output, document)
            if args.report is not None:
                publish_new_text(args.report, summary)
            print(document, end="")
            return 0 if result.verification_status == "passed" else 1
        if args.command == "review-descriptive":
            from eacbp.application.descriptive_review import review_descriptive_bundle, publish_descriptive_review
            result = review_descriptive_bundle(args.bundle, expected_bundle_sha256=args.bundle_sha256)
            publish_descriptive_review(result, output=args.output, report=args.report)
            print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
            return 0
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
        if args.command in {"run", "resume", "analyze"}:
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
