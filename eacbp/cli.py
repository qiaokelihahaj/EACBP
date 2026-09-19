"""Command line interface for reproducible EACBP study runs.

The command line surface is deliberately a thin package boundary around the
same orchestrator, artifact registry, planning, maintenance, and reporting
APIs used by Python callers.  It does not import a data file for ``plan`` or
``resume`` and it never removes a failed run directory.
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager
from copy import deepcopy
from datetime import date, datetime
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Optional, Sequence

from eacbp._atomic_json import atomic_write_json

CLI_SCHEMA_VERSION = 1

# Values in these fields are interpreted by capabilities as filesystem
# resources.  They are made absolute at the start of a run so a later
# ``resume`` launched from another working directory has the same meaning.
_PATH_CONFIG_KEYS = {
    "work_dir",
    "index_path",
    "t2g_path",
    "genome_dir",
    "whitelist_path",
    "gtf_path",
    "checkpoint",
    "checkpoint_path",
    "network_path",
    "resource_path",
    "output_dir",
    "run_cwd",
}


def _normalise_config_paths(value: Any, *, base_dir: Path, parent_key: str = "") -> Any:
    """Resolve configured resource paths against the fresh-run CWD.

    Only fields that are known to represent paths (or nested keys ending in
    ``_path``/``_dir``) are touched.  Command names such as ``STAR`` remain
    command names; URLs and URI-like artifact references are never rewritten.
    """

    if isinstance(value, Mapping):
        if parent_key in {"marker_reference", "terminal_states", "cellrank_terminal_states"}:
            return deepcopy(value)
        checkpoint_base = base_dir
        if isinstance(value.get("run_cwd"), str):
            checkpoint_base = Path(_normalise_config_paths(value["run_cwd"], base_dir=base_dir, parent_key="run_cwd"))
        return {
            str(key): _normalise_config_paths(item, base_dir=checkpoint_base if key in {"extra_args", "cli_args"} else base_dir, parent_key=str(key))
            for key, item in value.items()
        }
    if isinstance(value, list):
        if parent_key.casefold() not in {"extra_args", "cli_args"}:
            return [_normalise_config_paths(item, base_dir=base_dir, parent_key=parent_key) for item in value]
        result = []
        checkpoint_next = False
        for item in value:
            if isinstance(item, str) and item.startswith(("--checkpoint=", "--checkpoint-path=")):
                flag, raw_path = item.split("=", 1)
                resolved = _normalise_config_paths(raw_path, base_dir=base_dir, parent_key="checkpoint")
                result.append(f"{flag}={resolved}")
                checkpoint_next = False
                continue
            if checkpoint_next and isinstance(item, str):
                result.append(_normalise_config_paths(item, base_dir=base_dir, parent_key="checkpoint"))
                checkpoint_next = False
                continue
            result.append(_normalise_config_paths(item, base_dir=base_dir, parent_key=parent_key))
            checkpoint_next = isinstance(item, str) and item in {"--checkpoint", "--checkpoint-path"}
        return result
    if isinstance(value, tuple):
        return tuple(_normalise_config_paths(item, base_dir=base_dir, parent_key=parent_key) for item in value)
    if not isinstance(value, str):
        return value
    key = parent_key.casefold()
    path_like = (
        key in _PATH_CONFIG_KEYS
        or key.endswith("_path")
        or key.endswith("_dir")
        or key == "run_cwd"
        or key == "executable"
        or key == "filepath"
    )
    if not path_like or "://" in value or value.startswith("adata:") or value.startswith("table:"):
        return value
    # A bare executable (``STAR``/``kb``) is resolved by the runtime PATH and
    # must retain that meaning.  An explicit relative/absolute executable
    # path is anchored to the run's original CWD.
    if key == "executable" and "/" not in value and "\\" not in value and not Path(value).is_absolute():
        return value
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = base_dir / candidate
    return str(candidate.resolve(strict=False))


def _json_default(value: Any) -> Any:
    """Serialize common scientific values without hiding unsupported objects."""

    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if hasattr(value, "value") and not isinstance(value, (str, bytes)):
        try:
            return value.value
        except Exception:
            pass
    # numpy and pandas are optional dependencies of the command boundary.
    try:
        import numpy as np

        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, np.ndarray):
            return value.tolist()
    except ImportError:
        pass
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _jsonable(value: Any) -> Any:
    """Return a detached JSON-compatible representation of *value*."""

    if hasattr(value, "model_dump"):
        return value.model_dump(mode="json")
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    try:
        return json.loads(json.dumps(value, default=_json_default, allow_nan=False))
    except (TypeError, ValueError):
        return value


def _read_json_value(value: str | Path, *, label: str) -> Any:
    """Read a JSON file, while also accepting a JSON object on the command line."""

    raw = str(value)
    path = Path(raw).expanduser()
    try:
        is_file = path.is_file()
    except OSError:
        is_file = False
    if is_file:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Unable to read {label} JSON {path}: {exc}") from exc
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"{label} must be a JSON file or JSON value: {exc}") from exc


def _read_json_object(value: str | Path, *, label: str, allow_none: bool = False) -> dict[str, Any]:
    if value is None and allow_none:
        return {}
    payload = _read_json_value(value, label=label)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    """Atomically write one UTF-8 JSON document below an existing run path."""

    atomic_write_json(
        path,
        _jsonable(payload),
        create_parent=True,
        dump_kwargs={"ensure_ascii": False, "indent": 2, "allow_nan": False},
    )


def _hash_file(path: Path) -> str:
    """Return the SHA-256 digest of a source file as lowercase hex."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_manifest(value: str | Path):
    from eacbp.schemas.study import StudyManifest

    payload = _read_json_object(value, label="manifest")
    try:
        return StudyManifest.model_validate(payload)
    except Exception as exc:
        raise ValueError(f"Invalid study manifest: {exc}") from exc


def _load_config(value: str | Path | None) -> dict[str, Any]:
    if value is None:
        return {}
    return _read_json_object(value, label="config")


def _resolve_run_dir(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if path.name in {"", ".", ".."}:
        raise ValueError("run-dir must name a concrete directory")
    return path


def _prepare_new_run_dir(value: str | Path) -> Path:
    path = _resolve_run_dir(value)
    if path.exists():
        raise ValueError(f"run-dir already exists; choose a new isolated directory: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    # Atomic directory creation claims the run path before any metadata or
    # registry is written.  A concurrent invocation cannot share this run.
    path.mkdir(exist_ok=False)
    return path


@contextmanager
def _run_lock(path: Path):
    """Serialize CLI metadata, computation and snapshot publication together."""
    if not path.is_dir():
        raise FileNotFoundError(f"run directory does not exist: {path}")
    from eacbp.artifact.storage import ArtifactStorageBackend
    with ArtifactStorageBackend(str(path)).lock():
        yield


def _manifest_for_input(manifest, data, source_path: Path, source_hash: str):
    """Fill the raw URI and observed spatial flag without changing user intent."""

    result = manifest.model_copy(deep=True)
    if not result.data.raw_artifact_uri:
        result.data.raw_artifact_uri = f"adata://{result.study_id}/raw/v1"
    if not result.data.has_spatial_coordinates and "spatial" in getattr(data, "obsm", {}):
        result.data.has_spatial_coordinates = True
    return result


def _import_h5ad(path: Path):
    """Load h5ad once and convert it to the package's stable SCData payload."""

    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - depends on optional bio extra
        raise RuntimeError("run requires the optional 'anndata' dependency") from exc
    from eacbp.capabilities.sc_data import SCData

    payload = ad.read_h5ad(str(path))
    return SCData.from_anndata(payload)


def _make_registry(path: Path):
    from eacbp.artifact.registry import ArtifactRegistry

    return ArtifactRegistry(str(path))


def _make_orchestrator(registry):
    from eacbp.orchestrator.loop import ScientificOrchestrator

    return ScientificOrchestrator(artifact_registry=registry)


def _source_record(path: Path, *, before: str | None = None, after: str | None = None) -> dict[str, Any]:
    record = {
        "path": str(path.resolve()),
        "sha256": after or before,
        "sha256_before_import": before,
        "sha256_after_import": after,
    }
    return record


def _write_run_config(
    run_dir: Path,
    *,
    manifest,
    config: Mapping[str, Any],
    source: Mapping[str, Any] | None,
    status: str = "created",
    error: str | None = None,
) -> dict[str, Any]:
    """Persist the user inputs and provenance before execution begins."""

    existing: dict[str, Any] = {}
    config_path = run_dir / "run_config.json"
    if config_path.is_file():
        existing = _read_json_object(config_path, label="run_config")
    payload: dict[str, Any] = {
        **existing,
        "schema_version": CLI_SCHEMA_VERSION,
        "status": status,
        "study_id": manifest.study_id,
        "manifest": _jsonable(manifest),
        "config": deepcopy(dict(config)),
        "source": deepcopy(dict(source)) if source is not None else None,
        "registry_dir": str((run_dir / "artifacts").resolve()),
        "import_completed": bool(existing.get("import_completed", False)),
    }
    if source is not None:
        payload["source_path"] = source.get("path")
        payload["source_sha256"] = source.get("sha256")
        payload["source_sha256_before_import"] = source.get("sha256_before_import")
        payload["source_sha256_after_import"] = source.get("sha256_after_import")
    payload.setdefault("run_cwd", existing.get("run_cwd") or str(Path.cwd().resolve()))
    if error:
        payload["error"] = error
    _write_json(run_dir / "run_config.json", payload)
    # These short names are intentionally kept as first-class files: they are
    # convenient for external inspection and preserve a stable package CLI
    # contract independent of the internal run-config envelope.
    _write_json(run_dir / "manifest.json", manifest)
    _write_json(run_dir / "config.json", dict(config))
    if source is not None:
        _write_json(run_dir / "source.json", source)
    return payload


def _update_run_config(run_dir: Path, **updates: Any) -> dict[str, Any]:
    path = run_dir / "run_config.json"
    payload = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    payload.update(_jsonable(updates))
    _write_json(path, payload)
    if "source" in updates and updates["source"] is not None:
        _write_json(run_dir / "source.json", updates["source"])
    return payload


def _summary_failure(run_dir: Path, exc: BaseException) -> dict[str, Any]:
    return {
        "schema_version": CLI_SCHEMA_VERSION,
        "status": "failed",
        "run_dir": str(run_dir),
        "error_type": type(exc).__name__,
        "error": str(exc),
    }


def _persist_snapshot(run_dir: Path, manifest, config: Mapping[str, Any], summary: Mapping[str, Any], orchestrator) -> Path:
    """Persist one typed evidence snapshot through the shared snapshot API."""

    from eacbp.evidence.snapshot import write_study_snapshot

    from uuid import UUID
    run_id = UUID(str(summary["run_id"])).hex
    archive = run_dir / "snapshots" / f"{run_id}.json"
    # Snapshot failures intentionally propagate: a run without its evidence
    # receipt cannot claim complete CLI delivery.
    write_study_snapshot(
        archive,
        manifest=manifest,
        config=dict(config),
        summary=dict(summary),
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=orchestrator.artifact_registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    path = run_dir / "snapshot.json"
    _write_json(path, json.loads(archive.read_text(encoding="utf-8")))
    return path


def _record_failure(run_dir: Path, exc: BaseException, *, source: Mapping[str, Any] | None = None) -> dict[str, Any]:
    if (run_dir / "run_config.json").exists():
        _update_run_config(run_dir, status="failed", error_type=type(exc).__name__, error=str(exc), source=source)
    else:
        _write_json(run_dir / "error.json", _summary_failure(run_dir, exc))
    failure = _summary_failure(run_dir, exc)
    _write_json(run_dir / "summary.json", failure)
    return failure


def _execute_run(
    *,
    run_dir: Path,
    manifest,
    config: Mapping[str, Any],
    source: Mapping[str, Any] | None,
    resume: bool,
    registry=None,
) -> tuple[dict[str, Any], Any]:
    registry = registry if registry is not None else _make_registry(run_dir / "artifacts")
    orchestrator = _make_orchestrator(registry)
    execution_config = deepcopy(dict(config))
    if resume:
        execution_config["resume"] = True
    else:
        execution_config.pop("resume", None)
    summary = orchestrator.run_study(manifest, execution_config)
    summary = dict(_jsonable(summary))
    summary["run_dir"] = str(run_dir)
    summary["report"] = str(run_dir / "report.md")
    summary["snapshot"] = str(run_dir / "snapshot.json")
    _write_json(run_dir / "summary.json", summary)
    _write_json(run_dir / "summaries" / f"{summary.get('run_id') or 'latest'}.json", summary)
    _update_run_config(run_dir, status=summary.get("status", "unknown"), last_run_id=summary.get("run_id"))
    _persist_snapshot(run_dir, manifest, config, summary, orchestrator)
    # Rebuild through the persisted snapshot, which also validates all
    # registry-backed output hashes before generating Markdown.
    report_path = _render_report(run_dir)
    summary["report"] = str(report_path)
    _write_json(run_dir / "summary.json", summary)
    return summary, orchestrator


def _run_study(*, manifest, config: Mapping[str, Any] | None, data: str | Path, run_dir: str | Path) -> dict[str, Any]:
    """Run one fresh study and return its persisted summary.

    This public helper is useful to callers that want the CLI lifecycle but do
    not need argument parsing.  A failed run raises after writing its failure
    record, allowing callers to inspect the exact run directory.
    """

    raw_config = dict(config or {})
    run_cwd = Path.cwd().resolve()
    config = _normalise_config_paths(raw_config, base_dir=run_cwd)
    from eacbp.schemas.runtime import RunConfig

    # Persist the exact typed/validated options used by ExecutionState.  This
    # rejects invalid keys and prevents stringly typed booleans from changing
    # resume semantics.
    config = RunConfig.from_mapping(config).to_mapping()
    run_path = _resolve_run_dir(run_dir)
    source_path = Path(data).expanduser().resolve()
    if not source_path.is_file():
        exc = FileNotFoundError(f"data must be an existing h5ad file: {source_path}")
        _write_json(run_path / "error.json", _summary_failure(run_path, exc))
        _write_json(run_path / "summary.json", _summary_failure(run_path, exc))
        raise exc
    source_before = _hash_file(source_path)
    source = _source_record(source_path, before=source_before)
    manifest = _manifest_for_input(manifest, None, source_path, source_before)
    _write_run_config(run_path, manifest=manifest, config=config, source=source)
    _update_run_config(run_path, config_cwd=str(run_cwd), run_cwd=str(run_cwd))
    try:
        payload = _import_h5ad(source_path)
        source_after = _hash_file(source_path)
        source = _source_record(source_path, before=source_before, after=source_after)
        _update_run_config(run_path, source=source, import_completed=False)
        if source_before != source_after:
            raise ValueError("Input h5ad changed during import")
        manifest = _manifest_for_input(manifest, payload, source_path, source_before)
        _write_run_config(run_path, manifest=manifest, config=config, source=source, status="ready")
        registry = _make_registry(run_path / "artifacts")
        raw_uri = manifest.data.raw_artifact_uri or f"adata://{manifest.study_id}/raw/v1"
        manifest.data.raw_artifact_uri = raw_uri
        registry.register(
            uri_str=raw_uri,
            payload=payload,
            artifact_type=__import__("eacbp.schemas.artifact", fromlist=["ArtifactType"]).ArtifactType.ANNDATA,
            study_id=manifest.study_id,
            created_by_task="task_000_ingest",
            operation="h5ad_import",
            parameters={"source_file": str(source_path), "source_sha256": source_before},
            summary_metrics={"n_obs": payload.n_obs, "n_vars": payload.n_vars},
        )
        _update_run_config(run_path, manifest=manifest, source=source, import_completed=True)
        summary, _ = _execute_run(
            run_dir=run_path, manifest=manifest, config=config, source=source,
            resume=False, registry=registry,
        )
        return summary
    except Exception as exc:
        _record_failure(run_path, exc, source=source)
        raise


def _resume_study(*, run_dir: str | Path) -> dict[str, Any]:
    """Resume a saved run without reading or importing its original h5ad."""

    run_path = _resolve_run_dir(run_dir)
    config_path = run_path / "run_config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"run directory has no saved run_config.json: {run_path}")
    envelope = _read_json_object(config_path, label="run_config")
    if envelope.get("schema_version") != CLI_SCHEMA_VERSION:
        raise ValueError(
            f"unsupported run_config schema_version {envelope.get('schema_version')!r}; "
            f"expected {CLI_SCHEMA_VERSION}"
        )
    if envelope.get("import_completed") is not True:
        raise ValueError("run was not fully imported; resume requires import_completed=true")
    registry_dir = run_path / "artifacts"
    index_path = registry_dir / ".artifact_registry.json"
    if not registry_dir.is_dir() or not index_path.is_file():
        raise ValueError(f"run directory has no persisted artifact registry index: {index_path}")
    manifest_payload = envelope.get("manifest")
    if not isinstance(manifest_payload, dict):
        raise ValueError(f"run_config.json has no manifest object: {run_path}")
    from eacbp.schemas.study import StudyManifest

    manifest = StudyManifest.model_validate(manifest_payload)
    config = envelope.get("config", {})
    if not isinstance(config, dict):
        raise ValueError("saved run configuration must be a JSON object")
    from eacbp.schemas.runtime import RunConfig

    config = RunConfig.from_mapping(config).to_mapping()
    # Constructing ArtifactRegistry verifies the persisted registry index.  No
    # source path is resolved or opened, so changing/removing the original
    # h5ad cannot trigger an implicit re-import.
    registry = _make_registry(registry_dir)
    try:
        summary, _ = _execute_run(
            run_dir=run_path,
            manifest=manifest,
            config=config,
            source=envelope.get("source"),
            resume=True,
            registry=registry,
        )
        return summary
    except Exception as exc:
        _record_failure(run_path, exc, source=envelope.get("source"))
        raise


def _planning_preview(manifest, config: Mapping[str, Any] | None) -> dict[str, Any]:
    from eacbp.orchestrator.planning import preview_study_plan

    # Deliberately leave capability_registry=None.  The planner owns its
    # static registry assembly and no execution, data import, or run directory
    # is created by this code path.
    value = preview_study_plan(manifest, config=config, capability_registry=None)
    result = _jsonable(value)
    if not isinstance(result, dict):
        raise ValueError("preview_study_plan must return a JSON object")
    return result


def _render_report(run_dir: Path, output: str | Path | None = None) -> Path:
    snapshot_path = run_dir / "snapshot.json"
    if not snapshot_path.is_file():
        raise FileNotFoundError(f"run directory has no snapshot.json: {run_dir}")
    from eacbp.evidence.snapshot import render_snapshot_report

    registry_dir = run_dir / "artifacts"
    index_path = registry_dir / ".artifact_registry.json"
    if not registry_dir.is_dir() or not index_path.is_file():
        raise ValueError(f"run directory has no persisted artifact registry index: {index_path}")
    registry = _make_registry(registry_dir)
    # The shared snapshot implementation verifies artifact hashes and audit
    # admission against this registry before rendering the report.
    markdown = render_snapshot_report(snapshot_path, registry)
    if not isinstance(markdown, str):
        raise ValueError("render_snapshot_report must return Markdown text")
    target = Path(output).expanduser().resolve() if output is not None else run_dir / "report.md"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(markdown, encoding="utf-8")
    return target


def run_study(*, manifest, config=None, data, run_dir):
    path = _prepare_new_run_dir(run_dir)
    with _run_lock(path):
        return _run_study(manifest=manifest, config=config, data=data, run_dir=path)


def resume_study(*, run_dir):
    path = _resolve_run_dir(run_dir)
    with _run_lock(path):
        return _resume_study(run_dir=path)


def report_study(*, run_dir: str | Path, output: str | Path | None = None) -> dict[str, Any]:
    run_path = _resolve_run_dir(run_dir)
    with _run_lock(run_path):
        path = _render_report(run_path, output=output)
    return {"status": "success", "run_dir": str(run_path), "report": str(path)}


def _inspect_events(path: Path) -> dict[str, Any]:
    from eacbp.orchestrator.events import read_run_events

    events = read_run_events(path)
    counts = Counter(event.kind for event in events)
    return {
        "status": "success",
        "event_log": str(path.resolve()),
        "run_id": events[0].run_id if events else None,
        "study_id": events[0].study_id if events else None,
        "last_event": events[-1].kind if events else None,
        "elapsed_seconds": events[-1].elapsed_seconds if events else None,
        "event_counts": dict(counts),
        "issues": [
            event.model_dump(mode="json")
            for event in events
            if event.kind in {"task_failed", "task_blocked", "audit_rejected", "run_failed"}
        ],
    }


def inspect_run(*, run_dir: str | Path | None = None, events: str | Path | None = None) -> dict[str, Any]:
    if events is not None:
        return _inspect_events(Path(events).expanduser().resolve())
    if run_dir is None:
        raise ValueError("inspect requires --run-dir or --events")
    path = _resolve_run_dir(run_dir)
    summary_path = path / "summary.json"
    event_path: Path | None = None
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            candidate = summary.get("event_log")
            if isinstance(candidate, str) and Path(candidate).is_file():
                event_path = Path(candidate)
        except (OSError, json.JSONDecodeError):
            pass
    if event_path is None:
        candidates = sorted(path.rglob("*.jsonl"), key=lambda item: item.stat().st_mtime_ns)
        if candidates:
            event_path = candidates[-1]
    if event_path is None:
        raise FileNotFoundError(f"no event JSONL log found below {path}")
    return _inspect_events(event_path)


def cleanup_artifacts(*, artifact_root: str | Path, older_than_days: int = 30, apply: bool = False) -> dict[str, Any]:
    from eacbp.artifact.maintenance import cleanup_artifacts as _cleanup

    report = _cleanup(artifact_root, apply=apply, older_than_days=older_than_days)
    result = _jsonable(report.as_dict() if hasattr(report, "as_dict") else report)
    if isinstance(result, dict):
        result["blocked"] = bool(getattr(report, "blocked", False))
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="eacbp", description="Evidence-aware Agentic Computational Biology Platform")
    parser.add_argument("--version", action="version", version="eacbp 0.1.0")
    commands = parser.add_subparsers(dest="command", required=True)

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
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
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
