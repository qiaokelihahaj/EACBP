#!/usr/bin/env python3
"""Run a real-data acceptance check on the public Kang IFN-stimulation study.

The source H5AD is never modified. The script selects every cell of one exact
cell-type label, retains all genes, copies the raw X matrix to ``layers/counts``,
and writes every output beneath the required external work directory.
"""

from __future__ import annotations

import argparse
from datetime import date, datetime
import hashlib
import importlib.metadata
import json
import math
import platform
import re
import sys
import tempfile
from pathlib import Path, PurePath, PureWindowsPath
from typing import Any

SOURCE_URL = "https://exampledata.scverse.org/pertpy/kang_2018.h5ad"
SOURCE_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"
SOURCE_CITATION = (
    "Kang et al. (2018), Multiplexed droplet single-cell RNA-sequencing using "
    "natural genetic variation, Nat Biotechnol 36:89-94, DOI:10.1038/nbt.4042; "
    "GEO:GSE96583; distributed as Figshare article 19397624, version 2."
)
SOURCE_LICENSE = "CC BY 4.0 (Figshare article 19397624, version 2)"
MARKER_REFERENCE_URL = "https://satijalab.org/seurat/articles/integration_introduction.html"
MARKER_GENES = ("ISG15", "IFI6", "LY6E", "MX1", "IFIT1", "IFIT2", "ISG20")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
LOCKFILE_DEFAULT = (
    Path(__file__).resolve().parents[1]
    / "docker"
    / "requirements-linux-py312-scientific.lock"
)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_external_path(work_dir: PurePath, *, project_root: PurePath) -> PurePath:
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute path on the external validation volume")
    if isinstance(work_dir, PureWindowsPath) != isinstance(project_root, PureWindowsPath):
        raise ValueError("--work-dir and repository root must use the same path style")
    try:
        work_dir.relative_to(project_root)
    except ValueError:
        pass
    else:
        raise ValueError("--work-dir must be outside the repository/workspace so generated files stay on the validation volume")
    if isinstance(work_dir, PureWindowsPath) and work_dir.drive.casefold() == "c:":
        raise ValueError("--work-dir must not use the C: drive; choose an external validation volume such as F:")
    return work_dir


def _validate_external_work_dir(work_dir: Path, *, project_root: Path = PROJECT_ROOT) -> Path:
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute path on the external validation volume")
    resolved = work_dir.expanduser().resolve()
    _validate_external_path(resolved, project_root=project_root.resolve())
    return resolved


def _source_integrity_report(path: Path, sha256_before: str) -> dict[str, Any]:
    try:
        sha256_after = _sha256_file(path)
    except OSError as exc:
        return {
            "path": str(path),
            "sha256_before": sha256_before,
            "sha256_after": None,
            "unchanged": False,
            "integrity_error": f"{type(exc).__name__}: {exc}",
        }
    return {
        "path": str(path),
        "sha256_before": sha256_before,
        "sha256_after": sha256_after,
        "unchanged": sha256_before == sha256_after,
    }


def _apply_source_integrity(summary: dict[str, Any], path: Path, sha256_before: str) -> dict[str, Any]:
    integrity = _source_integrity_report(path, sha256_before)
    source_summary = summary.setdefault("source", {})
    source_summary.update(integrity)
    if not integrity["unchanged"]:
        summary["overall_status"] = "failed"
        source_summary["integrity_failure"] = (
            "The original source H5AD changed during the run."
            if integrity["sha256_after"] is not None
            else "The original source H5AD could not be re-read for the final SHA256 check."
        )
    return integrity


def _parse_lock_packages(text: str) -> dict[str, str]:
    """Read exact ``name==version`` pins from the pip-compatible uv lock."""
    package_pattern = re.compile(r"^([A-Za-z0-9][A-Za-z0-9._-]*)==([^\s\\]+)\s*\\?\s*$")
    packages: dict[str, str] = {}
    normalized: dict[str, str] = {}
    for line_number, line in enumerate(text.splitlines(), start=1):
        match = package_pattern.fullmatch(line.strip())
        if not match:
            continue
        name, version = match.groups()
        key = re.sub(r"[-_.]+", "-", name).lower()
        if key in normalized:
            raise ValueError(f"Duplicate distribution pin in lockfile at line {line_number}: {name}")
        normalized[key] = name
        packages[name] = version
    if not packages:
        raise ValueError("No exact name==version pins were found in the scientific lockfile.")
    return packages


def _read_lock_packages(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Scientific lockfile does not exist: {path}")
    return _parse_lock_packages(path.read_text(encoding="utf-8"))


def _validate_raw_counts(matrix: Any) -> None:
    """Reject transformed or invalid values without densifying sparse counts."""
    import numpy as np
    from scipy import sparse

    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).reshape(-1)
    if not np.issubdtype(values.dtype, np.integer):
        if not np.isfinite(values).all() or not np.equal(values, np.floor(values)).all():
            raise ValueError("Source X contains non-finite or non-integer values; raw counts are required.")
    if (values < 0).any():
        raise ValueError("Source X contains negative values; non-negative raw counts are required.")


def _safe_slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
    if not slug:
        raise ValueError(f"Cell type label cannot be converted to a study identifier: {value!r}")
    return slug


def _json_safe(value: Any) -> Any:
    import numpy as np

    if hasattr(value, "model_dump"):
        return _json_safe(value.model_dump(mode="json"))
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _write_json(payload: dict[str, Any], path: Path) -> None:
    path.write_text(
        json.dumps(_json_safe(payload), indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _marker_direction_report(table: Any) -> dict[str, Any]:
    """Check pre-registered sign only; statistical significance is not implied."""
    import pandas as pd

    if "gene" not in table.columns or "log2_fold_change" not in table.columns:
        return {
            "status": "failed",
            "criterion": "finite log2_fold_change > 0 for every pre-registered marker",
            "reason": "PyDESeq2 output lacks gene or log2_fold_change columns.",
            "markers": [],
        }

    by_symbol: dict[str, list[pd.Series]] = {}
    for _, row in table.iterrows():
        by_symbol.setdefault(str(row["gene"]).upper(), []).append(row)

    observations: list[dict[str, Any]] = []
    for marker in MARKER_GENES:
        matching_rows = by_symbol.get(marker.upper(), [])
        if not matching_rows:
            observations.append({"gene": marker, "status": "missing_gene"})
            continue
        if len(matching_rows) != 1:
            observations.append({"gene": marker, "status": "ambiguous_gene_symbol", "matches": len(matching_rows)})
            continue
        row = matching_rows[0]
        raw_effect = pd.to_numeric(pd.Series([row["log2_fold_change"]]), errors="coerce").iloc[0]
        if pd.isna(raw_effect) or not math.isfinite(float(raw_effect)):
            status = "not_estimable"
            effect = None
        elif float(raw_effect) > 0:
            status = "positive_direction"
            effect = float(raw_effect)
        else:
            status = "direction_not_positive"
            effect = float(raw_effect)
        observation = {"gene": marker, "status": status, "log2_fold_change": effect}
        for column, key in (("p_value", "p_value"), ("fdr_q_value", "fdr_q_value")):
            if column in row.index:
                value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
                observation[key] = None if pd.isna(value) or not math.isfinite(float(value)) else float(value)
        observations.append(observation)

    passed = all(item["status"] == "positive_direction" for item in observations)
    return {
        "status": "passed" if passed else "failed",
        "contrast": {"tested": "stim", "reference": "ctrl"},
        "expected_direction": "positive log2 fold change for stim versus ctrl",
        "criterion": "finite log2_fold_change > 0 for every pre-registered marker",
        "interpretation": "Direction check only; it does not assert significance or biological validity.",
        "reference": MARKER_REFERENCE_URL,
        "reference_context": (
            "The Seurat tutorial highlights these genes in aggregate control-versus-stimulated plots; "
            "the positive stim-versus-ctrl sign is the pre-registered expectation for this check."
        ),
        "markers": observations,
    }


def _runtime_packages(lock_path: Path) -> tuple[dict[str, str], str]:
    packages = _read_lock_packages(lock_path)
    # The package is intentionally excluded from the dependency lock, so record
    # the installed EACBP distribution version as an additional exact runtime pin.
    packages["eacbp"] = importlib.metadata.version("eacbp")
    return packages, _sha256_file(lock_path)


def _build_manifest(
    *,
    study_id: str,
    data_name: str,
    data_sha256: str,
    lock_path: Path,
    lock_sha256: str,
    packages: dict[str, str],
) -> dict[str, Any]:
    from eacbp.acceptance.manifest import AcceptanceManifest

    manifest = AcceptanceManifest.model_validate(
        {
            "schema_version": 1,
            "study_id": study_id,
            "dataset_kind": "real",
            "data_path": data_name,
            "sha256": data_sha256,
            "source": SOURCE_CITATION,
            "source_url": SOURCE_URL,
            "license": SOURCE_LICENSE,
            "species": "Homo sapiens",
            "counts_layer": "counts",
            "contrast": {"column": "label", "tested": "stim", "reference": "ctrl"},
            "donor": {
                "column": "replicate",
                "biological_unit": "donor",
                "paired": True,
                "minimum_per_condition": 2,
            },
            "batch_column": None,
            "covariates": [],
            "reference_environment": {
                "python_version": platform.python_version(),
                "packages": packages,
                "lockfile_path": str(lock_path),
                "lockfile_sha256": lock_sha256,
            },
        }
    )
    return manifest.model_dump(mode="json")


def _run_py_deseq2(adata: Any, manifest: Any, work_dir: Path) -> tuple[dict[str, Any], Any | None]:
    import pandas as pd

    from eacbp.acceptance.runner import _parameters
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.capabilities.advanced_statistics import PyDESeq2PseudobulkCapability
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskContract, TaskStatus

    result_summary: dict[str, Any]
    table: Any | None = None
    with tempfile.TemporaryDirectory(prefix="kang-artifacts-", dir=work_dir) as artifact_dir:
        registry = ArtifactRegistry(artifact_dir)
        input_uri = f"adata://{manifest.study_id}/public_counts/v1"
        output_uri = f"table://{manifest.study_id}/pydeseq2_markers/v1"
        registry.register(
            input_uri,
            adata,
            ArtifactType.ANNDATA,
            manifest.study_id,
            "verify_public_kang",
            "public_count_matrix_with_counts_layer",
            summary_metrics={"source_sha256": SOURCE_SHA256, "dataset_kind": "real"},
        )
        contract = TaskContract(
            task_id=f"{manifest.study_id}_pydeseq2",
            capability="deg",
            method="pydeseq2_pseudobulk_v1",
            input_artifacts=[input_uri],
            expected_outputs=[output_uri],
            validation_requirements=[
                "advanced_statistics_integrity",
                "multiple_testing_correction",
                "pseudoreplication_audit",
            ],
            parameters={**_parameters(manifest), "study_id": manifest.study_id},
        )
        result = PyDESeq2PseudobulkCapability().execute(contract, registry)
        result_summary = {
            "status": result.status.value,
            "method": result.method_used,
            "output_artifacts": list(result.output_artifacts),
            "error_type": result.error_type.value if result.error_type else None,
            "error_message": result.error_message,
            "inference_contract_id": result.inference_contract_id,
            "metrics": _json_safe(result.metrics),
            "scientific_result": _json_safe(result.scientific_result) if result.scientific_result else None,
        }
        if result.status != TaskStatus.SUCCESS:
            return result_summary, None
        table = registry.load_payload(output_uri)
        if not isinstance(table, pd.DataFrame):
            raise TypeError(f"PyDESeq2 capability returned {type(table).__name__}, expected pandas.DataFrame.")
        result_summary["n_rows"] = int(len(table))
        result_summary["columns"] = [str(column) for column in table.columns]
    return result_summary, table


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path, help="Downloaded public Kang H5AD; it is read only.")
    parser.add_argument(
        "--work-dir",
        required=True,
        type=Path,
        help="Absolute external directory for all generated outputs (not the repository or C: drive).",
    )
    parser.add_argument("--cell-type", default="CD14+ Monocytes", help="Exact label in obs['cell_type'] (default: CD14+ Monocytes).")
    parser.add_argument("--lockfile", type=Path, default=LOCKFILE_DEFAULT, help="Hashed scientific requirements lockfile.")
    return parser.parse_args(argv)


def _run_validation(
    args: argparse.Namespace,
    data_path: Path,
    work_dir: Path,
    lock_path: Path,
    source_sha: str,
) -> tuple[dict[str, Any], Path]:
    cell_type_slug = _safe_slug(args.cell_type)
    study_id = f"kang_2018_{cell_type_slug}"
    stem = f"{study_id}"
    work_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "data": work_dir / f"{stem}.h5ad",
        "manifest": work_dir / f"{stem}.manifest.json",
        "acceptance": work_dir / f"{stem}.acceptance.json",
        "table": work_dir / f"{stem}.pydeseq2.tsv",
        "markers": work_dir / f"{stem}.marker-direction.json",
        "summary": work_dir / f"{stem}.summary.json",
    }
    existing = [str(path) for path in paths.values() if path.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite existing outputs: " + ", ".join(existing))

    lock_packages, lock_sha = _runtime_packages(lock_path)
    import anndata as ad

    source = ad.read_h5ad(data_path)
    required_columns = {"cell_type", "label", "replicate"}
    missing = sorted(required_columns.difference(source.obs.columns))
    if missing:
        raise ValueError(f"Kang source is missing required obs columns: {missing}")
    if bool(source.uns.get("is_simulated", False)):
        raise ValueError("The source AnnData is marked is_simulated; refusing real-data acceptance.")
    _validate_raw_counts(source.X)

    labels = source.obs["cell_type"].astype("string")
    mask = labels.eq(args.cell_type).fillna(False).to_numpy(dtype=bool)
    selected_cells = int(mask.sum())
    if selected_cells == 0:
        observed = sorted(str(value) for value in source.obs["cell_type"].dropna().unique())
        raise ValueError(f"No cells have cell_type={args.cell_type!r}; available labels: {observed}")
    subset = source[mask, :].copy()
    _validate_raw_counts(subset.X)
    if subset.n_vars != source.n_vars or not subset.var_names.equals(source.var_names):
        raise RuntimeError("Cell-type selection changed the gene set or order.")
    if subset.n_obs != selected_cells or not subset.obs_names.equals(source.obs_names[mask]):
        raise RuntimeError("Cell-type selection did not retain every matching cell in source order.")
    if subset.obs[["label", "replicate"]].isna().any().any():
        raise ValueError("Selected cells contain missing condition or replicate metadata.")
    observed_conditions = set(subset.obs["label"].astype(str).unique())
    if observed_conditions != {"ctrl", "stim"}:
        raise ValueError(f"Expected only ctrl/stim labels for this contrast; observed {sorted(observed_conditions)}.")

    # The source X is the original integer count matrix. Donor, condition,
    # batch, genes, and cells are retained without relabeling or resampling.
    subset.layers["counts"] = subset.X.copy()
    subset.write_h5ad(paths["data"], compression="gzip")
    derived_sha = _sha256_file(paths["data"])
    manifest_dict = _build_manifest(
        study_id=study_id,
        data_name=paths["data"].name,
        data_sha256=derived_sha,
        lock_path=lock_path,
        lock_sha256=lock_sha,
        packages=lock_packages,
    )
    from eacbp.acceptance.manifest import load_manifest
    from eacbp.acceptance.runner import acceptance_exit_code, run_acceptance, write_result

    _write_json(manifest_dict, paths["manifest"])
    manifest, manifest_dir = load_manifest(paths["manifest"])
    acceptance = run_acceptance(manifest, manifest_dir)
    write_result(acceptance, paths["acceptance"])

    # Re-open the exact derived file used by acceptance, then execute the
    # public PyDESeq2 capability once more to persist the gene-level table.
    capability_summary: dict[str, Any]
    marker_report: dict[str, Any] | None = None
    result_table: Any | None = None
    try:
        analyzed = ad.read_h5ad(paths["data"])
        capability_summary, result_table = _run_py_deseq2(analyzed, manifest, work_dir)
        if result_table is not None:
            result_table.to_csv(paths["table"], sep="\t", index=False, na_rep="NA")
            marker_report = _marker_direction_report(result_table)
            _write_json(marker_report, paths["markers"])
    except Exception as exc:
        capability_summary = {"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)}

    summary = {
        "study_id": study_id,
        "dataset_kind": "real",
        "source": {
            "url": SOURCE_URL,
            "citation": SOURCE_CITATION,
            "license": SOURCE_LICENSE,
            "expected_sha256": SOURCE_SHA256,
            "actual_sha256": source_sha,
        },
        "subset": {
            "cell_type_column": "cell_type",
            "cell_type": args.cell_type,
            "n_cells": int(subset.n_obs),
            "n_genes": int(subset.n_vars),
            "condition_column": "label",
            "condition_counts": {str(k): int(v) for k, v in subset.obs["label"].value_counts().sort_index().items()},
            "donor_column": "replicate",
            "donor_counts_by_condition": {
                str(condition): int(subset.obs.loc[subset.obs["label"].astype(str) == str(condition), "replicate"].nunique())
                for condition in ("ctrl", "stim")
            },
            "counts_layer": "counts",
            "all_genes_and_matching_cells_retained": True,
            "batch_column": None,
        },
        "lockfile": {"path": str(lock_path), "sha256": lock_sha, "pinned_distributions": len(lock_packages)},
        "files": {key: {"path": str(path), "sha256": _sha256_file(path) if path.is_file() else None} for key, path in paths.items()},
        "acceptance": {
            "statuses": acceptance.get("statuses"),
            "research_acceptance": acceptance.get("research_acceptance"),
            "exit_code": acceptance_exit_code(acceptance),
            "report": str(paths["acceptance"]),
            "failure": acceptance.get("failure"),
        },
        "pydeseq2_capability": capability_summary,
        "marker_direction_check": marker_report,
        "overall_status": "passed" if (
            acceptance_exit_code(acceptance) == 0
            and capability_summary.get("status") == "success"
            and marker_report is not None
            and marker_report.get("status") == "passed"
        ) else "failed",
    }
    return summary, paths["summary"]


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    work_dir = _validate_external_work_dir(args.work_dir)
    data_path = args.data.expanduser().resolve()
    lock_path = args.lockfile.expanduser().resolve()
    if not data_path.is_file():
        raise FileNotFoundError(f"Kang source H5AD does not exist: {data_path}")

    source_sha = _sha256_file(data_path)
    study_id = f"kang_2018_{_safe_slug(args.cell_type)}"
    summary_path = work_dir / f"{study_id}.summary.json"
    work_dir.mkdir(parents=True, exist_ok=True)
    if source_sha != SOURCE_SHA256:
        integrity = _source_integrity_report(data_path, source_sha)
        if not summary_path.exists():
            _write_json(
                {
                    "study_id": study_id,
                    "dataset_kind": "real",
                    "source": {
                        "url": SOURCE_URL,
                        "citation": SOURCE_CITATION,
                        "license": SOURCE_LICENSE,
                        "expected_sha256": SOURCE_SHA256,
                        **integrity,
                    },
                    "overall_status": "failed",
                    "error": {
                        "type": "ValueError",
                        "message": f"Source H5AD SHA256 mismatch: expected {SOURCE_SHA256}, got {source_sha}.",
                    },
                },
                summary_path,
            )
        raise ValueError(f"Source H5AD SHA256 mismatch: expected {SOURCE_SHA256}, got {source_sha}.")

    try:
        summary, summary_path = _run_validation(args, data_path, work_dir, lock_path, source_sha)
    except Exception as exc:
        integrity = _source_integrity_report(data_path, source_sha)
        if not summary_path.exists() and not isinstance(exc, FileExistsError):
            failure_summary = {
                "study_id": study_id,
                "dataset_kind": "real",
                "source": {
                    "url": SOURCE_URL,
                    "citation": SOURCE_CITATION,
                    "license": SOURCE_LICENSE,
                    "expected_sha256": SOURCE_SHA256,
                    **integrity,
                },
                "overall_status": "failed",
                "error": {"type": type(exc).__name__, "message": str(exc)},
            }
            _write_json(failure_summary, summary_path)
        raise

    _apply_source_integrity(summary, data_path, source_sha)
    _write_json(summary, summary_path)
    print(json.dumps(_json_safe(summary), indent=2, ensure_ascii=False, allow_nan=False), flush=True)
    return 0 if summary["overall_status"] == "passed" else 1


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"verify_public_kang failed: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(2)
