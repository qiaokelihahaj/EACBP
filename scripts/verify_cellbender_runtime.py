"""Verify a real CellBender runtime through the EACBP adapter.

This is an explicit WSL smoke verifier, rather than a pytest test.  It runs
the existing CellBender capability (which invokes the supplied CLI through a
real subprocess), reorders a target cell/gene subset to exercise ID mapping,
audits the published artifacts with :class:`ScientificAuditor`, and writes a
JSON receipt.  A successful CLI and finite matrix never imply biological
quality certification: the receipt keeps that flag false, while actual
quality warnings produce a non-zero exit code.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd
from scipy import sparse


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_array(value: Any) -> np.ndarray:
    if sparse.issparse(value):
        return value.toarray()
    return np.asarray(value)


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (pd.Timestamp,)):
        return value.isoformat()
    return value


def _warning_strings(value: Any) -> list[str]:
    """Flatten sidecar warnings after registry serialization round-trips."""

    if isinstance(value, np.ndarray):
        # AnnData may deserialize an uns list as an ndarray.  Convert each
        # element recursively so the receipt keeps one warning per item
        # instead of serializing the complete array repr as one string.
        if value.ndim == 0:
            return _warning_strings(value.item())
        flattened: list[str] = []
        for item in value:
            flattened.extend(_warning_strings(item))
        return flattened
    if isinstance(value, np.generic):
        return _warning_strings(value.item())
    if isinstance(value, (list, tuple, set)):
        flattened: list[str] = []
        for item in value:
            flattened.extend(_warning_strings(item))
        return flattened
    if isinstance(value, str):
        text = value.strip()
        if text.startswith("[") and text.endswith("]"):
            try:
                parsed = ast.literal_eval(text)
            except (SyntaxError, ValueError):
                parsed = None
            if isinstance(parsed, (list, tuple, set)):
                return _warning_strings(parsed)
        return [value]
    return [str(value)]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Unfiltered CellBender input .h5ad")
    parser.add_argument("--executable", required=True, help="CellBender executable path")
    parser.add_argument("--output", required=True, help="CellBender output .h5 path")
    parser.add_argument("--workdir", required=True, help="Fresh writable verification directory")
    parser.add_argument("--checkpoint", required=True, help="Existing CellBender checkpoint tarball")
    parser.add_argument("--target-cells", type=int, default=12)
    parser.add_argument("--target-genes", type=int, default=8)
    parser.add_argument("--expected-cells", type=int, default=500)
    parser.add_argument("--total-droplets-included", type=int, default=2000)
    parser.add_argument("--cpu-threads", type=int, default=4)
    parser.add_argument("--force-use-checkpoint", action="store_true")
    return parser


def _resolve_executable(raw: str) -> Path:
    candidate = Path(raw).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    located = shutil.which(raw)
    if located:
        return Path(located).resolve()
    raise FileNotFoundError(f"CellBender executable is unavailable: {raw}")


def _build_target(input_path: Path, workdir: Path, target_cells: int, target_genes: int):
    import anndata as ad

    source = ad.read_h5ad(input_path)
    if source.n_obs < target_cells or source.n_vars < target_genes:
        raise ValueError(
            f"Input shape {source.shape} is smaller than target subset "
            f"({target_cells}, {target_genes})"
        )
    if target_cells < 2 or target_genes < 2:
        raise ValueError("target-cells and target-genes must each be at least 2")
    # Reverse both axes so the adapter cannot pass an order-only check.  The
    # external output still contains the full source IDs and is mapped back by
    # the CellBender capability.
    cell_positions = np.arange(target_cells - 1, -1, -1, dtype=int)
    gene_positions = np.arange(target_genes - 1, -1, -1, dtype=int)
    subset = source[cell_positions, gene_positions].copy()
    cells = np.asarray(subset.obs_names.astype(str))
    genes = np.asarray(subset.var_names.astype(str))
    if len(set(cells)) != len(cells) or len(set(genes)) != len(genes):
        raise ValueError("Selected target cell and gene IDs must be unique")
    counts = subset.X.copy() if sparse.issparse(subset.X) else np.array(subset.X, copy=True)
    values = _matrix_array(counts)
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values)):
        raise ValueError("The official tiny input does not contain integer-valued nonnegative raw counts")
    obs = subset.obs.copy()
    var = subset.var.copy()
    obs["cell_id"] = cells
    var["gene_name"] = genes
    target = {
        "X": counts,
        "obs": obs,
        "var": var,
        "layers": {"counts": counts.copy() if sparse.issparse(counts) else np.array(counts, copy=True)},
        "obsm": {},
        "uns": {"cellbender_runtime_verification": {"source_input": str(input_path)}},
    }
    # Keep a readable target snapshot beside the receipt; the adapter input is
    # the registry artifact built from the same object.
    target_path = workdir / "target_subset.h5ad"
    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy() if sparse.issparse(counts) else np.array(counts, copy=True)
    adata.write_h5ad(target_path)
    return target, {
        "target_snapshot": str(target_path),
        "target_cells": cells.tolist(),
        "target_genes": genes.tolist(),
        "target_shape": [int(target_cells), int(target_genes)],
        "source_shape": [int(source.n_obs), int(source.n_vars)],
        "cell_order_reversed": True,
        "gene_order_reversed": True,
    }


def _write_receipt(workdir: Path, receipt: Mapping[str, Any]) -> Path:
    path = workdir / "verify_cellbender_runtime.json"
    path.write_text(
        json.dumps(_json_value(receipt), indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )
    return path


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    workdir = Path(args.workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    receipt: dict[str, Any] = {
        "status": "failed",
        "scientific_quality_certified": False,
        "arguments": vars(args),
    }
    try:
        input_path = Path(args.input).expanduser().resolve()
        output_path = Path(args.output).expanduser().resolve()
        checkpoint = Path(args.checkpoint).expanduser().resolve()
        executable = _resolve_executable(args.executable)
        for label, path in (("input", input_path), ("checkpoint", checkpoint)):
            if not path.is_file():
                raise FileNotFoundError(f"{label} does not exist: {path}")
        if output_path.suffix.casefold() != ".h5":
            raise ValueError("--output must be a CellBender .h5 path")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.exists():
            raise FileExistsError(f"Refusing to overwrite existing CellBender output: {output_path}")
        if args.cpu_threads < 1 or args.expected_cells < 1 or args.total_droplets_included < 1:
            raise ValueError("CellBender resource arguments must be positive")
        # The adapter launches the explicit executable but CellBender invokes
        # Jupyter for its HTML report by name.  Put the supplied runtime first
        # in PATH so report generation is tested in the same environment.
        os.environ["PATH"] = str(executable.parent) + os.pathsep + os.environ.get("PATH", "")

        # Import the project from this checkout even when the runtime's
        # editable install points at a different working directory.
        root = _repo_root()
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))
        from eacbp.artifact.registry import ArtifactRegistry
        from eacbp.capabilities.advanced_qc import CellBenderBackgroundRemovalCapability
        from eacbp.capabilities.sc_data import SCData
        from eacbp.auditor import ScientificAuditor
        from eacbp.orchestrator.resources import pin_resource_files
        from eacbp.schemas.artifact import ArtifactType
        from eacbp.schemas.task import TaskContract

        target_payload, target_info = _build_target(
            input_path,
            workdir,
            args.target_cells,
            args.target_genes,
        )
        target_data = SCData.from_dict(target_payload)
        registry = ArtifactRegistry(str(workdir / "artifacts"))
        study_id = "cellbender_runtime_verification"
        input_uri = f"adata://{study_id}/target_input/v1"
        registry.register(
            input_uri,
            target_data.to_dict(),
            ArtifactType.ANNDATA,
            study_id,
            "verify_cellbender_runtime_input",
            "prepare_reordered_target_subset",
            summary_metrics=target_info,
        )
        checkpoint_hash = _sha256(checkpoint)
        input_hash = _sha256(input_path)
        executable_hash = _sha256(executable)
        extra_args = [
            "--checkpoint", str(checkpoint),
            "--expected-cells", str(args.expected_cells),
            "--total-droplets-included", str(args.total_droplets_included),
            "--cpu-threads", str(args.cpu_threads),
        ]
        if args.force_use_checkpoint:
            extra_args.append("--force-use-checkpoint")
        parameters = pin_resource_files({
            "unfiltered_input_path": str(input_path),
            "executable": str(executable),
            "output_path": str(output_path),
            "run_cwd": str(workdir),
            "extra_args": extra_args,
            "cell_id_key": "cell_id",
            "gene_id_key": "gene_name",
        })
        task = TaskContract(
            task_id="verify_cellbender_runtime",
            capability="background_removal",
            method="cellbender_cli_v1",
            input_artifacts=[input_uri],
            expected_outputs=[
                f"adata://{study_id}/corrected/v1",
                f"json://{study_id}/run_report/v1",
            ],
            validation_requirements=["advanced_qc_integrity"],
            parameters=parameters,
        )
        result = CellBenderBackgroundRemovalCapability().execute(task, registry)
        if len(result.output_artifacts) != 2:
            raise RuntimeError(f"Expected AnnData and JSON outputs, received {result.output_artifacts}")
        output_meta, output_payload = registry.get(result.output_artifacts[0])
        output = output_payload if isinstance(output_payload, SCData) else SCData.from_dict(output_payload)
        _, report_payload = registry.get(result.output_artifacts[1])
        info = output.uns.get("cellbender", {})
        alignment = info.get("cell_gene_alignment", {}) if isinstance(info, Mapping) else {}
        raw_equal = np.array_equal(_matrix_array(target_data.layers["counts"]), _matrix_array(output.layers["counts"]))
        if not raw_equal:
            raise RuntimeError("Adapter output did not preserve the target raw counts layer")
        if alignment.get("cell_order_exact") is not False or alignment.get("gene_order_exact") is not False:
            raise RuntimeError("ID mapping smoke did not exercise reordered cell and gene identifiers")
        if alignment.get("extra_external_cells", 0) <= 0 or alignment.get("extra_external_genes", 0) <= 0:
            raise RuntimeError("CellBender output did not contain the expected full-matrix extras")

        audit = ScientificAuditor().audit_task(task, result, registry)
        quality_warnings = _warning_strings(info.get("quality_warnings", [])) if isinstance(info, Mapping) else []
        convergence_warnings = [
            str(item) for item in quality_warnings
            if any(token in str(item).casefold() for token in ("convergence", "elbo", "suboptimal", "could be"))
        ]
        report_available = bool(info.get("report_available")) if isinstance(info, Mapping) else False
        certified = bool(info.get("scientific_quality_certified")) if isinstance(info, Mapping) else False
        runtime_audit_passed = bool(audit.overall_passed)
        # The adapter intentionally leaves scientific_quality_certified=False.
        # A clean runtime can therefore return zero while the receipt still
        # says that biology was not certified.  Actual warnings remain a
        # non-zero quality caveat.
        runtime_quality_evidence_clean = bool(runtime_audit_passed and report_available and not quality_warnings)
        receipt.update({
            "status": "runtime_verified_quality_not_certified" if runtime_audit_passed else "audit_failed",
            "task_result": result.model_dump(mode="json"),
            "output_artifacts": result.output_artifacts,
            "output_metadata": output_meta.model_dump(mode="json"),
            "adapter_report": report_payload,
            "target": target_info,
            "alignment": alignment,
            "raw_counts_preserved": raw_equal,
            "audit_overall_passed": runtime_audit_passed,
            "audit_checks": [check.model_dump(mode="json") for check in audit.checks],
            "quality_status": info.get("quality_status") if isinstance(info, Mapping) else None,
            "quality_warnings": quality_warnings,
            "convergence_warnings": convergence_warnings,
            "report_available": report_available,
            "scientific_quality_certified": certified,
            "runtime_quality_evidence_clean": runtime_quality_evidence_clean,
            "scientific_quality_passed": False,
            "input_sha256": input_hash,
            "checkpoint_sha256": checkpoint_hash,
            "executable_sha256": executable_hash,
        })
        receipt_path = _write_receipt(workdir, receipt)
        print(json.dumps({
            "status": receipt["status"],
            "audit_overall_passed": runtime_audit_passed,
            "quality_status": receipt["quality_status"],
            "convergence_warnings": convergence_warnings,
            "receipt": str(receipt_path),
        }, indent=2, ensure_ascii=False))
        # A successful subprocess and structural audit are useful runtime
        # evidence, but this script never presents a zero exit as a
        # biological-quality certification.  The receipt always keeps that
        # flag false; warnings are surfaced as code 2 for CI/manual review.
        return 0 if runtime_quality_evidence_clean else (1 if not runtime_audit_passed else 2)
    except Exception as exc:
        receipt.update({
            "status": "failed",
            "error_type": type(exc).__name__,
            "error": str(exc),
        })
        receipt_path = _write_receipt(workdir, receipt)
        print(f"CellBender runtime verification failed: {type(exc).__name__}: {exc}", file=sys.stderr)
        print(f"Receipt: {receipt_path}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
