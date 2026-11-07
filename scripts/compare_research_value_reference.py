#!/usr/bin/env python3
"""Compare research-value full and donor-LOO tables with direct PyDESeq2 fits.

This audit recomputes the full model and every complete donor deletion from
the immutable cell-level subset. It measures independent reference-statistics
cost only; it is not a full GUI workflow benchmark or a measure of human time.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path, PurePath, PureWindowsPath
from statistics import NormalDist
from typing import Any, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DESIGN = "~condition + donor"
FULL_MODEL_KEY = "__full_model__"
PRIMARY_METRICS = {
    "effect": "log2_fold_change",
    "ci_low": "ci_low",
    "ci_high": "ci_high",
    "p_value": "p_value",
    "fdr_q_value": "fdr_q_value",
}


def _validate_external_path(work_dir: PurePath, project_root: PurePath) -> PurePath:
    """Validate an absolute output path without depending on the host OS."""
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute path")
    if isinstance(work_dir, PureWindowsPath) != isinstance(project_root, PureWindowsPath):
        raise ValueError("--work-dir and project root must use the same path style")
    try:
        work_dir.relative_to(project_root)
    except ValueError:
        pass
    else:
        raise ValueError("--work-dir must be outside the repository")
    if isinstance(work_dir, PureWindowsPath) and work_dir.drive.casefold() == "c:":
        raise ValueError("--work-dir must not use the C: drive")
    return work_dir


def _prepare_work_dir(work_dir: Path, project_root: Path = PROJECT_ROOT) -> Path:
    if not work_dir.is_absolute():
        raise ValueError("--work-dir must be an absolute path")
    resolved = work_dir.expanduser().resolve(strict=False)
    _validate_external_path(resolved, project_root.resolve())
    if resolved.exists():
        if not resolved.is_dir():
            raise NotADirectoryError(f"--work-dir is not a directory: {resolved}")
        if any(resolved.iterdir()):
            raise FileExistsError(f"--work-dir must be new or empty; refusing to overwrite: {resolved}")
    else:
        resolved.mkdir(parents=True, exist_ok=False)
    return resolved


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def _read_tsv(path: Path, required_columns: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        columns = set(reader.fieldnames or ())
        missing = sorted(required_columns - columns)
        if missing:
            raise ValueError(f"{path} is missing required TSV columns: {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"TSV contains no data rows: {path}")
    if any(None in row for row in rows):
        raise ValueError(f"TSV row has more fields than its header: {path}")
    if any(value is None for row in rows for value in row.values()):
        raise ValueError(f"TSV row has fewer fields than its header: {path}")
    return rows


def _as_float(value: Any, *, field: str) -> float:
    if value is None:
        return math.nan
    text = str(value).strip()
    if not text or text.casefold() in {"na", "nan", "null", "none", "<na>"}:
        return math.nan
    try:
        return float(text)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Non-numeric value for {field}: {value!r}") from exc


def _mask_kind(value: float) -> str:
    if math.isnan(value):
        return "nan"
    if math.isinf(value):
        return "positive_infinity" if value > 0 else "negative_infinity"
    return "finite"


def _index_rows(rows: Sequence[Mapping[str, Any]], label: str) -> dict[str, Mapping[str, Any]]:
    indexed: dict[str, Mapping[str, Any]] = {}
    for row in rows:
        gene = str(row.get("gene", "")).strip()
        if not gene:
            raise ValueError(f"{label} contains an empty gene identifier")
        if gene in indexed:
            raise ValueError(f"{label} contains duplicate gene identifier: {gene}")
        indexed[gene] = row
    return indexed


def _compare_gene_tables(
    primary_rows: Sequence[Mapping[str, Any]],
    reference_rows: Sequence[Mapping[str, Any]],
    *,
    rtol: float,
    atol: float,
) -> dict[str, Any]:
    """Compare every gene, retaining NaN and infinity masks in the check."""
    primary = _index_rows(primary_rows, "primary table")
    reference = _index_rows(reference_rows, "direct reference")
    primary_genes = set(primary)
    reference_genes = set(reference)
    missing_from_reference = sorted(primary_genes - reference_genes)
    extra_in_reference = sorted(reference_genes - primary_genes)
    if missing_from_reference or extra_in_reference:
        return {
            "status": "failed",
            "n_primary_genes": len(primary_genes),
            "n_reference_genes": len(reference_genes),
            "missing_from_reference": missing_from_reference[:25],
            "extra_in_reference": extra_in_reference[:25],
            "reason": "Gene sets differ; no rows were silently dropped.",
        }

    metrics: dict[str, Any] = {}
    for metric, column in PRIMARY_METRICS.items():
        n_finite = 0
        n_nan_primary = 0
        n_nan_reference = 0
        n_positive_inf_primary = 0
        n_positive_inf_reference = 0
        n_negative_inf_primary = 0
        n_negative_inf_reference = 0
        mask_mismatches: list[str] = []
        value_mismatches: list[str] = []
        max_absolute_difference = 0.0
        for gene in sorted(primary):
            left = _as_float(primary[gene].get(column), field=f"primary.{column}")
            right = _as_float(reference[gene].get(column), field=f"reference.{column}")
            left_kind, right_kind = _mask_kind(left), _mask_kind(right)
            n_nan_primary += left_kind == "nan"
            n_nan_reference += right_kind == "nan"
            n_positive_inf_primary += left_kind == "positive_infinity"
            n_positive_inf_reference += right_kind == "positive_infinity"
            n_negative_inf_primary += left_kind == "negative_infinity"
            n_negative_inf_reference += right_kind == "negative_infinity"
            if left_kind != right_kind:
                if len(mask_mismatches) < 25:
                    mask_mismatches.append(gene)
                continue
            if left_kind != "finite":
                continue
            n_finite += 1
            delta = abs(left - right)
            max_absolute_difference = max(max_absolute_difference, delta)
            if not math.isclose(left, right, rel_tol=rtol, abs_tol=atol):
                if len(value_mismatches) < 25:
                    value_mismatches.append(gene)
        passed = not mask_mismatches and not value_mismatches
        metrics[metric] = {
            "status": "passed" if passed else "failed",
            "n_genes": len(primary),
            "n_finite_compared": n_finite,
            "n_nan_primary": n_nan_primary,
            "n_nan_reference": n_nan_reference,
            "n_positive_infinity_primary": n_positive_inf_primary,
            "n_positive_infinity_reference": n_positive_inf_reference,
            "n_negative_infinity_primary": n_negative_inf_primary,
            "n_negative_infinity_reference": n_negative_inf_reference,
            "finite_nan_infinity_masks_equal": not mask_mismatches,
            "mask_mismatch_gene_examples": mask_mismatches,
            "value_mismatch_gene_examples": value_mismatches,
            "max_absolute_difference": max_absolute_difference if n_finite else None,
            "rtol": rtol,
            "atol": atol,
        }
    return {
        "status": "passed" if all(m["status"] == "passed" for m in metrics.values()) else "failed",
        "n_genes_compared": len(primary),
        "gene_sets_equal": True,
        "metrics": metrics,
    }


def _reference_rows(reference: Any, alpha: float) -> list[dict[str, Any]]:
    """Normalize direct PyDESeq2 output and derive the configured Wald CI."""
    if not reference.index.is_unique:
        raise ValueError("Direct PyDESeq2 result has duplicate gene identifiers")
    required = {"log2FoldChange", "lfcSE", "pvalue", "padj"}
    missing = sorted(required - set(reference.columns))
    if missing:
        raise ValueError(f"Direct PyDESeq2 result lacks columns: {missing}")
    z = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    rows = []
    for gene, row in reference.iterrows():
        effect = _as_float(row["log2FoldChange"], field="reference.log2FoldChange")
        se = _as_float(row["lfcSE"], field="reference.lfcSE")
        rows.append(
            {
                "gene": str(gene),
                "log2_fold_change": effect,
                "ci_low": effect - z * se,
                "ci_high": effect + z * se,
                "p_value": _as_float(row["pvalue"], field="reference.pvalue"),
                "fdr_q_value": _as_float(row["padj"], field="reference.padj"),
            }
        )
    return rows


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(_json_safe(dict(payload)), indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _format_tsv_value(value: Any) -> Any:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
    return value


def _write_reference_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    columns = [
        "fold_id", "left_out_donor", "gene", "log2_fold_change", "ci_low",
        "ci_high", "p_value", "fdr_q_value",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _format_tsv_value(row.get(key, "")) for key in columns})


def _check_primary_fold(rows: Sequence[Mapping[str, str]], fold_id: str) -> tuple[list[dict[str, str]], dict[str, Any]]:
    selected = [row for row in rows if row["left_out_donor"] == fold_id]
    expected_status = "full_model" if fold_id == FULL_MODEL_KEY else "estimated"
    statuses = sorted({row["fit_status"] for row in selected})
    is_full = fold_id == FULL_MODEL_KEY
    expected_fields = {
        "removed_sample_count": "0" if is_full else "2",
        "design": DESIGN,
        "n_pseudobulk_samples": "16" if is_full else "14",
        "n_donors_condition_a": "8" if is_full else "7",
        "n_donors_condition_b": "8" if is_full else "7",
        "paired": "true",
        "condition_a": "stim",
        "condition_b": "ctrl",
    }
    observed_fields = {
        field: sorted({str(row[field]).strip().casefold() for row in selected})
        for field in expected_fields
    }
    field_mismatches = {
        field: {"expected": expected.casefold(), "observed": observed_fields[field]}
        for field, expected in expected_fields.items()
        if observed_fields[field] != [expected.casefold()]
    }
    status_check = {
        "status": "passed" if selected and statuses == [expected_status] and not field_mismatches else "failed",
        "expected_fit_status": expected_status,
        "observed_fit_statuses": statuses,
        "n_rows": len(selected),
        "expected_fold_fields": expected_fields,
        "fold_field_mismatches": field_mismatches,
    }
    return selected, status_check


def _validate_primary_loo(rows: Sequence[Mapping[str, str]], donors: Sequence[str]) -> None:
    for row in rows:
        if row["left_out_donor"] != FULL_MODEL_KEY and row["left_out_donor"] not in donors:
            raise ValueError(f"LOO primary contains unknown donor fold: {row['left_out_donor']!r}")


def _fold_table_check(
    fold_rows: Sequence[Mapping[str, str]],
    reference_rows: Sequence[Mapping[str, Any]],
    *,
    fold_id: str,
    rtol: float,
    atol: float,
) -> dict[str, Any]:
    if not fold_rows:
        return {
            "status": "failed",
            "reason": f"Primary table has no rows for required fold {fold_id!r}; no reference genes were dropped.",
            "n_primary_genes": 0,
            "n_reference_genes": len(reference_rows),
        }
    return _compare_gene_tables(fold_rows, reference_rows, rtol=rtol, atol=atol)


def run_comparison(args: argparse.Namespace) -> dict[str, Any]:
    """Run the nine independent direct-reference models and persist evidence."""
    work_dir = _prepare_work_dir(args.work_dir)
    summary_path = work_dir / "run_summary.json"
    started = time.perf_counter()
    report: dict[str, Any] = {
        "schema_version": 1,
        "status": "failed",
        "study_id": None,
        "design": DESIGN,
        "folds": [],
        "direct_reference_total_seconds": 0.0,
        "human_elapsed_seconds": None,
        "efficiency_scope": (
            "Independent PyDESeq2 reference computation only; excludes GUI/workflow execution, "
            "data preparation outside this audit, and human time. No labor-savings claim is made."
        ),
        "interpretation_boundary": (
            "Donor-deletion fits are sensitivity analyses, not external validation or independent replication. "
            "Research interpretation remains pending domain review."
        ),
    }
    reference_output_rows: list[dict[str, Any]] = []
    direct_seconds = 0.0
    try:
        from eacbp.acceptance.manifest import load_manifest
        from eacbp.acceptance.runner import _check_reference_environment, _direct_reference

        manifest_path = args.manifest.expanduser().resolve()
        manifest, manifest_dir = load_manifest(manifest_path)
        if manifest.dataset_kind != "real":
            raise ValueError("Research-value reference comparison requires a real-data manifest")
        if manifest.contrast.tested != "stim" or manifest.contrast.reference != "ctrl":
            raise ValueError("Expected the preregistered Kang contrast: stim versus ctrl")
        if manifest.donor.paired is not True:
            raise ValueError("Expected a paired donor design")

        data_path = args.data.expanduser().resolve()
        manifest_data_path = Path(manifest.data_path).expanduser()
        if not manifest_data_path.is_absolute():
            manifest_data_path = manifest_dir / manifest_data_path
        manifest_data_path = manifest_data_path.resolve()
        if data_path != manifest_data_path:
            raise ValueError("--data must be the subset H5AD referenced by the manifest")
        if not data_path.is_file():
            raise FileNotFoundError(data_path)

        full_de_path = args.full_de.expanduser().resolve()
        loo_primary_path = args.loo_primary.expanduser().resolve()
        source_summary_path = args.source_summary.expanduser().resolve()
        source_summary = _read_json(source_summary_path)
        if source_summary.get("overall_status") != "passed":
            raise ValueError("Source Kang acceptance summary is not passed")
        if source_summary.get("study_id") != manifest.study_id:
            raise ValueError("Source summary and manifest study_id differ")
        capability = source_summary.get("pydeseq2_capability", {})
        capability_metrics = capability.get("metrics", {}) if isinstance(capability, dict) else {}
        if capability.get("status") != "success":
            raise ValueError("Source Kang PyDESeq2 full-model capability is not successful")
        if capability_metrics.get("design") != DESIGN:
            raise ValueError("Source Kang capability did not use the preregistered design")
        if int(capability_metrics.get("n_pseudobulk_samples", -1)) != 16:
            raise ValueError("Source Kang full model must contain 16 donor-condition pseudobulk samples")

        source_files = source_summary.get("files", {})
        expected_file_hashes = {
            "data": _sha256(data_path),
            "manifest": _sha256(manifest_path),
            "table": _sha256(full_de_path),
        }
        for key, actual_hash in expected_file_hashes.items():
            recorded = (source_files.get(key) or {}).get("sha256")
            if not recorded or recorded != actual_hash:
                raise ValueError(f"Source summary {key} SHA256 does not match the supplied file")
        if manifest.sha256 != expected_file_hashes["data"]:
            raise ValueError("Subset H5AD SHA256 does not match the manifest")

        environment = _check_reference_environment(manifest, manifest_dir)
        report["environment_check"] = environment
        if environment["status"] != "passed":
            raise ValueError("Active runtime does not match the manifest lockfile environment")

        primary_required = {
            "gene", *PRIMARY_METRICS.values(), "left_out_donor", "fit_status",
            "removed_sample_count", "design", "n_pseudobulk_samples",
            "n_donors_condition_a", "n_donors_condition_b", "paired",
            "condition_a", "condition_b",
        }
        loo_rows = _read_tsv(loo_primary_path, primary_required)
        full_rows = _read_tsv(full_de_path, {"gene", *PRIMARY_METRICS.values()})

        import anndata as ad

        adata = ad.read_h5ad(data_path)
        if manifest.counts_layer not in adata.layers:
            raise ValueError(f"Manifest counts layer {manifest.counts_layer!r} is missing")
        if bool(adata.uns.get("is_simulated", False)):
            raise ValueError("Subset H5AD is marked simulated")
        if manifest.contrast.column not in adata.obs or manifest.donor.column not in adata.obs:
            raise ValueError("Subset is missing contrast or donor metadata")
        if adata.obs[[manifest.contrast.column, manifest.donor.column]].isna().any().any():
            raise ValueError("Subset contains missing condition or donor identifiers")

        condition = adata.obs[manifest.contrast.column].astype(str)
        donor_values = adata.obs[manifest.donor.column].astype(str)
        donors = sorted(set(donor_values.tolist()))
        if len(donors) != 8:
            raise ValueError(f"Expected exactly 8 complete Kang donors, found {len(donors)}")
        observed_pairs = set(zip(donor_values.tolist(), condition.tolist()))
        expected_conditions = {manifest.contrast.reference, manifest.contrast.tested}
        incomplete = [
            donor for donor in donors
            if any((donor, level) not in observed_pairs for level in expected_conditions)
        ]
        if incomplete:
            raise ValueError(f"Every LOO unit must have both condition samples; incomplete donors: {incomplete}")
        observed_conditions = set(condition.tolist())
        if observed_conditions != expected_conditions:
            raise ValueError(f"Unexpected condition labels: {sorted(observed_conditions)}")
        _validate_primary_loo(loo_rows, donors)

        report.update(
            {
                "study_id": manifest.study_id,
                "cell_type": source_summary.get("subset", {}).get("cell_type"),
                "input": {
                    "data_path": str(data_path),
                    "sha256": expected_file_hashes["data"],
                    "n_cells": int(adata.n_obs),
                    "n_genes": int(adata.n_vars),
                    "donors": donors,
                    "conditions": sorted(observed_conditions),
                    "samples_full": len(observed_pairs),
                },
                "primary_inputs": {
                    "full_de_path": str(full_de_path),
                    "full_de_sha256": expected_file_hashes["table"],
                    "loo_primary_path": str(loo_primary_path),
                    "source_summary_path": str(source_summary_path),
                },
                "fit_count_expected": 1 + len(donors),
                "folds": [],
            }
        )
        rtol, atol = manifest.statistics.rtol, manifest.statistics.atol

        def fit_and_check(
            *,
            fold_id: str,
            left_out_donor: str | None,
            primary_source: str,
            expected_primary_rows: Sequence[Mapping[str, str]],
        ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
            nonlocal direct_seconds
            if left_out_donor is None:
                fold_adata = adata
            else:
                keep = donor_values != left_out_donor
                fold_adata = adata[keep.to_numpy(), :].copy()
            fit_started = time.perf_counter()
            try:
                direct = _direct_reference(fold_adata, manifest, DESIGN)
            except Exception as exc:
                elapsed = max(0.0, time.perf_counter() - fit_started)
                direct_seconds += elapsed
                return (
                    {
                        "fold_id": fold_id,
                        "left_out_donor": left_out_donor,
                        "n_donors": 8 if left_out_donor is None else 7,
                        "n_samples": 16 if left_out_donor is None else 14,
                        "status": "failed",
                        "elapsed_seconds": elapsed,
                        "error": f"{type(exc).__name__}: {exc}",
                    },
                    [],
                )
            elapsed = max(0.0, time.perf_counter() - fit_started)
            direct_seconds += elapsed
            normalized = _reference_rows(direct, manifest.statistics.alpha)
            comparison = _fold_table_check(
                expected_primary_rows,
                normalized,
                fold_id=fold_id,
                rtol=rtol,
                atol=atol,
            )
            fold_report = {
                "fold_id": fold_id,
                "left_out_donor": left_out_donor,
                "n_donors": 8 if left_out_donor is None else 7,
                "n_samples": 16 if left_out_donor is None else 14,
                "status": comparison["status"],
                "elapsed_seconds": elapsed,
                "primary_source": primary_source,
                "comparison": comparison,
            }
            for row in normalized:
                reference_output_rows.append(
                    {
                        "fold_id": fold_id,
                        "left_out_donor": left_out_donor or FULL_MODEL_KEY,
                        **row,
                    }
                )
            return fold_report, normalized

        full_loo_rows, full_status = _check_primary_fold(loo_rows, FULL_MODEL_KEY)
        full_report, full_reference = fit_and_check(
            fold_id=FULL_MODEL_KEY,
            left_out_donor=None,
            primary_source="full-de.tsv",
            expected_primary_rows=full_rows,
        )
        full_report["loo_primary_full_status"] = full_status
        full_report["loo_primary_full_comparison"] = _fold_table_check(
            full_loo_rows,
            full_reference,
            fold_id=FULL_MODEL_KEY,
            rtol=rtol,
            atol=atol,
        ) if full_reference else {"status": "not_run"}
        if full_status["status"] != "passed" or full_report.get("loo_primary_full_comparison", {}).get("status") != "passed":
            full_report["status"] = "failed"
        report["folds"].append(full_report)

        for donor in donors:
            fold_rows, fit_status = _check_primary_fold(loo_rows, donor)
            fold_report, _ = fit_and_check(
                fold_id=f"leave_out:{donor}",
                left_out_donor=donor,
                primary_source="loo_primary.tsv",
                expected_primary_rows=fold_rows,
            )
            fold_report["primary_fit_status"] = fit_status
            if fit_status["status"] != "passed":
                fold_report["status"] = "failed"
            report["folds"].append(fold_report)

        _write_reference_tsv(work_dir / "direct_reference_results.tsv", reference_output_rows)
        report["direct_reference_total_seconds"] = direct_seconds
        report["machine_elapsed_seconds"] = max(0.0, time.perf_counter() - started)
        report["human_elapsed_seconds"] = None
        report["status"] = "passed" if all(fold["status"] == "passed" for fold in report["folds"]) else "failed"
        report["exit_code"] = 0 if report["status"] == "passed" else 1
        report["output"] = {
            "run_summary": str(summary_path),
            "direct_reference_results": str(work_dir / "direct_reference_results.tsv"),
        }
    except Exception as exc:
        report["status"] = "failed"
        report["exit_code"] = 1
        report["error"] = {"type": type(exc).__name__, "message": str(exc)}
        report["direct_reference_total_seconds"] = direct_seconds
        report["machine_elapsed_seconds"] = max(0.0, time.perf_counter() - started)
    _write_json(summary_path, report)
    return report


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path, help="subset H5AD referenced by the acceptance manifest")
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--full-de", required=True, type=Path, help="Full-model primary DEG TSV")
    parser.add_argument("--loo-primary", required=True, type=Path, help="Research-value full/LOO primary TSV")
    parser.add_argument("--source-summary", required=True, type=Path, help="Passed public Kang summary with file hashes")
    parser.add_argument("--work-dir", required=True, type=Path, help="New or empty absolute external output directory")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    result = run_comparison(args)
    print(json.dumps(_json_safe(result), indent=2, ensure_ascii=False, allow_nan=False))
    return int(result.get("exit_code", 1))


if __name__ == "__main__":
    raise SystemExit(main())
