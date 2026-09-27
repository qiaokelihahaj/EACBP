"""Executable acceptance checks that reuse the existing PyDESeq2 capability."""

from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import platform
import tempfile
import threading
import time
from pathlib import Path
from typing import Any, Iterator, Optional

import numpy as np
import pandas as pd
from scipy import sparse, stats

from eacbp.acceptance.manifest import AcceptanceManifest


_RESEARCH_PACKAGES = ("pydeseq2", "numpy", "pandas", "scipy", "formulaic", "anndata")


class AcceptanceFailure(ValueError):
    """A data or design prerequisite rejected before biological interpretation."""


class _RssSampler:
    def __init__(self, interval_seconds: float = 0.01):
        self.interval_seconds = interval_seconds
        self.peak_bytes: Optional[int] = None
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        try:
            import psutil
            self._process = psutil.Process()
        except ImportError:
            self._process = None

    def _sample(self):
        if self._process is None:
            return
        while not self._stop.is_set():
            try:
                rss = int(self._process.memory_info().rss)
                self.peak_bytes = max(self.peak_bytes or 0, rss)
            except Exception:
                return
            self._stop.wait(self.interval_seconds)

    def __enter__(self):
        if self._process is not None:
            self._thread = threading.Thread(target=self._sample, daemon=True)
            self._thread.start()
        return self

    def __exit__(self, exc_type, exc, tb):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=1.0)
        return False

    def as_dict(self) -> dict[str, Any]:
        return {
            "rss_peak_observed_bytes": self.peak_bytes,
            "rss_sampling_interval_seconds": self.interval_seconds if self._process else None,
            "rss_measurement": "sampled current-process RSS; short peaks may be missed" if self._process else "unavailable (psutil is not installed)",
        }


@contextmanager
def _measure() -> Iterator[dict[str, Any]]:
    started = time.perf_counter()
    with _RssSampler() as rss:
        metrics: dict[str, Any] = {}
        try:
            yield metrics
        finally:
            metrics["elapsed_seconds"] = max(0.0, time.perf_counter() - started)
            metrics.update(rss.as_dict())


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _installed_versions(extra_packages=()) -> dict[str, str]:
    result = {"python": platform.python_version()}
    packages = {"eacbp", "pydeseq2", "numpy", "pandas", "scipy", "formulaic", "anndata", "psutil"}
    packages.update(str(package) for package in extra_packages)
    for package in sorted(packages):
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = "unavailable"
    return result


def _check_reference_environment(manifest: AcceptanceManifest, manifest_dir: Path) -> dict[str, Any]:
    expected = manifest.reference_environment
    if manifest.dataset_kind == "synthetic":
        actual = _installed_versions()
        return {"status": "informational_only", "expected": None, "actual": actual,
                "note": "Synthetic regressions record the active runtime but do not establish a research lock."}
    assert expected is not None
    actual = _installed_versions(expected.packages.keys())
    required = {"python": expected.python_version, **expected.packages}
    missing = [name for name in _RESEARCH_PACKAGES if not required.get(name) or not _filled(str(required.get(name, "")))]
    mismatches = {
        name: {"expected": version, "actual": actual.get(name, "unavailable")}
        for name, version in required.items()
        if actual.get(name, "unavailable") != version
    }
    lock_path = Path(expected.lockfile_path or "")
    if not lock_path.is_absolute():
        lock_path = manifest_dir / lock_path
    lock_info: dict[str, Any] = {"path": str(lock_path), "expected_sha256": expected.lockfile_sha256}
    if lock_path.is_file():
        lock_info["actual_sha256"] = _sha256_file(lock_path)
        if lock_info["actual_sha256"] != expected.lockfile_sha256:
            mismatches["lockfile_sha256"] = {"expected": expected.lockfile_sha256,
                                              "actual": lock_info["actual_sha256"]}
    else:
        mismatches["lockfile_path"] = {"expected": str(lock_path), "actual": "missing"}
    return {"status": "passed" if not missing and not mismatches else "failed",
            "expected": required, "actual": actual, "missing_pins": missing,
            "mismatches": mismatches, "lockfile": lock_info}


def _parameters(manifest: AcceptanceManifest) -> dict[str, Any]:
    settings = manifest.statistics
    donor = manifest.donor
    contrast = manifest.contrast
    params: dict[str, Any] = {
        "condition_col": contrast.column,
        "condition_a": contrast.reference,
        "condition_b": contrast.tested,
        "contrast": [contrast.column, contrast.tested, contrast.reference],
        "donor_col": donor.column,
        "paired": donor.paired,
        "min_donors": donor.minimum_per_condition,
        "counts_layer": manifest.counts_layer,
        "covariates": list(manifest.covariates),
        "alpha": settings.alpha,
        "cooks_filter": settings.cooks_filter,
        "independent_filter": settings.independent_filter,
        "fit_type": settings.fit_type,
        "size_factors_fit_type": settings.size_factors_fit_type,
        "n_cpus": 1,
    }
    if manifest.design_formula:
        params["design_formula"] = manifest.design_formula
    return params


def _input_metrics(adata: Any, manifest: AcceptanceManifest, input_size: int) -> dict[str, Any]:
    matrix = adata.layers[manifest.counts_layer]
    if sparse.issparse(matrix):
        nnz = int(np.count_nonzero(matrix.data))
    else:
        nnz = int(np.count_nonzero(np.asarray(matrix)))
    obs = adata.obs
    condition = obs[manifest.contrast.column].astype("string")
    donor = obs[manifest.donor.column].astype("string")
    mask = condition.isin([manifest.contrast.tested, manifest.contrast.reference])
    selected = pd.DataFrame({"condition": condition[mask], "donor": donor[mask]})
    groups = selected.drop_duplicates().shape[0]
    donor_counts = {
        manifest.contrast.reference: int(selected.loc[selected.condition == manifest.contrast.reference, "donor"].nunique()),
        manifest.contrast.tested: int(selected.loc[selected.condition == manifest.contrast.tested, "donor"].nunique()),
    }
    batch_counts = {}
    if manifest.batch_column and manifest.batch_column in obs.columns:
        batch_counts = {str(key): int(value) for key, value in obs.loc[mask, manifest.batch_column].value_counts(dropna=False).items()}
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_nonzero_counts": nnz,
        "input_bytes": input_size,
        "counts_dtype": str(matrix.dtype),
        "n_donor_condition_groups_observed": int(groups),
        "donors_per_condition": donor_counts,
        "cells_per_batch": batch_counts,
    }


def _check_resource_limits(metrics: dict[str, Any], manifest: AcceptanceManifest) -> None:
    limits = manifest.resource_limits
    values = {
        "max_cells": metrics["n_cells"], "max_genes": metrics["n_genes"],
        "max_nonzero_counts": metrics["n_nonzero_counts"], "max_input_bytes": metrics["input_bytes"],
    }
    over = [f"{name}={value} exceeds limit {limit}" for name, value in values.items()
            if (limit := getattr(limits, name)) is not None and value > limit]
    if over:
        raise AcceptanceFailure("Resource preflight failed: " + "; ".join(over))


def _direct_reference(adata: Any, manifest: AcceptanceManifest, design: str) -> pd.DataFrame:
    """Independently aggregate cells then call the documented PyDESeq2 API."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    condition_column = manifest.contrast.column
    donor_column = manifest.donor.column
    observed_conditions = adata.obs[condition_column].astype(str)
    donor_values = adata.obs[donor_column].astype(str)
    selected_mask = observed_conditions.isin([manifest.contrast.tested, manifest.contrast.reference]).to_numpy()
    selected_positions = np.flatnonzero(selected_mask)
    selected_obs = pd.DataFrame({
        "donor": donor_values.iloc[selected_positions].to_numpy(),
        "condition": observed_conditions.iloc[selected_positions].to_numpy(),
        **{name: adata.obs[name].iloc[selected_positions].to_numpy() for name in manifest.covariates},
    })
    source_counts = adata.layers[manifest.counts_layer]
    gene_names = adata.var["gene_name"].astype(str).tolist() if "gene_name" in adata.var else adata.var.index.astype(str).tolist()
    count_rows: list[np.ndarray] = []
    metadata_rows: list[dict[str, Any]] = []
    for (donor, condition), positions in selected_obs.groupby(["donor", "condition"], sort=True, observed=True).indices.items():
        cell_positions = selected_positions[np.asarray(positions, dtype=int)]
        row = source_counts[cell_positions].sum(axis=0)
        count_rows.append(np.asarray(row).reshape(-1).astype(np.int64, copy=False))
        record: dict[str, Any] = {
            "donor": str(donor), "condition": str(condition),
            donor_column: str(donor), condition_column: str(condition),
        }
        for name in manifest.covariates:
            values = pd.unique(selected_obs.iloc[positions][name])
            if len(values) != 1 or pd.isna(values[0]):
                raise AcceptanceFailure(f"Covariate {name!r} is missing or changes within donor-condition group {(donor, condition)}")
            record[name] = values[0]
        metadata_rows.append(record)
    sample_ids = [f"{row['donor']}__{row['condition']}" for row in metadata_rows]
    count_frame = pd.DataFrame(count_rows, index=sample_ids, columns=gene_names, dtype=np.int64)
    metadata = pd.DataFrame(metadata_rows, index=sample_ids)
    n_cpus = 1
    settings = manifest.statistics
    dds = DeseqDataSet(
        counts=count_frame,
        metadata=metadata,
        design=design,
        min_replicates=manifest.donor.minimum_per_condition,
        n_cpus=n_cpus,
        quiet=True,
        fit_type=settings.fit_type,
        size_factors_fit_type=settings.size_factors_fit_type,
    )
    dds.deseq2()
    result = DeseqStats(
        dds,
        contrast=["condition", manifest.contrast.tested, manifest.contrast.reference],
        alpha=settings.alpha,
        cooks_filter=settings.cooks_filter,
        independent_filter=settings.independent_filter,
        quiet=True,
        n_cpus=n_cpus,
    )
    result.summary()
    return result.results_df.copy()


def _compare(eacbp_table: pd.DataFrame, reference: pd.DataFrame, manifest: AcceptanceManifest) -> dict[str, Any]:
    settings = manifest.statistics
    result = eacbp_table.set_index("gene").sort_index()
    reference = reference.copy()
    reference.index = reference.index.astype(str)
    reference = reference.sort_index()
    if not result.index.equals(reference.index):
        return {"status": "failed", "reason": "Gene sets or ordering differ between EACBP and direct PyDESeq2 reference."}
    z = float(stats.norm.ppf(1.0 - settings.alpha / 2.0))
    reference_low = pd.to_numeric(reference["log2FoldChange"], errors="coerce") - z * pd.to_numeric(reference["lfcSE"], errors="coerce")
    reference_high = pd.to_numeric(reference["log2FoldChange"], errors="coerce") + z * pd.to_numeric(reference["lfcSE"], errors="coerce")
    comparisons = {
        "effect": (result["log2_fold_change"], reference["log2FoldChange"]),
        "wald_ci_low": (result["ci_low"], reference_low),
        "wald_ci_high": (result["ci_high"], reference_high),
        "fdr": (result["fdr_q_value"], reference["padj"]),
    }
    metrics = {}
    failures = []
    for name, (left, right) in comparisons.items():
        left_values = pd.to_numeric(left, errors="coerce").to_numpy(dtype=float)
        right_values = pd.to_numeric(right, errors="coerce").to_numpy(dtype=float)
        left_nan, right_nan = np.isnan(left_values), np.isnan(right_values)
        finite_mask_equal = np.array_equal(np.isfinite(left_values), np.isfinite(right_values))
        positive_inf_equal = np.array_equal(np.isposinf(left_values), np.isposinf(right_values))
        negative_inf_equal = np.array_equal(np.isneginf(left_values), np.isneginf(right_values))
        equal_nan = np.array_equal(left_nan, right_nan)
        common = np.isfinite(left_values) & np.isfinite(right_values)
        if common.any():
            delta = np.abs(left_values[common] - right_values[common])
            max_abs = float(delta.max())
            close = np.isclose(left_values[common], right_values[common], rtol=settings.rtol, atol=settings.atol)
        else:
            max_abs, close = 0.0, np.asarray([], dtype=bool)
        nonfinite_masks_equal = finite_mask_equal and equal_nan and positive_inf_equal and negative_inf_equal
        passed = nonfinite_masks_equal and (not common.any() or bool(np.all(close)))
        metrics[name] = {"passed": passed, "n_finite_compared": int(common.sum()),
                         "n_nan_eacbp": int(left_nan.sum()), "n_nan_reference": int(right_nan.sum()),
                         "finite_masks_equal": finite_mask_equal,
                         "positive_infinity_masks_equal": positive_inf_equal,
                         "negative_infinity_masks_equal": negative_inf_equal,
                         "max_absolute_difference": max_abs,
                         "rtol": settings.rtol, "atol": settings.atol}
        if not passed:
            failures.append(name)
    return {"status": "passed" if not failures else "failed", "metrics": metrics,
            "failed_metrics": failures,
            "interval_definition": f"Wald normal interval: log2FoldChange +/- z_(1-alpha/2) * lfcSE, alpha={settings.alpha}"}


def _research_status(manifest: AcceptanceManifest) -> dict[str, Any]:
    if manifest.dataset_kind == "synthetic":
        return {"status": "not_applicable_synthetic", "biological_validity_claimed": False,
                "reason": "Synthetic fixtures verify software behavior only."}
    return {"status": "pending_domain_review", "biological_validity_claimed": False,
            "reason": "A matching implementation does not establish biological validity, study design appropriateness, or scientific interpretation."}


def run_acceptance(manifest: AcceptanceManifest, manifest_dir: Path) -> dict[str, Any]:
    """Run hash, input/design, EACBP, and direct-reference checks."""
    out: dict[str, Any] = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "study_id": manifest.study_id,
        "dataset_kind": manifest.dataset_kind,
        "evidence_class": "synthetic_software_regression" if manifest.dataset_kind == "synthetic" else "real_dataset_candidate",
        "statuses": {"engineering": "not_run", "method": "not_run"},
        "research_acceptance": _research_status(manifest),
        "checks": [],
        "runtime_versions": _installed_versions(),
    }
    data_path = Path(manifest.data_path).expanduser()
    if not data_path.is_absolute():
        data_path = manifest_dir / data_path
    data_path = data_path.resolve()
    out["input"] = {"path": str(data_path), "expected_sha256": manifest.sha256}
    try:
        if not data_path.is_file():
            raise AcceptanceFailure(f"Input file does not exist: {data_path}")
        input_size = data_path.stat().st_size
        out["input"]["input_bytes"] = input_size
        max_input_bytes = manifest.resource_limits.max_input_bytes
        if max_input_bytes is not None and input_size > max_input_bytes:
            raise AcceptanceFailure(
                f"Input file size {input_size:,} bytes exceeds max_input_bytes={max_input_bytes:,}; "
                "rejected before loading h5ad."
            )
        input_sha = _sha256_file(data_path)
        out["input"]["actual_sha256"] = input_sha
        if input_sha != manifest.sha256:
            raise AcceptanceFailure("Input SHA256 does not match the manifest.")
        out["checks"].append({"id": "input_sha256", "status": "passed"})
        try:
            import anndata as ad
        except ImportError as exc:
            raise AcceptanceFailure("Reading h5ad requires the optional anndata dependency.") from exc
        with _measure() as read_metrics:
            adata = ad.read_h5ad(data_path)
        out["input"]["read"] = read_metrics
        if manifest.counts_layer not in adata.layers:
            raise AcceptanceFailure(f"Raw counts layer {manifest.counts_layer!r} is missing.")
        if manifest.contrast.column not in adata.obs:
            raise AcceptanceFailure(f"Condition column {manifest.contrast.column!r} is missing.")
        if manifest.donor.column not in adata.obs:
            raise AcceptanceFailure(f"Donor column {manifest.donor.column!r} is missing.")
        if manifest.batch_column and manifest.batch_column not in adata.obs:
            raise AcceptanceFailure(f"Batch column {manifest.batch_column!r} is missing.")
        if manifest.species_metadata_key:
            species_value = adata.uns.get(manifest.species_metadata_key)
            if species_value is None:
                raise AcceptanceFailure(f"Species metadata key {manifest.species_metadata_key!r} is missing from AnnData.uns.")
            if str(species_value).strip().lower() != manifest.species.strip().lower():
                raise AcceptanceFailure(f"Species metadata says {species_value!r}, manifest says {manifest.species!r}.")
            out["checks"].append({"id": "species_metadata", "status": "passed", "value": str(species_value)})
        else:
            out["checks"].append({"id": "species_metadata", "status": "attested_in_manifest", "value": manifest.species})
        metrics = _input_metrics(adata, manifest, input_size)
        out["input"].update(metrics)
        _check_resource_limits(metrics, manifest)
        out["checks"].append({"id": "input_metadata_and_resource_limits", "status": "passed"})
        if bool(adata.uns.get("is_simulated", False)) and manifest.dataset_kind == "real":
            raise AcceptanceFailure("AnnData is marked is_simulated; it cannot enter real research acceptance.")
        out["checks"].append({"id": "dataset_kind_guard", "status": "passed" if manifest.dataset_kind == "real" else "synthetic_fixture"})

        env = _check_reference_environment(manifest, manifest_dir)
        out["reference_environment"] = env
        if env["status"] == "failed":
            out["checks"].append({"id": "fixed_reference_environment", "status": "failed", "details": env})
            out["statuses"]["engineering"] = "passed"
            out["statuses"]["method"] = "blocked_environment"
            if manifest.dataset_kind == "real":
                out["research_acceptance"] = {"status": "blocked", "biological_validity_claimed": False,
                                               "reason": "The declared fixed reference environment did not match the runtime."}
            return out
        out["checks"].append({"id": "fixed_reference_environment", "status": env["status"]})

        from eacbp.capabilities.advanced_statistics import (
            _prepare_pseudobulk,
            _validate_design_and_make_dds,
            PyDESeq2PseudobulkCapability,
        )
        from eacbp.artifact.registry import ArtifactRegistry
        from eacbp.schemas.artifact import ArtifactType
        from eacbp.schemas.task import TaskContract, TaskStatus

        params = _parameters(manifest)
        with _measure() as design_metrics:
            prepared = _prepare_pseudobulk(adata, params)
            _validate_design_and_make_dds(prepared, params)
        out["design"] = {
            **design_metrics,
            "formula": prepared.design,
            "rank": prepared.design_rank,
            "columns": list(prepared.design_columns),
            "n_pseudobulk_samples": int(len(prepared.metadata)),
            "n_donors_per_condition": {
                manifest.contrast.reference: len(prepared.donor_ids_a),
                manifest.contrast.tested: len(prepared.donor_ids_b),
            },
            "paired": manifest.donor.paired,
            "covariates": list(manifest.covariates),
        }
        working_limit = manifest.resource_limits.max_pseudobulk_working_bytes
        if working_limit is not None and prepared.estimated_aggregation_working_bytes > working_limit:
            raise AcceptanceFailure(
                "Estimated pseudobulk working memory exceeds configured resource limit: "
                f"{prepared.estimated_aggregation_working_bytes} > {working_limit} bytes."
            )
        out["design"]["estimated_aggregation_working_bytes"] = int(prepared.estimated_aggregation_working_bytes)
        out["checks"].append({"id": "donor_qualification_and_design_rank", "status": "passed"})
        out["statuses"]["engineering"] = "passed"

        with tempfile.TemporaryDirectory(prefix="eacbp-research-acceptance-") as artifact_dir:
            registry = ArtifactRegistry(artifact_dir)
            input_uri = f"adata://{manifest.study_id}/acceptance_input/v1"
            output_uri = f"table://{manifest.study_id}/acceptance_pydeseq2/v1"
            registry.register(input_uri, adata, ArtifactType.ANNDATA, manifest.study_id, "research_acceptance", "verified_h5ad_input")
            contract = TaskContract(
                task_id="research_acceptance_pydeseq2",
                capability="deg",
                method="pydeseq2_pseudobulk_v1",
                input_artifacts=[input_uri],
                expected_outputs=[output_uri],
                validation_requirements=["advanced_statistics_integrity", "multiple_testing_correction", "pseudoreplication_audit"],
                parameters={**params, "study_id": manifest.study_id},
            )
            with _measure() as fit_metrics:
                result = PyDESeq2PseudobulkCapability().execute(contract, registry)
            if result.status != TaskStatus.SUCCESS:
                raise AcceptanceFailure(f"EACBP capability returned {result.status.value}.")
            eacbp_table = registry.load_payload(output_uri)
            from eacbp.auditor import ScientificAuditor
            audit = ScientificAuditor().audit_task(contract, result, registry)
            scientific_passed = bool(audit.overall_passed)
            out["scientific_result"] = {
                "task_id": result.task_id,
                "status": result.status.value,
                "method": result.method_used,
                "output_artifacts": list(result.output_artifacts),
                "inference_contract_id": result.inference_contract_id,
                "scientific_result": result.scientific_result.model_dump(mode="json") if result.scientific_result else None,
                "metrics": result.metrics,
                "audit": {
                    "overall_passed": scientific_passed,
                    "checks": [check.model_dump(mode="json") for check in audit.checks],
                },
                **fit_metrics,
            }
        out["checks"].append({"id": "eacbp_pydeseq2_fit", "status": "passed"})
        out["checks"].append({"id": "independent_scientific_audit", "status": "passed" if scientific_passed else "failed",
                              "task_id": result.task_id,
                              "checks": [check.check_name for check in audit.checks if check.passed is False]})

        with _measure() as reference_metrics:
            reference = _direct_reference(adata, manifest, prepared.design)
        comparison = _compare(eacbp_table, reference, manifest)
        out["direct_reference"] = {
            "api": "pydeseq2.dds.DeseqDataSet + pydeseq2.ds.DeseqStats",
            "software_version": out["runtime_versions"].get("pydeseq2"),
            "design": prepared.design,
            "contrast": ["condition", manifest.contrast.tested, manifest.contrast.reference],
            "alpha": manifest.statistics.alpha,
            "cooks_filter": manifest.statistics.cooks_filter,
            "independent_filter": manifest.statistics.independent_filter,
            "metrics": reference_metrics,
            "comparison": comparison,
        }
        out["checks"].append({"id": "direct_pydeseq2_effect_interval_fdr_comparison", "status": comparison["status"],
                              "details": comparison})
        out["statuses"]["method"] = "passed" if comparison["status"] == "passed" and scientific_passed else "failed"
        out["research_acceptance"] = _research_status(manifest)
        out["checks"].append({"id": "research_interpretation", "status": out["research_acceptance"]["status"],
                              "details": out["research_acceptance"]})
    except Exception as exc:
        out["failure"] = {"type": type(exc).__name__, "message": str(exc)}
        out["checks"].append({"id": "acceptance_run", "status": "failed", "error": out["failure"]})
        if out["statuses"]["engineering"] == "not_run":
            out["statuses"]["engineering"] = "failed"
        if out["statuses"]["method"] == "not_run":
            out["statuses"]["method"] = "not_run"
        if manifest.dataset_kind == "real":
            out["research_acceptance"] = {"status": "blocked", "biological_validity_claimed": False,
                                           "reason": str(exc)}
    return out


def acceptance_exit_code(result: dict[str, Any]) -> int:
    return 0 if result.get("statuses", {}).get("engineering") == "passed" and result.get("statuses", {}).get("method") == "passed" else 1


def write_result(result: dict[str, Any], path: Optional[str | Path]) -> None:
    payload = json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if path is None:
        print(payload, end="")
        return
    target = Path(path).expanduser().resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(payload, encoding="utf-8")
