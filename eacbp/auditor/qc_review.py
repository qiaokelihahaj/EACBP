"""Independently reconstruct every raw-count QC review output.

The candidate branch records threshold exclusions; it does not diagnose damage,
restore cells, select markers, or certify reference identities. No executor
reconstruction helper is used by this auditor.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.schemas.artifact import ArtifactType


_REASONS = (
    "retained", "low_gene_count", "high_mito_percentage",
    "low_gene_count_and_high_mito_percentage",
)
_CELL_COLUMNS = (
    "cell_id", "n_genes", "total_counts", "mito_pct", "mito_pct_available",
    "low_gene_flag", "high_mito_flag", "retained", "reason",
)
_FEATURE_COLUMNS = (
    "reason", "gene_index", "gene_id", "gene", "n_cells", "n_detected",
    "total_counts", "detection_fraction", "mean_raw_count", "status",
)


def _same_json(actual, expected, label):
    # JSON spelling distinguishes a Boolean from an equal-valued integer.
    try:
        actual_json = json.dumps(actual, sort_keys=True, allow_nan=False,
                                 default=lambda value: value.item() if isinstance(value, np.generic) else value.tolist())
        expected_json = json.dumps(expected, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not finite JSON: {exc}") from exc
    if actual_json != expected_json:
        raise ValueError(f"{label} differs from independently reconstructed input")


def _same_counts(actual, expected, label, *, require_csr=False):
    if require_csr and not sparse.isspmatrix_csr(actual):
        raise ValueError(f"{label} must use CSR raw counts")
    if tuple(actual.shape) != tuple(expected.shape) or actual.dtype != expected.dtype:
        raise ValueError(f"{label} changed raw count shape or dtype")
    left, right = sparse.csr_matrix(actual, copy=True), sparse.csr_matrix(expected, copy=True)
    for matrix in (left, right):
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        matrix.sort_indices()
    if not (np.array_equal(left.indptr, right.indptr)
            and np.array_equal(left.indices, right.indices)
            and np.array_equal(left.data, right.data)):
        raise ValueError(f"{label} differs from the exact raw count subset")


def _numeric_column(table, name):
    return pd.to_numeric(table[name], errors="raise").to_numpy()


def _exact_column(table, name, expected):
    actual = _numeric_column(table, name)
    if not np.array_equal(actual, expected):
        raise ValueError(f"{name} differs from independently reconstructed counts")


def _bool_column(table, name, expected):
    actual = table[name]
    if actual.dtype != bool or not np.array_equal(actual.to_numpy(), expected):
        raise ValueError(f"{name} must contain the exact complete Boolean disposition")


def _float_column(table, name, expected):
    actual = _numeric_column(table, name).astype(float)
    if not np.allclose(actual, expected, rtol=1e-12, atol=1e-12, equal_nan=True):
        raise ValueError(f"{name} differs from independently reconstructed values")


class QCReviewValidator(BaseAuditor):
    def __init__(self):
        super().__init__("qc_review_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(
            auditor_name=self.auditor_name, target_task_id=contract.task_id,
            target_artifact_uri=result.output_artifacts[0] if result.output_artifacts else None,
        )
        if contract.capability != "qc_review":
            return report
        try:
            metrics = self._reconstruct(contract, result, registry)
            report.add_check(
                "qc_review_integrity", True, ValidationSeverity.ERROR,
                "Independently reconstructed all four outputs from original raw counts; "
                "candidate preservation is descriptive and restores no cells.", metrics,
            )
        except (ValueError, TypeError, KeyError, IndexError, AttributeError, OverflowError, AssertionError) as exc:
            report.add_check(
                "qc_review_integrity", False, ValidationSeverity.ERROR,
                f"Raw-count QC review cannot be independently verified: {exc}",
            )
        return report

    def _reconstruct(self, contract, result, registry):
        if result.method_used != "raw_count_qc_review_v1":
            raise ValueError("Unrecognized QC review method")
        if result.executed_operations != ["summarize_qc_candidates"]:
            raise ValueError("QC review cannot perform restoration, normalization or identity inference")
        if len(contract.input_artifacts) != 1 or result.input_artifacts != contract.input_artifacts:
            raise ValueError("Exactly one original raw-count input is required")
        if (contract.inference_contract is not None or result.inference_contract is not None
                or result.inference_contract_id is not None or result.scientific_result is not None):
            raise ValueError("QC review cannot carry an inferential or damage diagnosis result")
        input_uri = contract.input_artifacts[0]
        source_meta, source = registry.get_scdata(input_uri)
        if source_meta.type not in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA) or min(source.shape) == 0:
            raise ValueError("QC review requires a nonempty original single-cell input")
        cell_ids = source.obs.index.astype(str).to_numpy()
        if not source.obs.index.is_unique or len(set(cell_ids)) != source.n_obs:
            raise ValueError("Original cell IDs are duplicated or collide after string conversion")
        min_genes = contract.parameters.get("min_genes", 100)
        max_mito_pct = contract.parameters.get("max_mito_pct", 20.0)
        if isinstance(min_genes, bool) or not isinstance(min_genes, int) or min_genes < 0:
            raise ValueError("min_genes must be a nonnegative integer")
        if (isinstance(max_mito_pct, bool) or not isinstance(max_mito_pct, (float, int))
                or not np.isfinite(max_mito_pct) or not 0 <= max_mito_pct <= 100):
            raise ValueError("max_mito_pct must be a finite percentage")
        max_mito_pct = float(max_mito_pct)
        counts = source.layers["counts"] if "counts" in source.layers else source.X
        counts_source = "layers.counts" if "counts" in source.layers else "X"
        values = counts.data if sparse.issparse(counts) else np.asarray(counts)
        if (values.dtype.kind not in "biuf" or not np.isfinite(values).all()
                or np.any(values < 0) or not np.equal(values, np.floor(values)).all()):
            raise ValueError("Original counts must be finite nonnegative integers")
        if values.dtype.kind in "iu" and values.size:
            # Independently bound both reduction axes before a signed or
            # unsigned sum can wrap into another apparently finite value.
            largest = int(values.max())
            sum_limit = int(np.iinfo(np.uint64 if values.dtype == np.uint64 else np.int64).max)
            if largest * max(source.shape) > sum_limit:
                raise ValueError("Original raw count reductions could overflow")
        detected = np.asarray((counts > 0).sum(axis=1)).ravel()
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if not np.isfinite(totals).all() or (totals < 0).any():
            raise ValueError("Raw count reductions overflowed")
        n_mito = None
        if "percent_mito" in source.obs:
            mito_pct = pd.to_numeric(source.obs["percent_mito"], errors="raise").to_numpy(dtype=float)
            mito_available = np.ones(source.n_obs, dtype=bool)
            mito_source = "obs.percent_mito"
        elif "mt" in source.var:
            if source.var["mt"].dtype != bool or source.var["mt"].isna().any():
                raise ValueError("var.mt must contain complete Boolean annotations")
            mito_mask = source.var["mt"].to_numpy()
            mito_source = "var.mt"
        else:
            genes = source.var["gene_name"] if "gene_name" in source.var else source.var.index
            mito_mask = np.asarray(pd.Index(genes).astype(str).str.lower().str.startswith("mt-"))
            mito_source = "gene_symbol_prefix"
        if mito_source != "obs.percent_mito":
            n_mito = int(mito_mask.sum())
            if n_mito:
                mito_counts = np.asarray(counts[:, mito_mask].sum(axis=1)).ravel()
                mito_pct = np.divide(mito_counts * 100.0, totals,
                                     out=np.zeros(source.n_obs, dtype=float), where=totals > 0)
                mito_available = totals > 0
            else:
                mito_pct = np.zeros(source.n_obs, dtype=float)
                mito_available = np.zeros(source.n_obs, dtype=bool)
                mito_source = "unavailable"
        filtering_applied = mito_source != "unavailable"
        if (not np.isfinite(mito_pct).all() or (mito_pct < 0).any() or (mito_pct > 100).any()):
            raise ValueError("Mitochondrial percentages must be finite and between zero and one hundred")
        low = detected < min_genes
        high = mito_pct > max_mito_pct if filtering_applied else np.zeros(source.n_obs, dtype=bool)
        retained = ~(low | high)
        reasons = np.full(source.n_obs, "retained", dtype=object)
        reasons[low] = "low_gene_count"
        reasons[high] = "high_mito_percentage"
        reasons[low & high] = "low_gene_count_and_high_mito_percentage"
        candidate_mask = ~retained
        sid = source_meta.study_id
        outputs = {
            "summary": f"json://{sid}/qc_review/v1",
            "candidates": f"adata://{sid}/qc_candidates/v1",
            "cells": f"table://{sid}/qc_review_cells/v1",
            "features": f"table://{sid}/qc_review_features/v1",
        }
        uris = list(outputs.values())
        if contract.expected_outputs != uris or result.output_artifacts != uris:
            raise ValueError("All four declared QC review outputs are required in order")
        # A fixed interpretation statement is not an executor reconstruction helper.
        from eacbp.capabilities.qc_review import INTERPRETATION_WARNING
        expected_summary = {
            "schema_version": 1, "source_artifact_uri": input_uri,
            "source_sha256": source_meta.sha256_hash, "counts_source": counts_source,
            "min_genes": min_genes, "max_mito_pct": max_mito_pct,
            "initial_cells": source.n_obs, "candidate_cells": int(candidate_mask.sum()),
            "retained_cells": int(retained.sum()), "n_genes": source.n_vars,
            "mitochondrial_metric_source": mito_source, "mitochondrial_filtering_applied": filtering_applied,
            "n_mitochondrial_genes": n_mito,
            "cell_disposition_counts": {reason: int(np.count_nonzero(reasons == reason)) for reason in _REASONS},
            "zero_count_cells": int(np.count_nonzero(totals == 0)), "outputs": outputs,
            "restoration_performed": False, "damage_diagnosis": False, "reference_labels_used": False,
            "analysis_scope": "within_dataset_descriptive", "interpretation_warning": INTERPRETATION_WARNING,
        }
        result_summary = dict(result.metrics)
        if "execution_attempts" in result_summary:
            attempts = result_summary.pop("execution_attempts")
            if (isinstance(attempts, bool) or not isinstance(attempts, int)
                    or not 1 <= attempts <= contract.retry_policy.require_human_after):
                raise ValueError("TaskResult execution_attempts is outside the contracted retry limit")
        _same_json(result_summary, expected_summary, "TaskResult QC review metrics")
        kinds = (ArtifactType.JSON, ArtifactType.ANNDATA, ArtifactType.TABLE, ArtifactType.TABLE)
        operations = ("raw_count_qc_review_summary", "preserve_qc_candidates",
                      "qc_review_cell_dispositions", "qc_review_feature_counts")
        payloads = []
        for uri, kind, operation in zip(uris, kinds, operations):
            meta, payload = registry.get(uri)
            if (meta.type != kind or meta.parent_uris != [input_uri] or meta.created_by_task != contract.task_id
                    or meta.operation != operation or meta.study_id != sid):
                raise ValueError(f"QC review output {uri} has incorrect provenance")
            _same_json(meta.parameters, {"min_genes": min_genes, "max_mito_pct": max_mito_pct},
                       f"QC review output {uri} parameters")
            if set(meta.summary_metrics) - set(expected_summary) - {"is_simulated", "data_origin"}:
                raise ValueError(f"QC review output {uri} carries undeclared summary metrics")
            _same_json({key: meta.summary_metrics.get(key) for key in expected_summary}, expected_summary,
                       f"QC review output {uri} summary metrics")
            payloads.append(payload)
        summary, candidate_payload, cells, features = payloads
        _same_json(summary, expected_summary, "QC review summary")
        from eacbp.capabilities.sc_data import SCData
        candidate = SCData.from_payload(candidate_payload)
        if (candidate.raw is not None or candidate.obsm or candidate.obsp or candidate.varm or candidate.varp
                or set(candidate.layers) != {"counts"}):
            raise ValueError("Candidate branch must contain only raw-count slots")
        _same_json(candidate.uns, {"qc_review": expected_summary}, "Candidate branch provenance")
        pd.testing.assert_frame_equal(candidate.var, source.var, check_exact=True)
        expected_obs = source.obs.loc[candidate_mask].copy(deep=True)
        expected_obs["qc_review_reason"] = reasons[candidate_mask]
        expected_obs["qc_review_candidate"] = True
        expected_obs["qc_passed"] = False
        # H5AD can encode the newly emitted repeated string column as categorical.
        candidate_obs = candidate.obs.copy(deep=True)
        candidate_obs["qc_review_reason"] = candidate_obs["qc_review_reason"].astype(object)
        pd.testing.assert_frame_equal(candidate_obs, expected_obs, check_exact=True)
        expected_counts = counts[candidate_mask]
        _same_counts(candidate.X, expected_counts, "Candidate X", require_csr=True)
        _same_counts(candidate.layers["counts"], expected_counts, "Candidate counts layer")
        if not isinstance(cells, pd.DataFrame) or tuple(cells.columns) != _CELL_COLUMNS or len(cells) != source.n_obs:
            raise ValueError("Cell table must cover every original cell with the exact descriptive columns")
        if not np.array_equal(cells["cell_id"].astype(str).to_numpy(), cell_ids):
            raise ValueError("Cell table changed original cell identity or ordering")
        _exact_column(cells, "n_genes", detected)
        _exact_column(cells, "total_counts", totals)
        disposition_mito = mito_pct.copy()
        disposition_mito[~mito_available] = np.nan
        _float_column(cells, "mito_pct", disposition_mito)
        for name, expected in {"mito_pct_available": mito_available, "low_gene_flag": low,
                               "high_mito_flag": high, "retained": retained}.items():
            _bool_column(cells, name, expected)
        if not np.array_equal(cells["reason"].to_numpy(), reasons):
            raise ValueError("Cell dispositions differ from the frozen thresholds")
        if (not isinstance(features, pd.DataFrame) or tuple(features.columns) != _FEATURE_COLUMNS
                or len(features) != source.n_vars * len(_REASONS)
                or features.duplicated(["reason", "gene_index"]).any()
                or set(features["reason"]) != set(_REASONS)):
            raise ValueError("Feature table must cover every original feature for all four reasons")
        gene_ids = source.var.index.astype(str).to_numpy()
        genes = (source.var["gene_name"] if "gene_name" in source.var else source.var.index).astype(str).to_numpy()
        for reason in _REASONS:
            rows = features.loc[features["reason"] == reason].sort_values("gene_index")
            _exact_column(rows, "gene_index", np.arange(source.n_vars))
            if (not np.array_equal(rows["gene_id"].astype(str).to_numpy(), gene_ids)
                    or not np.array_equal(rows["gene"].astype(str).to_numpy(), genes)):
                raise ValueError("Feature table changed original feature IDs or symbols")
            group = counts[reasons == reason]
            size = group.shape[0]
            total = np.asarray(group.sum(axis=0)).ravel()
            n_detected = np.asarray((group > 0).sum(axis=0)).ravel()
            if not np.isfinite(total).all() or (total < 0).any():
                raise ValueError("Raw feature count reductions overflowed")
            _exact_column(rows, "n_cells", np.repeat(size, source.n_vars))
            _exact_column(rows, "n_detected", n_detected)
            _exact_column(rows, "total_counts", total)
            _float_column(rows, "detection_fraction", n_detected / size if size else np.full(source.n_vars, np.nan))
            _float_column(rows, "mean_raw_count", total / size if size else np.full(source.n_vars, np.nan))
            if not rows["status"].eq("described" if size else "no_cells").all():
                raise ValueError("Feature group status does not expose empty groups")
        return {"source_sha256": source_meta.sha256_hash, "initial_cells": source.n_obs,
                "candidate_cells": int(candidate_mask.sum()), "outputs_verified": 4,
                "feature_rows_verified": len(features), "reference_labels_used": False,
                "restoration_performed": False, "damage_diagnosis": False}
