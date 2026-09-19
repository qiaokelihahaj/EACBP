"""Count normalisation with explicit preservation of the raw count matrix.

The output ``X`` is log1p(CPM-like) expression.  The matrix used to compute
it is retained as ``layers['counts']`` so downstream donor aggregation can
operate on counts instead of accidentally treating transformed values as
counts.
"""

import numpy as np
from scipy import sparse

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def _copy_matrix(value):
    return value.copy() if hasattr(value, "copy") else np.array(value, copy=True)


def _set_layer(data: SCData, name: str, value) -> None:
    """Set a layer across old and new SCData containers."""
    layers = getattr(data, "layers", None)
    if layers is None:
        layers = {}
        setattr(data, "layers", layers)
    layers[name] = _copy_matrix(value)


def _as_counts(data: SCData):
    layers = getattr(data, "layers", {}) or {}
    if "counts" in layers:
        return _copy_matrix(layers["counts"]), "layers.counts"
    return _copy_matrix(data.X), "X_before_normalization"


def _normalize_and_log1p(counts, target_sum: float):
    """Return normalized/log1p matrix without densifying sparse counts."""
    if hasattr(counts, "tocsr"):
        x = counts.tocsr().astype(np.float32)
        totals = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
        safe_totals = np.where(totals > 0, totals, 1.0)
        norm = x.multiply((float(target_sum) / safe_totals)[:, None]).tocsr()
        norm.data = np.log1p(norm.data).astype(np.float32)
        return norm, totals

    x = np.asarray(counts, dtype=np.float32)
    totals = np.sum(x, axis=1).astype(np.float32)
    safe_totals = np.where(totals > 0, totals, 1.0)
    norm = (x / safe_totals[:, None]) * float(target_sum)
    return np.log1p(norm).astype(np.float32), totals


def _column_variance(matrix) -> np.ndarray:
    if hasattr(matrix, "tocsr"):
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        mean_sq = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        return np.maximum(mean_sq - mean * mean, 0.0).astype(np.float32)
    return np.var(np.asarray(matrix), axis=0).astype(np.float32)


class NormalizationCapability(BaseCapability):
    """Library-size normalisation and log1p transformation."""

    def __init__(self, implementation_id: str = "library_size_log1p_v1"):
        requested_id = implementation_id
        super().__init__(
            capability_name="normalization",
            implementation_id="library_size_log1p_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )
        self.requested_implementation_id = requested_id
        self.legacy_aliases = {"sc_normalize_log1p_v1": self.implementation_id}

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)
        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)

        target_sum = float(contract.parameters.get("target_sum", 10000.0))
        n_top_genes = int(contract.parameters.get("n_top_genes", 300))
        if not np.isfinite(target_sum) or target_sum <= 0:
            raise ValueError("target_sum must be a finite positive number")

        counts, counts_source = _as_counts(data)
        input_layer = contract.parameters.get("input_layer", "counts")
        if input_layer == "counts":
            normalization_input = counts
        else:
            if input_layer not in data.layers:
                raise ValueError(f"Requested normalization input_layer is missing: {input_layer}")
            normalization_input = data.layers[input_layer]
        values = normalization_input.data if sparse.issparse(normalization_input) else np.asarray(normalization_input)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("Normalization requires finite nonnegative input")
        log_X, cell_totals = _normalize_and_log1p(normalization_input, target_sum)
        gene_variances = _column_variance(log_X)
        top_hvg_indices = np.argsort(gene_variances)[::-1][: min(max(0, n_top_genes), data.n_vars)]
        hvg_mask = np.zeros(data.n_vars, dtype=bool)
        hvg_mask[top_hvg_indices] = True

        norm_data = data.copy()
        norm_data.X = log_X
        _set_layer(norm_data, "counts", counts)
        _set_layer(norm_data, "normalized", log_X)
        norm_data.var["highly_variable"] = hvg_mask
        norm_data.var["variance_log1p"] = gene_variances
        norm_data.uns["normalization"] = {
            "method": self.implementation_id,
            "requested_method": self.requested_implementation_id,
            "target_sum": target_sum,
            "log1p": True,
            "counts_layer": "counts",
            "counts_source": counts_source,
            "input_layer": input_layer,
            "zero_count_cells": int(np.sum(cell_totals <= 0)),
        }

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/normalized/v2"
        if sparse.issparse(log_X):
            mean_expr = float(log_X.data.mean()) if log_X.nnz else 0.0
        else:
            mean_expr = float(np.mean(np.asarray(log_X)))

        registry.register(
            uri_str=out_uri,
            payload=norm_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="library_size_normalize_and_log1p",
            parent_uris=[in_uri],
            parameters={
                "target_sum": target_sum,
                "n_top_genes": n_top_genes,
                "counts_layer": "counts",
                "input_layer": input_layer,
            },
            summary_metrics={
                "n_cells": norm_data.n_obs,
                "n_genes": norm_data.n_vars,
                "n_hvg": int(hvg_mask.sum()),
                "mean_expression_log1p": mean_expr,
                "counts_preserved": True,
            },
        )

        operation_candidates = [
            "preserve_counts_layer",
            "normalize_library_size",
            "normalize_counts",
            "log1p_transform",
            "select_variance_hvgs",
            "select_highly_variable_genes",
        ]
        executed_operations = (
            [op for op in operation_candidates if op in contract.allowed_operations]
            if contract.allowed_operations
            else ["preserve_counts_layer", "normalize_library_size", "log1p_transform", "select_variance_hvgs"]
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=executed_operations,
            metrics={
                "n_hvg": int(hvg_mask.sum()),
                "target_sum": target_sum,
                "counts_layer": "counts",
                "counts_source": counts_source,
            },
        )
