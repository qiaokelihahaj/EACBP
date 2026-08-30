"""
Normalization capability performing library size adjustment, log1p transformation, and HVG selection.
"""

import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


class NormalizationCapability(BaseCapability):
    """Normalizes single-cell raw counts to target CPM (e.g. 10,000) and log1p transforms."""

    def __init__(self):
        super().__init__(
            capability_name="normalization",
            implementation_id="sc_normalize_log1p_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        target_sum = contract.parameters.get("target_sum", 10000.0)
        n_top_genes = contract.parameters.get("n_top_genes", 300)

        # 1. Total count normalization & Log1p transformation
        if hasattr(data.X, "tocsr"):
            import scanpy as sc
            adata = data.to_anndata()
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars), subset=False)
            norm_data = SCData.from_anndata(adata)
            norm_data.uns["normalization"] = {"target_sum": target_sum, "log1p": True}
            hvg_mask = norm_data.var["highly_variable"].values if "highly_variable" in norm_data.var.columns else np.ones(data.n_vars, dtype=bool)
            mean_expr = float(norm_data.X.data.mean()) if hasattr(norm_data.X, "data") and len(norm_data.X.data) > 0 else 0.0
        else:
            counts_per_cell = data.X.sum(axis=1, keepdims=True)
            counts_per_cell[counts_per_cell == 0] = 1.0
            norm_X = (data.X / counts_per_cell) * target_sum

            # 2. Log1p transformation
            log_X = np.log1p(norm_X)

            # 3. Find Highly Variable Genes (HVG) based on variance
            gene_variances = np.var(log_X, axis=0)
            top_hvg_indices = np.argsort(gene_variances)[::-1][:min(n_top_genes, data.n_vars)]
            hvg_mask = np.zeros(data.n_vars, dtype=bool)
            hvg_mask[top_hvg_indices] = True

            norm_data = data.copy()
            norm_data.X = log_X
            norm_data.var["highly_variable"] = hvg_mask
            norm_data.var["variance"] = gene_variances
            norm_data.uns["normalization"] = {"target_sum": target_sum, "log1p": True}
            mean_expr = float(np.mean(log_X))

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/normalized/v2"

        registry.register(
            uri_str=out_uri,
            payload=norm_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="normalize_total_and_log1p",
            parent_uris=[in_uri],
            parameters={"target_sum": target_sum, "n_top_genes": n_top_genes},
            summary_metrics={
                "n_cells": norm_data.n_obs,
                "n_genes": norm_data.n_vars,
                "n_hvg": int(hvg_mask.sum()),
                "mean_expression": mean_expr,
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["normalize_counts", "log1p_transform", "select_highly_variable_genes"],
            metrics={
                "n_hvg": int(hvg_mask.sum()),
                "target_sum": target_sum,
            }
        )
