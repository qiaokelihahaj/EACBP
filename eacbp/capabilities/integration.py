"""
Integration capability evaluating batch effects and applying batch integration methods.
"""

from typing import Dict, Any, List, Optional
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def compute_pca(X: np.ndarray, n_components: int = 30) -> np.ndarray:
    """Standard SVD / PCA dimensionality reduction."""
    X_centered = X - np.mean(X, axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    k = min(n_components, X.shape[0], X.shape[1])
    return U[:, :k] * S[:k]


def calculate_batch_mixing_score(embedding: np.ndarray, batches: np.ndarray, k: int = 15) -> float:
    """Calculates average fraction of nearest neighbors from different batches."""
    if len(np.unique(batches)) <= 1 or embedding.shape[0] < k:
        return 1.0
    
    # Subsample if too large for speed
    n = embedding.shape[0]
    indices = np.random.choice(n, size=min(n, 300), replace=False) if n > 300 else np.arange(n)
    dists = cdist(embedding[indices], embedding)
    
    mixing_scores = []
    for i, idx in enumerate(indices):
        nn_indices = np.argsort(dists[i])[1:k+1]
        same_batch_count = (batches[nn_indices] == batches[idx]).sum()
        # Mixing score is higher when neighbors are diverse
        mixing_scores.append(1.0 - (same_batch_count / k))
    
    return float(np.mean(mixing_scores))


class IntegrationCapability(BaseCapability):
    """Evaluates batch effects and performs integration (e.g. Harmony-like or No-Correction)."""

    def __init__(self, implementation_id: str = "harmony"):
        super().__init__(
            capability_name="integration",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        # Use highly variable genes if available
        if "highly_variable" in data.var.columns:
            hvg_mask = data.var["highly_variable"].values
            X_use = data.X[:, hvg_mask]
        else:
            X_use = data.X

        n_comps = contract.parameters.get("n_components", 20)
        pca_emb = compute_pca(X_use, n_components=n_comps)

        batch_col = "batch" if "batch" in data.obs.columns else None
        batches = data.obs[batch_col].values if batch_col else np.array(["b1"] * data.n_obs)

        pre_mixing = calculate_batch_mixing_score(pca_emb, batches)

        if self.implementation_id == "harmony" and batch_col:
            # Batch alignment adjustment on PCA space
            adjusted_pca = pca_emb.copy()
            unique_batches = np.unique(batches)
            global_mean = np.mean(adjusted_pca, axis=0, keepdims=True)
            for b in unique_batches:
                mask = (batches == b)
                batch_mean = np.mean(adjusted_pca[mask], axis=0, keepdims=True)
                adjusted_pca[mask] -= 0.8 * (batch_mean - global_mean)
            integrated_emb = adjusted_pca
            method_desc = "harmony_pca_alignment"
        else:
            integrated_emb = pca_emb
            method_desc = "no_correction_baseline"

        post_mixing = calculate_batch_mixing_score(integrated_emb, batches)

        integrated_data = data.copy()
        integrated_data.obsm["X_pca"] = integrated_emb
        integrated_data.uns["integration"] = {
            "method": self.implementation_id,
            "pre_batch_mixing": pre_mixing,
            "post_batch_mixing": post_mixing,
        }

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/integrated/v3"

        registry.register(
            uri_str=out_uri,
            payload=integrated_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation=f"batch_integration_{self.implementation_id}",
            parent_uris=[in_uri],
            parameters={"method": self.implementation_id, "n_components": n_comps},
            summary_metrics={
                "n_cells": integrated_data.n_obs,
                "pre_batch_mixing": pre_mixing,
                "post_batch_mixing": post_mixing,
                "mixing_improvement": float(post_mixing - pre_mixing),
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["compute_pca", method_desc, "evaluate_batch_mixing"],
            metrics={
                "pre_batch_mixing": pre_mixing,
                "post_batch_mixing": post_mixing,
                "batch_correction_applied": self.implementation_id == "harmony",
            }
        )
