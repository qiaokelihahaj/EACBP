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


def compute_pca(X: Any, n_components: int = 30) -> np.ndarray:
    """Fast PCA dimensionality reduction compatible with sparse and dense matrices using scipy."""
    if hasattr(X, "tocsr"):
        X_sp = X.tocsr().astype(np.float32)
        k = min(n_components, max(1, X_sp.shape[0] - 2), max(1, X_sp.shape[1] - 2))
        try:
            from scipy.sparse.linalg import svds
            u, s, vt = svds(X_sp, k=k)
            return u[:, ::-1] * s[::-1]
        except Exception:
            X_dense = X_sp.toarray()
            X_centered = X_dense - np.mean(X_dense, axis=0, keepdims=True)
            u, s, vt = np.linalg.svd(X_centered, full_matrices=False)
            return u[:, :k] * s[:k]
    else:
        X_dense = np.asarray(X, dtype=np.float32)
        X_centered = X_dense - np.mean(X_dense, axis=0, keepdims=True)
        u, s, vt = np.linalg.svd(X_centered, full_matrices=False)
        k = min(n_components, X_dense.shape[0], X_dense.shape[1])
        return u[:, :k] * s[:k]


def calculate_batch_mixing_score(
    embedding: np.ndarray,
    batches: np.ndarray,
    k: int = 15,
    random_seed: int = 42,
) -> float:
    """Calculates average fraction of nearest neighbors from different batches."""
    if len(np.unique(batches)) <= 1 or embedding.shape[0] <= 1:
        return 1.0
    k = min(max(1, int(k)), embedding.shape[0] - 1)
    
    # Subsample if too large for speed
    n = embedding.shape[0]
    rng = np.random.default_rng(random_seed)
    indices = rng.choice(n, size=min(n, 300), replace=False) if n > 300 else np.arange(n)
    dists = cdist(embedding[indices], embedding)
    
    mixing_scores = []
    for i, idx in enumerate(indices):
        nn_indices = np.argsort(dists[i])[1:k+1]
        same_batch_count = (batches[nn_indices] == batches[idx]).sum()
        # Mixing score is higher when neighbors are diverse
        mixing_scores.append(1.0 - (same_batch_count / k))
    
    return float(np.mean(mixing_scores))


class IntegrationCapability(BaseCapability):
    """Evaluate batch mixing and apply deterministic batch mean centering.

    The previous implementation exposed the mean-shift calculation as
    ``harmony``.  Harmony is an iterative soft-clustering algorithm and is
    not implemented here.  ``harmony`` is therefore accepted only as a
    compatibility request alias; ``method_used`` is always the honest
    ``batch_mean_centering_v1`` identifier.
    """

    def __init__(self, implementation_id: str = "batch_mean_centering_v1"):
        requested_id = implementation_id
        canonical_id = "no_correction_v1" if implementation_id in {"no_correction", "no_correction_v1"} else "batch_mean_centering_v1"
        super().__init__(
            capability_name="integration",
            implementation_id=canonical_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )
        self.requested_implementation_id = requested_id
        self.legacy_aliases = {
            "harmony": "batch_mean_centering_v1",
            "harmony_v1": "batch_mean_centering_v1",
            "no_correction": "no_correction_v1",
        }

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

        seed = int(contract.parameters.get("random_seed", 42))
        pre_mixing = calculate_batch_mixing_score(pca_emb, batches, random_seed=seed)

        if self.implementation_id == "batch_mean_centering_v1" and batch_col:
            # Per-batch mean centering in PCA space.  This is deliberately a
            # simple deterministic correction, not Harmony.
            adjusted_pca = pca_emb.copy()
            unique_batches = np.unique(batches)
            global_mean = np.mean(adjusted_pca, axis=0, keepdims=True)
            for b in unique_batches:
                mask = (batches == b)
                batch_mean = np.mean(adjusted_pca[mask], axis=0, keepdims=True)
                adjusted_pca[mask] -= 0.8 * (batch_mean - global_mean)
            integrated_emb = adjusted_pca
            method_desc = "batch_mean_centering_pca"
        else:
            integrated_emb = pca_emb
            method_desc = "no_correction_baseline"

        post_mixing = calculate_batch_mixing_score(integrated_emb, batches, random_seed=seed)

        integrated_data = data.copy()
        integrated_data.obsm["X_pca"] = integrated_emb
        integrated_data.uns["integration"] = {
            "method": self.implementation_id,
            "requested_method": self.requested_implementation_id,
            "algorithm": "per_batch_pca_mean_centering" if self.implementation_id == "batch_mean_centering_v1" else "identity",
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
            parameters={"method": self.implementation_id, "requested_method": self.requested_implementation_id, "n_components": n_comps, "random_seed": seed},
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
                "batch_correction_applied": self.implementation_id == "batch_mean_centering_v1",
            }
        )
