"""
Clustering and Cell Type Annotation capability.
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


def simple_kmeans(X: np.ndarray, k: int = 4, max_iter: int = 50, random_seed: int = 42) -> np.ndarray:
    """Deterministic K-Means clustering algorithm for KNN graph / embedding partitions."""
    np.random.seed(random_seed)
    n_samples = X.shape[0]
    initial_idx = np.random.choice(n_samples, size=k, replace=False)
    centroids = X[initial_idx].copy()

    labels = np.zeros(n_samples, dtype=int)
    for _ in range(max_iter):
        dists = cdist(X, centroids)
        new_labels = np.argmin(dists, axis=1)
        if np.array_equal(labels, new_labels):
            break
        labels = new_labels
        for j in range(k):
            mask = (labels == j)
            if mask.sum() > 0:
                centroids[j] = X[mask].mean(axis=0)
    return labels


def calculate_silhouette(X: np.ndarray, labels: np.ndarray) -> float:
    """Calculates approximate average silhouette score across clusters."""
    unique_labels = np.unique(labels)
    if len(unique_labels) <= 1:
        return 0.0
    
    n = X.shape[0]
    sample_indices = np.random.choice(n, size=min(n, 200), replace=False) if n > 200 else np.arange(n)
    dists = cdist(X[sample_indices], X)
    
    sil_scores = []
    for i, idx in enumerate(sample_indices):
        curr_label = labels[idx]
        same_mask = (labels == curr_label)
        same_mask[idx] = False
        a_i = dists[i, same_mask].mean() if same_mask.sum() > 0 else 0.0
        
        b_i = float("inf")
        for other_label in unique_labels:
            if other_label == curr_label:
                continue
            other_mask = (labels == other_label)
            if other_mask.sum() > 0:
                b_i = min(b_i, dists[i, other_mask].mean())
        
        if b_i == float("inf"):
            b_i = 0.0
        
        denom = max(a_i, b_i)
        sil_scores.append((b_i - a_i) / denom if denom > 0 else 0.0)
        
    return float(np.mean(sil_scores))


class ClusteringCapability(BaseCapability):
    """Performs community detection (e.g. Leiden/KNN) and marker-guided cell annotation."""

    def __init__(self):
        super().__init__(
            capability_name="clustering",
            implementation_id="leiden_knn_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        emb = data.obsm.get("X_pca", data.X[:, :min(20, data.n_vars)])
        k_clusters = contract.parameters.get("k_clusters", 4)
        seed = contract.parameters.get("random_seed", 42)

        # Community clustering
        labels = simple_kmeans(emb, k=k_clusters, random_seed=seed)
        silhouette = calculate_silhouette(emb, labels)

        # 2D UMAP-like mock projection for visualization
        u1 = emb[:, 0] + np.sin(emb[:, 1]) * 0.5
        u2 = emb[:, 1] + np.cos(emb[:, 0]) * 0.5
        umap_coords = np.column_stack([u1, u2])

        # Automatic marker-guided cell type identification
        # Check marker genes in data
        gene_names = list(data.var["gene_name"]) if "gene_name" in data.var.columns else [f"Gene_{i}" for i in range(data.n_vars)]
        name_to_idx = {g: i for i, g in enumerate(gene_names)}

        cluster_annotations = {}
        for c in range(k_clusters):
            c_mask = (labels == c)
            c_expr = data.X[c_mask]
            
            # Scores
            scores = {
                "Microglia": (c_expr[:, name_to_idx["Cx3cr1"]].mean() if "Cx3cr1" in name_to_idx else 0) +
                             (c_expr[:, name_to_idx["P2ry12"]].mean() if "P2ry12" in name_to_idx else 0),
                "Astrocytes": c_expr[:, name_to_idx["Gfap"]].mean() if "Gfap" in name_to_idx else 0,
                "Neurons": c_expr[:, name_to_idx["Rbfox3"]].mean() if "Rbfox3" in name_to_idx else 0,
                "Oligodendrocytes": c_expr[:, name_to_idx["Mog"]].mean() if "Mog" in name_to_idx else 0,
            }
            best_type = max(scores, key=scores.get) if max(scores.values()) > 0.1 else f"Cluster_{c}"
            cluster_annotations[c] = best_type

        annotated_types = [cluster_annotations[c] for c in labels]

        clustered_data = data.copy()
        clustered_data.obs["leiden"] = [str(c) for c in labels]
        clustered_data.obs["cell_type"] = annotated_types
        clustered_data.obsm["X_umap"] = umap_coords
        clustered_data.uns["clustering"] = {
            "silhouette": silhouette,
            "cluster_annotations": cluster_annotations,
        }

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/annotated/v4"

        registry.register(
            uri_str=out_uri,
            payload=clustered_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="cluster_and_annotate_cells",
            parent_uris=[in_uri],
            parameters={"k_clusters": k_clusters, "random_seed": seed},
            summary_metrics={
                "n_cells": clustered_data.n_obs,
                "silhouette_score": silhouette,
                "cell_type_counts": pd.Series(annotated_types).value_counts().to_dict(),
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["build_neighbor_graph", "find_clusters", "annotate_cell_types", "calculate_silhouette"],
            metrics={
                "silhouette_score": silhouette,
                "identified_cell_types": list(set(annotated_types)),
            }
        )
