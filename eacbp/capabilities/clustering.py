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
    """Run a small, deterministic Lloyd K-means implementation.

    This helper is deliberately named after the algorithm that is actually
    implemented.  It does not construct a KNN graph and must not be reported
    as Leiden (the old public ``leiden_knn_v1`` name is handled as a caller
    compatibility alias by the registry).
    """
    X = np.asarray(X, dtype=np.float32)
    n_samples = X.shape[0]
    if n_samples == 0:
        return np.zeros(0, dtype=int)
    k = min(max(1, int(k)), n_samples)
    rng = np.random.default_rng(random_seed)
    initial_idx = rng.choice(n_samples, size=k, replace=False)
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


def calculate_silhouette(X: np.ndarray, labels: np.ndarray, random_seed: int = 42) -> float:
    """Calculates approximate average silhouette score across clusters."""
    unique_labels = np.unique(labels)
    if len(unique_labels) <= 1:
        return 0.0
    
    n = X.shape[0]
    rng = np.random.default_rng(random_seed)
    sample_indices = rng.choice(n, size=min(n, 200), replace=False) if n > 200 else np.arange(n)
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
    """K-means clustering plus marker-score annotation.

    Earlier releases exposed this simplified calculation as
    ``leiden_knn_v1`` and wrote ``leiden``/``X_umap`` fields.  Those names
    were scientifically misleading.  The constructor accepts that old ID as
    a compatibility alias, while ``implementation_id`` and ``method_used``
    always describe the calculation performed here.
    """

    def __init__(self, implementation_id: str = "kmeans_marker_embedding_v1"):
        super().__init__(
            capability_name="clustering",
            # ``leiden_knn_v1`` remains an accepted request alias in the
            # registry, but this object reports the real implementation.
            implementation_id="kmeans_marker_embedding_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )
        self.requested_implementation_id = implementation_id
        self.legacy_aliases = {"leiden_knn_v1": self.implementation_id}

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        emb = data.obsm.get("X_pca", None)
        if emb is None:
            emb = data.X[:, :min(20, data.n_vars)]
            if hasattr(emb, "toarray"):
                emb = emb.toarray()
        else:
            if hasattr(emb, "toarray"):
                emb = emb.toarray()
            emb = np.asarray(emb, dtype=np.float32)

        k_clusters = contract.parameters.get("k_clusters", 4)
        seed = contract.parameters.get("random_seed", 42)

        # This is K-means; no graph/community algorithm is involved.
        labels = simple_kmeans(emb, k=k_clusters, random_seed=seed)
        silhouette = calculate_silhouette(emb, labels, random_seed=seed)

        # A deterministic two-coordinate display projection.  It is not UMAP.
        u1 = emb[:, 0] + np.sin(emb[:, 1] if emb.shape[1] > 1 else emb[:, 0]) * 0.5
        u2 = emb[:, 1] if emb.shape[1] > 1 else emb[:, 0] + np.cos(emb[:, 0]) * 0.5
        display_coords = np.column_stack([u1, u2]).astype(np.float32)

        # Automatic marker-guided cell type identification
        gene_names = list(data.var["gene_name"]) if "gene_name" in data.var.columns else [f"Gene_{i}" for i in range(data.n_vars)]
        name_to_idx = {g: i for i, g in enumerate(gene_names)}

        cluster_annotations = {}
        for c in range(int(np.max(labels)) + 1 if len(labels) else 0):
            c_mask = (labels == c)
            c_expr = data.X[c_mask]

            def get_mean_expr(gene_sym: str) -> float:
                if gene_sym in name_to_idx:
                    return float(c_expr[:, name_to_idx[gene_sym]].mean())
                # Try capitalized / upper / lower
                for cand in [gene_sym.capitalize(), gene_sym.upper(), gene_sym.lower()]:
                    if cand in name_to_idx:
                        return float(c_expr[:, name_to_idx[cand]].mean())
                return 0.0

            scores = {
                "Microglia": get_mean_expr("Cx3cr1") + get_mean_expr("P2ry12") + get_mean_expr("Tmem119"),
                "Progenitors": get_mean_expr("Sox2") + get_mean_expr("Cdk1") + get_mean_expr("Top2a"),
                "Immature_Neurons": get_mean_expr("Dcx") + get_mean_expr("Tubb3"),
                "Neurons": get_mean_expr("Rbfox3") + get_mean_expr("Syt1") + get_mean_expr("Snap25"),
                "Astrocytes": get_mean_expr("Gfap") + get_mean_expr("Aqp4"),
                "Oligodendrocytes": get_mean_expr("Mog") + get_mean_expr("Mbp"),
            }
            best_type = max(scores, key=scores.get) if max(scores.values()) > 0.05 else f"Cluster_{c}"
            cluster_annotations[c] = best_type

        # Never use a ground-truth annotation as a prediction.  The marker
        # scores above are the only source of the inferred labels.
        annotated_types = [cluster_annotations[c] for c in labels]

        clustered_data = data.copy()
        clustered_data.obs["cluster"] = [str(c) for c in labels]
        clustered_data.obs["cell_type"] = annotated_types
        clustered_data.obsm["X_embedding_2d"] = display_coords
        clustered_data.uns["clustering"] = {
            "method": self.implementation_id,
            "cluster_algorithm": "lloyd_kmeans",
            "embedding_method": "pca_first_two_with_sine_display",
            "silhouette": silhouette,
            # HDF5/AnnData requires string mapping keys.
            "cluster_annotations": {str(k): v for k, v in cluster_annotations.items()},
        }

        # Legacy columns are opt-in and explicitly marked.  This lets old
        # consumers migrate without silently presenting K-means as Leiden or
        # the display projection as UMAP.
        if bool(contract.parameters.get("emit_legacy_aliases", False)):
            clustered_data.obs["leiden"] = clustered_data.obs["cluster"]
            clustered_data.obsm["X_umap"] = clustered_data.obsm["X_embedding_2d"]
            clustered_data.uns["clustering"]["legacy_aliases"] = {
                "leiden": "cluster",
                "X_umap": "X_embedding_2d",
            }

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/annotated/v4"

        registry.register(
            uri_str=out_uri,
            payload=clustered_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="kmeans_cluster_and_marker_annotate_cells",
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
            executed_operations=["pca_first_two_display", "lloyd_kmeans", "marker_score_annotation", "calculate_silhouette"],
            metrics={
                "silhouette_score": silhouette,
                "identified_cell_types": list(set(annotated_types)),
            }
        )
