"""
Spatial Domain Capability for EACBP.
Performs spatial coordinate validation, spatial neighborhood graph construction,
spatially regularized latent embedding, and spatial domain clustering.
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def validate_spatial_coordinates(spatial_coords: np.ndarray, expected_n_obs: int) -> np.ndarray:
    """
    Validates spatial coordinates array.
    
    Args:
        spatial_coords: Array of spatial coordinates.
        expected_n_obs: Expected number of observations/cells.
        
    Returns:
        Validated float32 numpy array of shape (N, 2) or (N, 3).
        
    Raises:
        ValueError: If coordinates are missing, non-finite, wrong shape, or degenerate.
    """
    if spatial_coords is None:
        raise ValueError("Spatial coordinates array is None.")
    
    coords = np.asarray(spatial_coords, dtype=np.float32)
    
    if coords.ndim != 2:
        raise ValueError(f"Spatial coordinates must be 2D array, got shape {coords.shape}.")
    
    n_obs, n_dims = coords.shape
    if n_obs != expected_n_obs:
        raise ValueError(
            f"Spatial coordinates cell count ({n_obs}) does not match expected observations ({expected_n_obs})."
        )
    
    if n_dims not in (2, 3):
        raise ValueError(
            f"Spatial coordinates must have 2 or 3 dimensions (x, y[, z]), found {n_dims} dimensions."
        )
    
    nan_count = int(np.isnan(coords).sum())
    inf_count = int(np.isinf(coords).sum())
    if nan_count > 0 or inf_count > 0:
        raise ValueError(
            f"Spatial coordinates contain non-finite values ({nan_count} NaNs, {inf_count} Infs)."
        )
    
    coord_var = np.var(coords, axis=0)
    if np.sum(coord_var) <= 1e-12:
        raise ValueError("Degenerate spatial coordinates: all cells have identical spatial positions (zero variance).")
    
    return coords


def build_spatial_neighborhood_graph(
    coords: np.ndarray,
    k_neighbors: int = 6,
    metric: str = "euclidean",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Constructs spatial k-NN adjacency, distance, and row-normalized weight matrices.
    
    Args:
        coords: Spatial coordinates array of shape (N, d).
        k_neighbors: Number of nearest neighbors per cell.
        metric: Distance metric for cdist.
        
    Returns:
        Tuple of (W_binary, D_pairwise, W_row_normalized).
    """
    n_obs = coords.shape[0]
    k = min(max(1, k_neighbors), max(1, n_obs - 1))
    
    # Pairwise Euclidean distances
    dists = cdist(coords, coords, metric=metric)
    np.fill_diagonal(dists, 0.0)
    
    W_binary = np.zeros((n_obs, n_obs), dtype=np.float32)
    
    for i in range(n_obs):
        neighbor_indices = np.argsort(dists[i])
        selected = [idx for idx in neighbor_indices if idx != i][:k]
        W_binary[i, selected] = 1.0
    
    # Symmetrize adjacency for undirected spatial contact
    W_sym = np.maximum(W_binary, W_binary.T)
    np.fill_diagonal(W_sym, 0.0)
    
    # Row normalization
    row_sums = W_sym.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    W_norm = W_sym / row_sums
    
    return W_sym, dists, W_norm


def compute_spatially_smoothed_embedding(
    Z: np.ndarray,
    W_norm: np.ndarray,
    smoothing_lambda: float = 0.3,
) -> np.ndarray:
    """
    Computes spatially smoothed latent embedding:
    Z_smooth = (1 - lambda) * Z + lambda * (W_norm @ Z)
    """
    smoothing_lambda = float(np.clip(smoothing_lambda, 0.0, 1.0))
    spatial_neighbor_avg = np.dot(W_norm, Z)
    Z_smooth = (1.0 - smoothing_lambda) * Z + smoothing_lambda * spatial_neighbor_avg
    return Z_smooth.astype(np.float32)


def simple_kmeans(X: np.ndarray, k: int = 4, max_iter: int = 50, random_seed: int = 42) -> np.ndarray:
    """Deterministic K-Means clustering algorithm."""
    np.random.seed(random_seed)
    n_samples = X.shape[0]
    k = min(k, n_samples)
    if k <= 1:
        return np.zeros(n_samples, dtype=int)
        
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


class SpatialDomainCapability(BaseCapability):
    """
    Spatial Domain Identification & Microenvironment Capability.
    """

    def __init__(
        self,
        capability_name: str = "spatial_domain",
        implementation_id: str = "spatial_domain_knn_v1",
    ):
        super().__init__(
            capability_name=capability_name,
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["spatial", "scRNA"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            suitable_for=["spatial_microenvironment", "spatial_domain_clustering", "tissue_segmentation"],
            output_types=[ArtifactType.SPATIAL_DATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        # Extract and validate spatial coordinates
        coords = None
        if "spatial" in data.obsm:
            coords = data.obsm["spatial"]
        elif "spatial_coords" in data.obsm:
            coords = data.obsm["spatial_coords"]
        elif "x_coord" in data.obs.columns and "y_coord" in data.obs.columns:
            if "z_coord" in data.obs.columns:
                coords = data.obs[["x_coord", "y_coord", "z_coord"]].values
            else:
                coords = data.obs[["x_coord", "y_coord"]].values
        elif "x" in data.obs.columns and "y" in data.obs.columns:
            coords = data.obs[["x", "y"]].values

        if coords is None:
            raise ValueError(
                "Spatial coordinates not found. Expected in .obsm['spatial'] or .obs[['x_coord', 'y_coord']]."
            )

        validated_coords = validate_spatial_coordinates(coords, data.n_obs)

        # Parameters
        k_neighbors = int(contract.parameters.get("k_neighbors", 6))
        n_domains = int(contract.parameters.get("n_domains", contract.parameters.get("k_clusters", 4)))
        smoothing_lambda = float(contract.parameters.get("smoothing_lambda", 0.3))
        random_seed = int(contract.parameters.get("random_seed", 42))

        # Build spatial neighborhood graph
        W_sym, D_pairwise, W_norm = build_spatial_neighborhood_graph(
            validated_coords, k_neighbors=k_neighbors
        )

        # Spatially smoothed latent embedding
        if "X_pca" in data.obsm:
            Z_base = data.obsm["X_pca"]
        else:
            X_sub = data.X[:, :min(30, data.n_vars)]
            X_mean = X_sub - X_sub.mean(axis=0)
            u, s, vt = np.linalg.svd(X_mean, full_matrices=False)
            Z_base = u * s

        Z_smooth = compute_spatially_smoothed_embedding(
            Z_base, W_norm, smoothing_lambda=smoothing_lambda
        )

        # Cluster spatial domains
        domain_labels = simple_kmeans(Z_smooth, k=n_domains, random_seed=random_seed)
        silhouette = calculate_silhouette(Z_smooth, domain_labels)

        # Build output spatial dataset
        spatial_data = data.copy()
        domain_names = [f"Domain_{d}" for d in domain_labels]
        spatial_data.obs["spatial_domain"] = domain_names
        spatial_data.obsm["spatial"] = validated_coords
        spatial_data.obsm["X_spatial_pca"] = Z_smooth
        spatial_data.obsm["spatial_connectivities"] = W_sym
        spatial_data.obsm["spatial_distances"] = D_pairwise
        
        # Attach in uns and obsp
        spatial_data.uns["spatial_connectivities"] = W_sym
        spatial_data.uns["spatial_distances"] = D_pairwise
        spatial_data.uns["spatial_domain_silhouette"] = silhouette
        spatial_data.obsp = {
            "spatial_connectivities": W_sym,
            "spatial_distances": D_pairwise,
        }

        # Resolve output URI
        uri_obj = ArtifactURI.parse(in_uri)
        if contract.expected_outputs:
            out_uri = contract.expected_outputs[0]
        else:
            out_uri = f"adata://{uri_obj.study_id}/spatial_domains/v1"

        domain_counts = pd.Series(domain_names).value_counts().to_dict()

        registry.register(
            uri_str=out_uri,
            payload=spatial_data.to_dict(),
            artifact_type=ArtifactType.SPATIAL_DATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="identify_spatial_domains",
            parent_uris=[in_uri],
            parameters={
                "k_neighbors": k_neighbors,
                "n_domains": n_domains,
                "smoothing_lambda": smoothing_lambda,
                "random_seed": random_seed,
            },
            summary_metrics={
                "n_cells": spatial_data.n_obs,
                "n_domains": n_domains,
                "silhouette_score": silhouette,
                "domain_distribution": domain_counts,
            },
        )

        all_ops = [
            "validate_spatial_coordinates",
            "build_spatial_connectivities",
            "build_spatial_knn",
            "spatially_smoothed_embedding",
            "cluster_spatial_domains",
            "calculate_silhouette",
            "identify_domains",
        ]
        if contract.allowed_operations:
            executed_ops = [op for op in all_ops if op in contract.allowed_operations]
        else:
            executed_ops = [
                "validate_spatial_coordinates",
                "build_spatial_connectivities",
                "spatially_smoothed_embedding",
                "cluster_spatial_domains",
                "calculate_silhouette",
            ]

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=executed_ops,
            metrics={
                "n_domains": n_domains,
                "silhouette_score": silhouette,
                "domain_counts": domain_counts,
            },
        )
