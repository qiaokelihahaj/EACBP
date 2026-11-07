"""Spatial kNN graph and domain capability.

Neighbourhood construction uses :class:`scipy.spatial.cKDTree` and stores
only kNN edges for large datasets.  For backwards compatibility small
inputs still receive dense arrays, while large inputs receive CSR matrices;
consumers must therefore use sparse-safe ``dot``/row-sum operations.
"""

from typing import Any, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import sparse

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.clustering import simple_kmeans, calculate_silhouette as _calculate_silhouette
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.biological_units import resolve_donor_role
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


DENSE_GRAPH_LIMIT = 1500


def calculate_silhouette(X: np.ndarray, labels: np.ndarray, random_seed: int = 42) -> float:
    """Return the shared silhouette score with the spatial float32 contract.

    The general clustering helper accepts the caller's numeric dtype.  The
    historical spatial helper converted its embedding to ``float32`` before
    calculating distances; retain that small public compatibility detail
    while keeping the actual implementation in :mod:`capabilities.clustering`.
    """

    return _calculate_silhouette(np.asarray(X, dtype=np.float32), labels, random_seed=random_seed)


def validate_spatial_coordinates(spatial_coords: np.ndarray, expected_n_obs: int) -> np.ndarray:
    """Validate finite 2D/3D coordinates with non-zero variance."""
    if spatial_coords is None:
        raise ValueError("Spatial coordinates array is None.")
    coords = np.asarray(spatial_coords, dtype=np.float32)
    if coords.ndim != 2:
        raise ValueError(f"Spatial coordinates must be 2D array, got shape {coords.shape}.")
    n_obs, n_dims = coords.shape
    if n_obs != expected_n_obs:
        raise ValueError(f"Spatial coordinates cell count ({n_obs}) does not match expected observations ({expected_n_obs}).")
    if n_dims not in (2, 3):
        raise ValueError(f"Spatial coordinates must have 2 or 3 dimensions (x, y[, z]), found {n_dims} dimensions.")
    if not np.isfinite(coords).all():
        raise ValueError("Spatial coordinates contain non-finite values (NaNs/Infs).")
    if float(np.sum(np.var(coords, axis=0))) <= 1e-12:
        raise ValueError("Degenerate spatial coordinates: all cells have identical spatial positions (zero variance).")
    return coords


def extract_spatial_coordinates(data: SCData) -> np.ndarray:
    """Require measured coordinates; cell order cannot supply a spatial layout."""
    for key in ("spatial", "spatial_coords"):
        if key in data.obsm:
            return validate_spatial_coordinates(data.obsm[key], data.n_obs)
    for x_col, y_col, z_col in (("x_coord", "y_coord", "z_coord"),
                              ("spatial_x", "spatial_y", "spatial_z"),
                              ("x", "y", "z"), ("X", "Y", "Z")):
        if x_col in data.obs and y_col in data.obs:
            columns = [x_col, y_col] + ([z_col] if z_col in data.obs else [])
            return validate_spatial_coordinates(data.obs[columns].to_numpy(), data.n_obs)
    raise ValueError("Measured spatial coordinates not found; a spatial conclusion requires observed coordinates.")


def _spatial_partition_labels(values: Any, n_obs: int) -> np.ndarray:
    if values is None:
        return np.zeros(n_obs, dtype=int)
    series = pd.Series(values)
    if len(series) != n_obs:
        raise ValueError("Spatial section labels must align with observations.")
    labels = series.astype("string").str.strip()
    if (labels.isna().any() or labels.eq("").any()
            or labels.str.casefold().isin({"nan", "none", "null", "<na>"}).any()):
        raise ValueError("Spatial section labels must be complete and nonempty.")
    return pd.factorize(labels, sort=False)[0]


def resolve_spatial_sections(data: SCData, parameters: Optional[Mapping[str, Any]] = None) -> Tuple[np.ndarray, dict]:
    """Resolve observed section identity, keeping donor-local section IDs separate.

    Donors are biological units, not inferred section identities. Multi-donor
    inputs without section metadata therefore cannot define spatial adjacency.
    """
    params = parameters or {}
    observed_roles = data.uns.get("autonomous_input_preparation", {}).get("obs_column_map", {})
    role_map = params.get("obs_column_map", observed_roles)
    if not isinstance(role_map, Mapping):
        raise ValueError("Spatial observation role map must be a mapping.")
    aliases = {
        "section": ("section", "section_id", "slide_id", "slide", "slice_id", "slice"),
    }
    columns = {}
    codes = {}
    donor_params = dict(params)
    if not donor_params.get("donor_col") and role_map.get("donor"):
        donor_params["donor_col"] = role_map["donor"]
    donor_role = resolve_donor_role(data.obs, donor_params)
    columns["donor"] = donor_role.column
    if donor_role.column:
        codes["donor"] = _spatial_partition_labels(data.obs[donor_role.column], data.n_obs)
    for role, names in aliases.items():
        explicit = params.get(f"{role}_col") or role_map.get(role)
        candidates = [explicit] if explicit else [name for name in names if name in data.obs]
        if explicit and explicit not in data.obs:
            raise ValueError(f"Declared spatial {role} column {explicit!r} was not observed.")
        if not candidates:
            columns[role] = None
            continue
        selected = candidates[0]
        selected_labels = data.obs[selected].astype("string").str.strip()
        codes[role] = _spatial_partition_labels(data.obs[selected], data.n_obs)
        for alternate in candidates[1:]:
            if not selected_labels.equals(data.obs[alternate].astype("string").str.strip()):
                raise ValueError(f"Ambiguous spatial {role} columns {selected!r} and {alternate!r} disagree.")
        columns[role] = selected
    if not columns["section"] and "donor" in codes and len(np.unique(codes["donor"])) > 1:
        raise ValueError("Multiple donors require observed section identities; a donor cannot be assumed to be one section.")
    identities = [tuple(codes[role][i] for role in ("donor", "section") if role in codes)
                  for i in range(data.n_obs)]
    sections = pd.factorize(pd.Series(identities, dtype=object), sort=False)[0] if codes else np.zeros(data.n_obs, dtype=int)
    return sections, {"section_col": columns["section"], **donor_role.as_metrics(),
                      "n_spatial_partitions": int(len(np.unique(sections))),
                      "identity_basis": "observed_donor_and_section" if columns["section"] else "single_coordinate_system"}


def validate_spatial_graph(W: Any, n_obs: int, section_labels: Any = None) -> None:
    """Reject invalid supplied contacts, including edges between sections."""
    if getattr(W, "shape", None) != (n_obs, n_obs):
        raise ValueError("Spatial connectivity graph must be observation-aligned and square.")
    sections = _spatial_partition_labels(section_labels, n_obs)
    if sparse.issparse(W):
        edges = W.tocoo()
        values = edges.data
        rows, cols = edges.row, edges.col
    else:
        matrix = np.asarray(W)
        values = matrix.ravel()
        rows, cols = np.nonzero(matrix)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Spatial connectivity weights must be finite and nonnegative.")
    if sparse.issparse(W):
        nonzero = values != 0
        rows, cols = rows[nonzero], cols[nonzero]
    if np.any(rows == cols):
        raise ValueError("Spatial connectivity graph must not include self contacts.")
    if np.any(sections[rows] != sections[cols]):
        raise ValueError("Spatial connectivity graph contains edges across observed sections.")


def build_spatial_knn_indices(coords: np.ndarray, k_neighbors: int = 6,
                              metric: str = "euclidean", section_labels: Any = None) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """Query neighbors inside each measured coordinate system, including small sections."""
    coords = np.asarray(coords, dtype=np.float32)
    if coords.ndim != 2 or not np.isfinite(coords).all():
        raise ValueError("coords must be a finite 2D array")
    n_obs = len(coords)
    sections = _spatial_partition_labels(section_labels, n_obs)
    neighbors = [np.array([], dtype=int) for _ in range(n_obs)]
    distances = [np.array([], dtype=np.float32) for _ in range(n_obs)]
    for section in np.unique(sections):
        positions = np.flatnonzero(sections == section)
        if len(positions) < 2:
            continue
        k = min(max(1, int(k_neighbors)), len(positions) - 1)
        tree = cKDTree(coords[positions])
        try:
            edge_distances, indices = tree.query(coords[positions], k=k + 1, p=_tree_p(metric), workers=1)
        except TypeError:  # older SciPy lacks workers
            edge_distances, indices = tree.query(coords[positions], k=k + 1, p=_tree_p(metric))
        for local, position in enumerate(positions):
            keep = np.asarray(indices[local]) != local
            neighbors[position] = positions[np.asarray(indices[local])[keep][:k]]
            distances[position] = np.asarray(edge_distances[local])[keep][:k].astype(np.float32)
    return neighbors, distances


def _tree_p(metric: str):
    if metric in {"euclidean", "l2"}:
        return 2
    if metric in {"cityblock", "manhattan", "l1"}:
        return 1
    if metric in {"chebyshev", "linf", "infinity"}:
        return np.inf
    raise ValueError(f"Metric '{metric}' is unsupported by sparse cKDTree graph construction")


def build_spatial_neighborhood_graph(
    coords: np.ndarray,
    k_neighbors: int = 6,
    metric: str = "euclidean",
    section_labels: Any = None,
) -> Tuple[Any, Any, Any]:
    """Build a symmetric binary kNN graph without an all-pairs distance matrix.

    For ``N <= DENSE_GRAPH_LIMIT`` the return values are dense for compatibility
    with existing callers.  Above that limit all three outputs are CSR sparse
    matrices and contain only graph edges (O(N*k) storage).
    """
    coords = np.asarray(coords, dtype=np.float32)
    if coords.ndim != 2:
        raise ValueError("coords must be a 2D array")
    n_obs = coords.shape[0]
    if n_obs == 0:
        empty = np.zeros((0, 0), dtype=np.float32)
        return empty, empty.copy(), empty.copy()
    neighbors, distances = build_spatial_knn_indices(coords, k_neighbors, metric, section_labels)
    rows = np.repeat(np.arange(n_obs), [len(indices) for indices in neighbors])
    nbrs = np.concatenate(neighbors)
    edge_dist = np.concatenate(distances)

    directed = sparse.coo_matrix((np.ones(len(rows), dtype=np.float32), (rows, nbrs)), shape=(n_obs, n_obs)).tocsr()
    directed_dist = sparse.coo_matrix((edge_dist, (rows, nbrs)), shape=(n_obs, n_obs)).tocsr()
    W_sparse = directed.maximum(directed.T).tocsr()
    # Distances are stored only on graph edges.  ``maximum`` recovers
    # one-sided edges without constructing a dense distance matrix; for
    # reciprocal edges both directions represent the same metric distance.
    D_sparse = directed_dist.maximum(directed_dist.T).tocsr()
    D_sparse.setdiag(0.0)
    D_sparse.eliminate_zeros()
    row_sums = np.asarray(W_sparse.sum(axis=1)).ravel()
    safe = np.where(row_sums > 0, row_sums, 1.0)
    W_norm_sparse = sparse.diags(1.0 / safe).dot(W_sparse).tocsr()

    if n_obs <= DENSE_GRAPH_LIMIT:
        W_dense = W_sparse.toarray().astype(np.float32)
        D_dense = D_sparse.toarray().astype(np.float32)
        W_norm_dense = W_norm_sparse.toarray().astype(np.float32)
        np.fill_diagonal(W_dense, 0.0)
        np.fill_diagonal(D_dense, 0.0)
        return W_dense, D_dense, W_norm_dense
    return W_sparse.astype(np.float32), D_sparse.astype(np.float32), W_norm_sparse.astype(np.float32)


def compute_spatially_smoothed_embedding(Z: np.ndarray, W_norm, smoothing_lambda: float = 0.3) -> np.ndarray:
    """Compute ``(1-lambda) Z + lambda W_norm Z`` for dense or sparse W."""
    smoothing_lambda = float(np.clip(smoothing_lambda, 0.0, 1.0))
    Z = np.asarray(Z, dtype=np.float32)
    neighbor_avg = W_norm.dot(Z) if hasattr(W_norm, "dot") else np.dot(W_norm, Z)
    return ((1.0 - smoothing_lambda) * Z + smoothing_lambda * np.asarray(neighbor_avg)).astype(np.float32)


class SpatialDomainCapability(BaseCapability):
    """Spatial kNN smoothing followed by K-means domain partitioning."""

    def __init__(self, capability_name: str = "spatial_domain", implementation_id: str = "spatial_knn_kmeans_v1"):
        super().__init__(
            capability_name=capability_name,
            implementation_id="spatial_knn_kmeans_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["spatial", "scRNA"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            suitable_for=["spatial_microenvironment", "spatial_domain_clustering", "tissue_segmentation"],
            output_types=[ArtifactType.SPATIAL_DATA],
        )
        self.requested_implementation_id = implementation_id
        self.legacy_aliases = {"spatial_domain_knn_v1": self.implementation_id}

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        _, data = registry.get_scdata(in_uri)
        validated_coords = extract_spatial_coordinates(data)

        params = contract.parameters
        k_neighbors = int(params.get("k_neighbors", 6))
        requested_domains = int(params.get("n_domains", params.get("k_clusters", 4)))
        n_domains = min(max(1, requested_domains), max(1, data.n_obs))
        smoothing_lambda = float(params.get("smoothing_lambda", 0.3))
        random_seed = int(params.get("random_seed", 42))
        sections, section_receipt = resolve_spatial_sections(data, params)
        W_sym, D_pairwise, W_norm = build_spatial_neighborhood_graph(validated_coords, k_neighbors=k_neighbors, section_labels=sections)

        if "X_pca" in data.obsm:
            Z_base = np.asarray(data.obsm["X_pca"], dtype=np.float32)
        else:
            X_sub = data.X[:, : min(30, data.n_vars)]
            if hasattr(X_sub, "toarray"):
                X_sub = X_sub.toarray()
            X_sub = np.asarray(X_sub, dtype=np.float32)
            X_centered = X_sub - X_sub.mean(axis=0)
            u, s, _ = np.linalg.svd(X_centered, full_matrices=False)
            Z_base = u * s
        Z_smooth = compute_spatially_smoothed_embedding(Z_base, W_norm, smoothing_lambda=smoothing_lambda)
        domain_labels = simple_kmeans(Z_smooth, k=n_domains, random_seed=random_seed)
        silhouette = calculate_silhouette(Z_smooth, domain_labels, random_seed=random_seed)

        spatial_data = data.copy()
        domain_names = [f"Domain_{d}" for d in domain_labels]
        spatial_data.obs["spatial_domain"] = domain_names
        spatial_data.obsm["spatial"] = validated_coords
        spatial_data.obsm["X_spatial_pca"] = Z_smooth
        # Keep small-input legacy locations; new SCData also preserves obsp.
        spatial_data.obsm["spatial_connectivities"] = W_sym
        spatial_data.obsm["spatial_distances"] = D_pairwise
        spatial_data.uns["spatial_connectivities"] = W_sym
        spatial_data.uns["spatial_distances"] = D_pairwise
        spatial_data.uns["spatial_domain_silhouette"] = silhouette
        spatial_data.uns["spatial_partition_receipt"] = section_receipt
        if not hasattr(spatial_data, "obsp") or getattr(spatial_data, "obsp", None) is None:
            spatial_data.obsp = {}
        spatial_data.obsp["spatial_connectivities"] = W_sym
        spatial_data.obsp["spatial_distances"] = D_pairwise

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = contract.expected_outputs[0] if contract.expected_outputs else f"adata://{uri_obj.study_id}/spatial_domains/v1"
        domain_counts = pd.Series(domain_names).value_counts().to_dict()
        registry.register(
            uri_str=out_uri,
            payload=spatial_data.to_dict(),
            artifact_type=ArtifactType.SPATIAL_DATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="spatial_knn_smooth_and_kmeans_domains",
            parent_uris=[in_uri],
            parameters={"k_neighbors": k_neighbors, "n_domains": n_domains, "smoothing_lambda": smoothing_lambda, "random_seed": random_seed, "graph_storage": "csr_for_large_inputs", **section_receipt},
            summary_metrics={"n_cells": spatial_data.n_obs, "n_domains": len(domain_counts), "silhouette_score": silhouette, "domain_distribution": domain_counts},
        )
        all_ops = ["validate_spatial_coordinates", "build_spatial_knn", "build_spatial_connectivities", "spatially_smoothed_embedding", "cluster_spatial_domains", "calculate_silhouette"]
        executed_ops = [op for op in all_ops if not contract.allowed_operations or op in contract.allowed_operations]
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=executed_ops,
            metrics={"n_domains": len(domain_counts), "silhouette_score": silhouette, "domain_counts": domain_counts, "graph_storage": "sparse" if sparse.issparse(W_sym) else "dense_compatibility", **section_receipt},
        )
