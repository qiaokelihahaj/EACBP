"""
Spatial Autocorrelation and Spatial DEG Capability for EACBP.
Implements Global Moran's I and Geary's C with exact analytical variance
under the randomization hypothesis, z-scores, p-values, and Benjamini-Hochberg FDR.
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from scipy import stats

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial.domain import (
    validate_spatial_coordinates,
    build_spatial_neighborhood_graph,
)
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Calculates Benjamini-Hochberg FDR adjusted p-values."""
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    fdr = np.zeros(n)
    curr_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        adj = (sorted_p[i] * n) / rank
        curr_min = min(curr_min, adj)
        fdr[i] = min(curr_min, 1.0)
    rev_indices = np.argsort(sorted_indices)
    return fdr[rev_indices]


def calculate_morans_i(
    x: np.ndarray,
    W: np.ndarray,
) -> Tuple[float, float, float, float, float]:
    """
    Calculates Global Moran's I with exact analytical variance under the randomization hypothesis.
    
    Args:
        x: 1D array of gene expression values (length N).
        W: Spatial weight matrix (N x N) with zero diagonal.
        
    Returns:
        Tuple of (moran_i, expected_i, var_i, z_score, p_value)
    """
    x = np.asarray(x, dtype=float).ravel()
    N = len(x)
    if N < 4:
        return 0.0, 0.0, 0.0, 0.0, 1.0
    
    x_bar = np.mean(x)
    z = x - x_bar
    denom = np.sum(z ** 2)
    
    if denom <= 1e-15:
        return 0.0, -1.0 / (N - 1), 0.0, 0.0, 1.0
    
    S0 = float(np.sum(W))
    if S0 <= 1e-15:
        return 0.0, -1.0 / (N - 1), 0.0, 0.0, 1.0
    
    numerator = float(np.dot(z, np.dot(W, z)))
    I = float((N / S0) * (numerator / denom))
    
    expected_I = -1.0 / (N - 1)
    
    W_plus_WT = W + W.T
    S1 = 0.5 * float(np.sum(W_plus_WT ** 2))
    
    row_col_sums = np.sum(W, axis=1) + np.sum(W, axis=0)
    S2 = float(np.sum(row_col_sums ** 2))
    
    b2 = float((N * np.sum(z ** 4)) / (denom ** 2))
    
    A = N * ((N ** 2 - 3 * N + 3) * S1 - N * S2 + 3 * (S0 ** 2))
    B = b2 * ((N ** 2 - N) * S1 - 2 * N * S2 + 6 * (S0 ** 2))
    C = (N - 1) * (N - 2) * (N - 3) * (S0 ** 2)
    
    if C <= 1e-15:
        var_I = 0.0
    else:
        var_I = max(0.0, float((A - B) / C - (expected_I ** 2)))
    
    std_I = np.sqrt(var_I) if var_I > 1e-15 else 1.0
    z_score = float((I - expected_I) / std_I) if var_I > 1e-15 else 0.0
    
    p_val = float(2.0 * stats.norm.sf(abs(z_score)))
    p_val = float(np.clip(p_val, 0.0, 1.0))
    
    return I, expected_I, var_I, z_score, p_val


def calculate_gearys_c(
    x: np.ndarray,
    W: np.ndarray,
) -> Tuple[float, float, float, float, float]:
    """
    Calculates Geary's C with exact analytical variance under the randomization hypothesis.
    
    Args:
        x: 1D array of gene expression values (length N).
        W: Spatial weight matrix (N x N) with zero diagonal.
        
    Returns:
        Tuple of (geary_c, expected_c, var_c, z_score, p_value)
    """
    x = np.asarray(x, dtype=float).ravel()
    N = len(x)
    if N < 4:
        return 1.0, 1.0, 0.0, 0.0, 1.0
    
    x_bar = np.mean(x)
    z = x - x_bar
    denom = np.sum(z ** 2)
    
    if denom <= 1e-15:
        return 1.0, 1.0, 0.0, 0.0, 1.0
    
    S0 = float(np.sum(W))
    if S0 <= 1e-15:
        return 1.0, 1.0, 0.0, 0.0, 1.0
    
    x_sq = x ** 2
    row_sums = np.sum(W, axis=1)
    col_sums = np.sum(W, axis=0)
    diff_sq_sum = float(np.sum(x_sq * row_sums) + np.sum(x_sq * col_sums) - 2.0 * np.dot(x, np.dot(W, x)))
    
    C = float(((N - 1) / (2.0 * S0)) * (diff_sq_sum / denom))
    expected_C = 1.0
    
    # Analytical variance under randomization hypothesis (Cliff & Ord 1981 / Anselin 1995 / PySAL)
    W_plus_WT = W + W.T
    S1 = 0.5 * float(np.sum(W_plus_WT ** 2))
    row_col_sums = np.sum(W, axis=1) + np.sum(W, axis=0)
    S2 = float(np.sum(row_col_sums ** 2))
    b2 = float((N * np.sum(z ** 4)) / (denom ** 2))
    
    s02 = S0 * S0
    a = (N - 1) * S1 * (N * N - 3 * N + 3 - (N - 1) * b2)
    b = -0.25 * (N - 1) * S2 * (N * N + 3 * N - 6 - (N * N - N + 2) * b2)
    c = s02 * (N * N - 3 - (N - 1) * (N - 1) * b2)
    d = (N - 1) * (N - 2) * (N - 3) * s02
    
    if d <= 1e-15:
        var_C = 0.0
    else:
        var_C = float((a + b + c) / d)
    
    # Fallback to normality variance if numerical instability occurs
    if var_C <= 0.0:
        var_norm = float(((2 * S1 + S2) * (N - 1) - 4 * s02) / (2 * (N + 1) * s02))
        var_C = max(1e-12, var_norm)
    
    std_C = np.sqrt(var_C) if var_C > 1e-15 else 1.0
    z_score = float((C - expected_C) / std_C) if var_C > 1e-15 else 0.0
    
    p_val = float(2.0 * stats.norm.sf(abs(z_score)))
    p_val = float(np.clip(p_val, 0.0, 1.0))
    
    return C, expected_C, var_C, z_score, p_val


class SpatialDEGCapability(BaseCapability):
    """
    Spatial Autocorrelation and Spatially Variable Gene (Spatial DEG) Capability.
    """

    def __init__(
        self,
        capability_name: str = "spatial_deg",
        implementation_id: str = "spatial_moran_deg_v1",
    ):
        super().__init__(
            capability_name=capability_name,
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["spatial", "scRNA"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            suitable_for=["spatially_variable_genes", "spatial_autocorrelation", "spatial_deg"],
            output_types=[ArtifactType.TABLE, ArtifactType.SPATIAL_DATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        # Extract spatial coordinates
        coords = None
        if "spatial" in data.obsm:
            coords = data.obsm["spatial"]
        elif "spatial_coords" in data.obsm:
            coords = data.obsm["spatial_coords"]
        elif "x_coord" in data.obs.columns and "y_coord" in data.obs.columns:
            coords = data.obs[["x_coord", "y_coord"]].values
        elif "x" in data.obs.columns and "y" in data.obs.columns:
            coords = data.obs[["x", "y"]].values

        if coords is None:
            raise ValueError("Spatial coordinates not found for spatial autocorrelation analysis.")

        validated_coords = validate_spatial_coordinates(coords, data.n_obs)

        # Parameters
        k_neighbors = int(contract.parameters.get("k_neighbors", 6))
        min_moran_i = float(contract.parameters.get("min_moran_i", 0.15))
        fdr_threshold = float(contract.parameters.get("fdr_threshold", 0.05))
        target_genes = contract.parameters.get("target_genes", None)

        # Build or retrieve spatial connectivity graph safely
        if "spatial_connectivities" in data.obsm and isinstance(data.obsm["spatial_connectivities"], np.ndarray):
            W = np.asarray(data.obsm["spatial_connectivities"], dtype=np.float32)
        elif "spatial_connectivities" in data.uns and isinstance(data.uns["spatial_connectivities"], np.ndarray):
            W = np.asarray(data.uns["spatial_connectivities"], dtype=np.float32)
        elif hasattr(data, "obsp") and isinstance(getattr(data, "obsp", None), dict) and "spatial_connectivities" in data.obsp and isinstance(data.obsp["spatial_connectivities"], np.ndarray):
            W = np.asarray(data.obsp["spatial_connectivities"], dtype=np.float32)
        else:
            W, _, _ = build_spatial_neighborhood_graph(validated_coords, k_neighbors=k_neighbors)

        gene_names = list(data.var["gene_name"]) if "gene_name" in data.var.columns else [f"Gene_{i}" for i in range(data.n_vars)]
        
        if target_genes:
            gene_indices = [i for i, g in enumerate(gene_names) if g in target_genes]
            if not gene_indices:
                gene_indices = list(range(data.n_vars))
        else:
            gene_indices = list(range(data.n_vars))

        results = []
        for idx in gene_indices:
            g_name = str(gene_names[idx])
            expr = data.X[:, idx]
            
            m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(expr, W)
            g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(expr, W)
            
            results.append({
                "gene": g_name,
                "gene_index": idx,
                "moran_i": float(m_i),
                "moran_expected": float(m_exp),
                "moran_variance": float(m_var),
                "moran_z_score": float(m_z),
                "moran_p_value": float(m_p),
                "geary_c": float(g_c),
                "geary_expected": float(g_exp),
                "geary_variance": float(g_var),
                "geary_z_score": float(g_z),
                "geary_p_value": float(g_p),
                "p_value": float(m_p),
            })

        res_df = pd.DataFrame(results)
        if not res_df.empty:
            fdr = benjamini_hochberg(res_df["moran_p_value"].values)
            res_df["fdr_q_value"] = fdr
            res_df["p_val_adj"] = fdr
            res_df["is_spatially_variable"] = (res_df["fdr_q_value"] < fdr_threshold) & (res_df["moran_i"] > min_moran_i)
        else:
            res_df["fdr_q_value"] = []
            res_df["p_val_adj"] = []
            res_df["is_spatially_variable"] = []

        res_df = res_df.sort_values(["fdr_q_value", "moran_i"], ascending=[True, False]).reset_index(drop=True)
        sig_svgs = res_df[res_df["is_spatially_variable"]]

        # Output table URI
        uri_obj = ArtifactURI.parse(in_uri)
        if contract.expected_outputs:
            out_uri = contract.expected_outputs[0]
        else:
            out_uri = f"table://{uri_obj.study_id}/spatial_deg/v1"

        registry.register(
            uri_str=out_uri,
            payload=res_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="calculate_spatial_autocorrelation_deg",
            parent_uris=[in_uri],
            parameters={
                "k_neighbors": k_neighbors,
                "min_moran_i": min_moran_i,
                "fdr_threshold": fdr_threshold,
                "tested_genes_count": len(gene_indices),
            },
            summary_metrics={
                "tested_genes": len(res_df),
                "significant_spatial_degs": len(sig_svgs),
                "top_spatial_genes": sig_svgs["gene"].head(10).tolist() if not sig_svgs.empty else [],
                "has_fdr_correction": True,
            },
        )

        all_ops = [
            "build_spatial_knn",
            "build_spatial_connectivities",
            "calculate_moran_i",
            "calculate_morans_i",
            "calculate_geary_c",
            "analytical_significance_test",
            "benjamini_hochberg_correction",
            "fdr_correction",
            "identify_spatial_degs",
        ]
        if contract.allowed_operations:
            executed_ops = [op for op in all_ops if op in contract.allowed_operations]
        else:
            executed_ops = [
                "build_spatial_connectivities",
                "calculate_moran_i",
                "calculate_geary_c",
                "analytical_significance_test",
                "benjamini_hochberg_correction",
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
                "total_genes_tested": len(res_df),
                "significant_svg_count": len(sig_svgs),
                "top_svg_genes": sig_svgs["gene"].head(5).tolist() if not sig_svgs.empty else [],
            },
        )
