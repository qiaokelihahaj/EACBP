"""
Genetic Perturbation Capability: In silico CRISPR knockout and overexpression simulation.
Implements Gene Regulatory Network (GRN) propagation:
    \\Delta x = (I - \\alpha A)^{-1} v
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from scipy import stats

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
)
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def construct_grn_adjacency_from_data(
    X: np.ndarray,
    threshold: float = 0.05,
    max_degree: Optional[int] = None,
) -> np.ndarray:
    """
    Constructs a row-normalized Gene Regulatory Network (GRN) adjacency matrix A from gene expression.
    
    Parameters
    ----------
    X : np.ndarray, shape (N_cells, N_genes)
        Normalized expression matrix.
    threshold : float
        Absolute correlation threshold below which edge weights are set to 0.
    max_degree : int, optional
        Maximum number of top outgoing edges to retain per gene.
        
    Returns
    -------
    A : np.ndarray, shape (N_genes, N_genes)
        Row-normalized adjacency matrix with zero diagonal and spectral radius bounded.
    """
    n_cells, n_genes = X.shape
    if n_genes <= 1:
        return np.zeros((n_genes, n_genes), dtype=np.float32)

    # Compute correlation matrix safely
    std = np.std(X, axis=0)
    valid_genes = std > 1e-8
    
    # Gene-gene Pearson correlation
    corr = np.zeros((n_genes, n_genes), dtype=np.float32)
    if np.any(valid_genes):
        valid_indices = np.where(valid_genes)[0]
        sub_X = X[:, valid_indices]
        # Standardize sub_X
        sub_centered = sub_X - np.mean(sub_X, axis=0)
        sub_norm = np.sqrt(np.sum(sub_centered ** 2, axis=0)) + 1e-8
        sub_std = sub_centered / sub_norm
        sub_corr = np.dot(sub_std.T, sub_std)
        
        # Place into full matrix
        for i_idx, g_i in enumerate(valid_indices):
            for j_idx, g_j in enumerate(valid_indices):
                corr[g_i, g_j] = sub_corr[i_idx, j_idx]

    # Zero diagonal (no self-loops in propagation matrix)
    np.fill_diagonal(corr, 0.0)

    # Apply soft/hard threshold
    adj = np.where(np.abs(corr) >= threshold, corr, 0.0)

    # Optionally prune to top k degrees per gene
    if max_degree is not None and max_degree < n_genes:
        for i in range(n_genes):
            row = adj[i]
            top_k_indices = np.argsort(np.abs(row))[-max_degree:]
            mask = np.zeros(n_genes, dtype=bool)
            mask[top_k_indices] = True
            adj[i, ~mask] = 0.0

    # Row-normalize adjacency matrix: sum(|A_ij|) <= 1
    row_sums = np.sum(np.abs(adj), axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    A = adj / row_sums

    return A.astype(np.float32)


def compute_grn_propagator(
    A: np.ndarray,
    alpha: float = 0.25,
) -> np.ndarray:
    """
    Computes the linear propagation operator (I - \\alpha A)^{-1}.
    
    Parameters
    ----------
    A : np.ndarray, shape (G, G)
        Row-normalized GRN adjacency matrix.
    alpha : float
        Network attenuation parameter, typically in [0.05, 0.40].
        
    Returns
    -------
    M : np.ndarray, shape (G, G)
        Inverted propagation matrix.
    """
    n_genes = A.shape[0]
    # Bound alpha to ensure stability
    alpha_clamped = max(0.0, min(0.9, float(alpha)))
    I = np.eye(n_genes, dtype=np.float32)
    
    if alpha_clamped == 0.0:
        return I

    # Invert (I - alpha * A)
    system_mat = I - alpha_clamped * A
    try:
        M = np.linalg.inv(system_mat)
    except np.linalg.LinAlgError:
        # Fallback to Neumann series expansion if singular
        M = I + alpha_clamped * A + (alpha_clamped ** 2) * np.dot(A, A)

    return M.astype(np.float32)


class GeneticPerturbationCapability(BaseCapability):
    """
    Simulates in silico genetic perturbations (CRISPR KO, Knockdown, Overexpression)
    using Gene Regulatory Network (GRN) propagation dynamics.
    
    Mathematical Formulation:
        \\Delta x_i = (I - \\alpha A)^{-1} v_i
        X_{perturbed} = max(0, X + \\Delta X)
    """

    def __init__(self, implementation_id: str = "in_silico_crispr_ko_v1"):
        super().__init__(
            capability_name="genetic_perturbation_simulation",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["scRNA", "spatial"],
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
            suitable_for=["crispr_knockout", "overexpression", "grn_propagation"],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """
        Executes in silico genetic perturbation simulation within TaskContract bounds.
        """
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        X = data.X.copy()
        n_cells, n_genes = X.shape

        # Resolve gene names
        if "gene_name" in data.var.columns:
            gene_names = np.array(data.var["gene_name"].values, dtype=str)
        elif data.var.index is not None and len(data.var.index) == n_genes:
            gene_names = np.array(data.var.index.values, dtype=str)
        else:
            gene_names = np.array([f"Gene_{i}" for i in range(n_genes)], dtype=str)

        gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Contract Parameters
        params = contract.parameters
        target_gene = params.get("target_gene")
        target_genes = params.get("target_genes", [target_gene] if target_gene else [])
        perturb_type = str(params.get("perturbation_type", "knockout")).lower()
        efficiency = float(params.get("efficiency", params.get("knockout_efficiency", 0.95)))
        overexpression_factor = float(params.get("overexpression_factor", 3.0))
        alpha = float(params.get("network_attenuation", 0.25))
        custom_adjacency = params.get("grn_adjacency", None)

        if not target_genes or target_genes == [None]:
            # Default to first marker or top expressed gene if not specified
            target_genes = [gene_names[0]]

        # Validate target genes exist
        target_indices = []
        for g in target_genes:
            if g not in gene_name_to_idx:
                raise KeyError(f"Target gene '{g}' not found in dataset variables. Available genes: {list(gene_names[:10])}...")
            target_indices.append(gene_name_to_idx[g])

        # Step 1: Construct or load GRN adjacency matrix
        if custom_adjacency is not None:
            A = np.asarray(custom_adjacency, dtype=np.float32)
            if A.shape != (n_genes, n_genes):
                raise ValueError(f"Custom GRN adjacency matrix shape {A.shape} does not match gene count ({n_genes}, {n_genes}).")
        else:
            A = construct_grn_adjacency_from_data(X, threshold=0.05)

        # Step 2: Build initial perturbation matrix V (N_cells x N_genes)
        V = np.zeros((n_cells, n_genes), dtype=np.float32)
        for t_idx in target_indices:
            orig_expr = X[:, t_idx]
            if perturb_type in ("knockout", "ko"):
                # Reduce expression by efficiency delta (default 95%-100%)
                eff = min(1.0, max(0.0, efficiency))
                V[:, t_idx] = -eff * orig_expr
            elif perturb_type in ("knockdown", "kd"):
                eff = min(1.0, max(0.0, efficiency if "efficiency" in params else 0.50))
                V[:, t_idx] = -eff * orig_expr
            elif perturb_type in ("overexpression", "oe"):
                gene_std = float(np.std(orig_expr))
                added_expr = overexpression_factor * (gene_std if gene_std > 1e-4 else 1.0)
                V[:, t_idx] = added_expr
            else:
                raise ValueError(f"Unsupported perturbation_type: '{perturb_type}'. Must be 'knockout', 'knockdown', or 'overexpression'.")

        # Step 3: Propagate perturbation through GRN: Delta X = V * (I - alpha * A)^{-1}
        M = compute_grn_propagator(A, alpha=alpha)
        # Delta X = V * M (where M_ij represents downstream effect of gene i on gene j)
        delta_X = np.dot(V, M)

        # Step 4: Compute simulated post-perturbation expression matrix
        X_perturbed = np.maximum(0.0, X + delta_X).astype(np.float32)

        # Step 5: Compute state shift vectors and top perturbed downstream genes
        mean_baseline = np.mean(X, axis=0)
        mean_perturbed = np.mean(X_perturbed, axis=0)
        mean_shift = mean_perturbed - mean_baseline
        abs_shift = np.abs(mean_shift)

        # Differential expression table for perturbation
        results = []
        for i, g_name in enumerate(gene_names):
            is_target = (i in target_indices)
            b_val = float(mean_baseline[i])
            p_val = float(mean_perturbed[i])
            shift_val = float(mean_shift[i])
            rel_pct = float((shift_val / (b_val + 1e-4)) * 100.0)
            
            # Statistical test across cells
            if np.std(X[:, i]) > 1e-6 or np.std(X_perturbed[:, i]) > 1e-6:
                ttest_res = stats.ttest_rel(X_perturbed[:, i], X[:, i])
                p_value = float(ttest_res.pvalue) if not np.isnan(ttest_res.pvalue) else 1.0
            else:
                p_value = 1.0

            results.append({
                "gene": str(g_name),
                "baseline_mean": b_val,
                "perturbed_mean": p_val,
                "expression_shift": shift_val,
                "relative_shift_pct": rel_pct,
                "is_target_gene": bool(is_target),
                "p_value": p_value,
            })

        perturb_df = pd.DataFrame(results).sort_values("expression_shift", key=abs, ascending=False).reset_index(drop=True)

        # Step 6: Compute Latent State Shift & Reversion Rate
        reversion_rate = 0.0
        obs = data.obs
        cond_col = "condition" if "condition" in obs.columns else None

        # Approximate PCA projection for latent state shift
        if "X_pca" in data.obsm:
            orig_pca = data.obsm["X_pca"]
        else:
            # Simple PCA projection
            centered = X - np.mean(X, axis=0)
            u, s, vt = np.linalg.svd(centered, full_matrices=False)
            orig_pca = u[:, :min(10, n_genes)] * s[:min(10, n_genes)]

        # Project delta_X to latent space: delta_Z = delta_X * V_pca
        if n_genes > 1:
            centered_base = X - np.mean(X, axis=0)
            _, _, vt = np.linalg.svd(centered_base, full_matrices=False)
            pca_components = vt[:min(orig_pca.shape[1], vt.shape[0])].T
            delta_pca = np.dot(delta_X, pca_components)
            perturbed_pca = orig_pca[:, :pca_components.shape[1]] + delta_pca
        else:
            perturbed_pca = orig_pca.copy()

        # If disease and control conditions exist, calculate reversion rate
        if cond_col and len(obs[cond_col].unique()) >= 2:
            conditions = obs[cond_col].unique()
            cond_disease = "AD" if "AD" in conditions else conditions[0]
            cond_ctrl = "control" if "control" in conditions else conditions[1]

            mask_dis = (obs[cond_col] == cond_disease).values
            mask_ctrl = (obs[cond_col] == cond_ctrl).values

            if np.sum(mask_dis) > 0 and np.sum(mask_ctrl) > 0:
                center_dis_orig = np.mean(orig_pca[mask_dis], axis=0)
                center_ctrl = np.mean(orig_pca[mask_ctrl], axis=0)
                center_dis_perturbed = np.mean(perturbed_pca[mask_dis], axis=0)

                dist_baseline = float(np.linalg.norm(center_dis_orig - center_ctrl))
                dist_perturbed = float(np.linalg.norm(center_dis_perturbed - center_ctrl))

                if dist_baseline > 1e-6:
                    reversion_rate = float((dist_baseline - dist_perturbed) / dist_baseline)
                    reversion_rate = max(-1.0, min(1.0, reversion_rate))

        # Build output SCData
        res_data = data.copy()
        res_data.X = X_perturbed
        res_data.obsm["X_pca"] = perturbed_pca
        res_data.uns["perturbation"] = {
            "target_genes": target_genes,
            "perturbation_type": perturb_type,
            "efficiency": efficiency,
            "network_attenuation": alpha,
            "reversion_rate": reversion_rate,
            "mean_shift": mean_shift,
        }

        # Format URIs
        uri_obj = ArtifactURI.parse(in_uri)
        target_suffix = "_".join(target_genes[:2])
        if contract.expected_outputs and len(contract.expected_outputs) >= 1:
            out_adata_uri = contract.expected_outputs[0]
            out_table_uri = contract.expected_outputs[1] if len(contract.expected_outputs) >= 2 else f"table://{uri_obj.study_id}/perturbation_summary_{target_suffix}/v1"
        else:
            out_adata_uri = f"adata://{uri_obj.study_id}/perturbation_ko_{target_suffix}/v1" if perturb_type in ("knockout", "ko") else f"adata://{uri_obj.study_id}/perturbation_{perturb_type}_{target_suffix}/v1"
            out_table_uri = f"table://{uri_obj.study_id}/perturbation_summary_{target_suffix}/v1"

        # Increment version if already registered (Invariant 2)
        adata_uri_obj = ArtifactURI.parse(out_adata_uri)
        while registry.exists(adata_uri_obj.to_string()):
            adata_uri_obj = adata_uri_obj.next_version()
        out_adata_uri = adata_uri_obj.to_string()

        table_uri_obj = ArtifactURI.parse(out_table_uri)
        while registry.exists(table_uri_obj.to_string()):
            table_uri_obj = table_uri_obj.next_version()
        out_table_uri = table_uri_obj.to_string()

        # Register Artifacts
        registry.register(
            uri_str=out_adata_uri,
            payload=res_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="simulate_genetic_perturbation",
            parent_uris=[in_uri],
            parameters={
                "target_genes": target_genes,
                "perturbation_type": perturb_type,
                "efficiency": efficiency,
                "network_attenuation": alpha,
            },
            summary_metrics={
                "target_genes": target_genes,
                "reversion_rate": reversion_rate,
                "perturbed_cells": n_cells,
                "mean_absolute_shift": float(np.mean(abs_shift)),
            }
        )

        registry.register(
            uri_str=out_table_uri,
            payload=perturb_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="summarize_perturbation_shifts",
            parent_uris=[in_uri],
            summary_metrics={
                "target_genes": target_genes,
                "top_downstream_genes": perturb_df[~perturb_df["is_target_gene"]]["gene"].head(5).tolist(),
            }
        )

        top_downstream = perturb_df[~perturb_df["is_target_gene"]]["gene"].head(5).tolist()

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_adata_uri, out_table_uri],
            executed_operations=[
                "construct_grn_adjacency",
                "simulate_genetic_perturbation",
                "propagate_network_shift",
                "compute_state_reversion",
            ],
            metrics={
                "target_genes": target_genes,
                "perturbation_type": perturb_type,
                "reversion_rate": reversion_rate,
                "mean_absolute_shift": float(np.mean(abs_shift)),
                "top_perturbed_genes": perturb_df["gene"].head(10).tolist(),
                "top_downstream_genes": top_downstream,
            }
        )


def generate_genetic_perturbation_evidence(
    contract: TaskContract,
    result: TaskResult,
    target_gene: str,
    reversion_rate: float,
    perturbation_type: str = "knockout",
) -> EvidenceNode:
    """
    Generates a calibrated EvidenceType.PERTURBATION node from genetic simulation results.
    Strictly caps in silico causal confidence score at <= 0.50.
    """
    out_uris = result.output_artifacts
    task_id = contract.task_id
    
    # Quantitative score calibrated and capped at 0.50
    normalized_score = max(0.10, min(0.50, float(abs(reversion_rate))))

    strength = (
        EvidenceStrength.STRONG if abs(reversion_rate) >= 0.40
        else (EvidenceStrength.MODERATE if abs(reversion_rate) >= 0.15 else EvidenceStrength.WEAK)
    )

    action_label = "knockout" if perturbation_type in ("knockout", "ko") else ("overexpression" if perturbation_type in ("overexpression", "oe") else "perturbation")
    direction_label = "attenuates" if reversion_rate > 0 else "exacerbates"

    summary = (
        f"In silico {action_label} simulation of {target_gene} predicts {abs(reversion_rate)*100:.1f}% "
        f"state reversion ({direction_label} disease phenotype) via gene regulatory network propagation."
    )

    return EvidenceNode(
        evidence_id=f"E_perturb_{target_gene}_{task_id}",
        type=EvidenceType.PERTURBATION,
        polarity=EvidencePolarity.SUPPORTING if reversion_rate >= 0 else EvidencePolarity.CONTRADICTING,
        strength=strength,
        score=normalized_score,
        summary=summary,
        source_task_id=task_id,
        source_artifact_uris=out_uris,
        metrics={
            "target_gene": target_gene,
            "perturbation_type": perturbation_type,
            "reversion_rate": reversion_rate,
            "in_silico_confidence_cap": 0.50,
            "model": "grn_linear_propagation",
        },
        biological_context={
            "target_gene": target_gene,
            "causal_status": "in_silico_perturbed",
        },
    )
