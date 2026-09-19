"""
Genetic Perturbation Capability: In silico CRISPR knockout and overexpression simulation.
Implements Gene Regulatory Network (GRN) propagation:
    \\Delta x = (I - \\alpha A)^{-1} v
"""

from typing import Dict, Any, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse

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


MAX_DENSE_GRN_GENES = 512
MAX_PERTURBATION_ELEMENTS = 20_000_000
MAX_GRN_CORRELATION_ELEMENTS = 20_000_000


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
    if sparse.issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)
    n_cells, n_genes = X.shape
    if int(n_cells) * int(n_genes) > MAX_GRN_CORRELATION_ELEMENTS:
        raise ValueError(
            f"GRN correlation is limited to {MAX_GRN_CORRELATION_ELEMENTS:,} input elements; provide a precomputed sparse adjacency"
        )
    if n_genes <= 1:
        return np.zeros((n_genes, n_genes), dtype=np.float32)

    # A bounded dense path keeps the simple implementation usable for small
    # demos.  Larger matrices use blockwise top-k correlations and return CSR
    # so an accidental GxG allocation cannot occur.
    std = np.std(X, axis=0)
    valid_genes = std > 1e-8
    
    valid_indices = np.where(valid_genes)[0]
    if len(valid_indices) == 0:
        return np.zeros((n_genes, n_genes), dtype=np.float32) if n_genes <= MAX_DENSE_GRN_GENES else sparse.csr_matrix((n_genes, n_genes), dtype=np.float32)

    centered = X[:, valid_indices] - np.mean(X[:, valid_indices], axis=0, keepdims=True)
    normalized = centered / (np.sqrt(np.sum(centered ** 2, axis=0, keepdims=True)) + 1e-8)
    degree = min(max(1, int(max_degree or 32)), max(1, n_genes - 1))

    if n_genes <= MAX_DENSE_GRN_GENES:
        corr = np.zeros((n_genes, n_genes), dtype=np.float32)
        corr_sub = normalized.T @ normalized
        corr[np.ix_(valid_indices, valid_indices)] = corr_sub
        np.fill_diagonal(corr, 0.0)
        adj = np.where(np.abs(corr) >= threshold, corr, 0.0)
        if max_degree is not None and degree < n_genes:
            for i in range(n_genes):
                keep = np.argpartition(np.abs(adj[i]), -degree)[-degree:]
                mask = np.zeros(n_genes, dtype=bool)
                mask[keep] = True
                adj[i, ~mask] = 0.0
        row_sums = np.sum(np.abs(adj), axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0
        return (adj / row_sums).astype(np.float32)

    # Blockwise top-k correlation edges.  The largest temporary is
    # (block_size x number_of_genes), never a dense GxG matrix.
    rows, cols, vals = [], [], []
    block_size = 64
    for start in range(0, len(valid_indices), block_size):
        block_gene_idx = valid_indices[start : start + block_size]
        scores = normalized[:, start : start + len(block_gene_idx)].T @ normalized
        for local_i, gene_i in enumerate(block_gene_idx):
            scores[local_i, gene_i] = 0.0
            candidate = np.flatnonzero(np.abs(scores[local_i]) >= threshold)
            if len(candidate) > degree:
                candidate = candidate[np.argpartition(np.abs(scores[local_i, candidate]), -degree)[-degree:]]
            for gene_j in candidate:
                rows.append(int(gene_i)); cols.append(int(gene_j)); vals.append(float(scores[local_i, gene_j]))
    adj = sparse.coo_matrix((vals, (rows, cols)), shape=(n_genes, n_genes), dtype=np.float32).tocsr()
    row_sums = np.asarray(np.abs(adj).sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    return sparse.diags(1.0 / row_sums).dot(adj).tocsr().astype(np.float32)


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
    if sparse.issparse(A):
        A_work = A.tocsr().astype(np.float32)
        n_genes = A_work.shape[0]
        sparse_output = True
    else:
        A_work = np.asarray(A, dtype=np.float32)
        if A_work.ndim != 2 or A_work.shape[0] != A_work.shape[1]:
            raise ValueError("GRN adjacency must be a square matrix")
        n_genes = A_work.shape[0]
        if n_genes > MAX_DENSE_GRN_GENES:
            raise ValueError(
                f"Dense GRN propagation is limited to {MAX_DENSE_GRN_GENES} genes; provide a CSR adjacency for larger networks"
            )
        sparse_output = False
    alpha_clamped = max(0.0, min(0.9, float(alpha)))
    if sparse_output:
        identity = sparse.eye(n_genes, format="csr", dtype=np.float32)
    else:
        identity = np.eye(n_genes, dtype=np.float32)
    if alpha_clamped == 0.0:
        return identity

    # Neumann truncation avoids an unbounded dense inverse.  Row L1-normalised
    # A and alpha<1 make the omitted tail geometrically bounded.
    terms = 12
    result = identity.copy()
    power = identity.copy()
    for order in range(1, terms + 1):
        power = power.dot(A_work) if sparse_output else power @ A_work
        result = result + (alpha_clamped ** order) * power
    return result.astype(np.float32)


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
        n_cells, n_genes = data.X.shape
        if int(n_cells) * int(n_genes) > MAX_PERTURBATION_ELEMENTS:
            raise ValueError(
                f"Genetic perturbation dense cell-state simulation is limited to {MAX_PERTURBATION_ELEMENTS:,} matrix elements; "
                "subset cells/genes or provide a sparse/local model"
            )
        X = data.X.toarray() if hasattr(data.X, "toarray") else np.asarray(data.X, dtype=np.float32).copy()

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
            raise ValueError(
                "Genetic perturbation requires an explicit target_gene or target_genes; "
                "a default marker would make the perturbation hypothesis implicit"
            )

        # Validate target genes exist
        target_indices = []
        for g in target_genes:
            if g not in gene_name_to_idx:
                raise KeyError(f"Target gene '{g}' not found in dataset variables. Available genes: {list(gene_names[:10])}...")
            target_indices.append(gene_name_to_idx[g])

        # Step 1: Construct or load GRN adjacency matrix
        network_source = None
        network_verified = False
        if custom_adjacency is not None:
            if sparse.issparse(custom_adjacency):
                A = custom_adjacency.tocsr().astype(np.float32)
            else:
                if n_genes > MAX_DENSE_GRN_GENES:
                    raise ValueError(
                        f"Dense custom GRN adjacency is limited to {MAX_DENSE_GRN_GENES} genes; pass scipy.sparse.csr_matrix for a larger local graph"
                    )
                A = np.asarray(custom_adjacency, dtype=np.float32)
            if A.shape != (n_genes, n_genes):
                raise ValueError(f"Custom GRN adjacency matrix shape {A.shape} does not match gene count ({n_genes}, {n_genes}).")
            if not np.isfinite(A.data if sparse.issparse(A) else A).all():
                raise ValueError("Custom GRN adjacency contains non-finite values")
            network_source = str(params.get("grn_source", "user_provided_unverified"))
            network_verified = bool(params.get("grn_verified", False))
        else:
            A = construct_grn_adjacency_from_data(X, threshold=0.05)
            # Correlations are useful for a counterfactual demo, but they are
            # not a validated regulatory network.  Record that distinction in
            # every output so downstream reports cannot present the result as
            # an experimentally established GRN effect.
            network_source = "expression_correlation_inferred"
            network_verified = False

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
        delta_X = V.dot(M)
        if sparse.issparse(delta_X):
            delta_X = delta_X.toarray()
        delta_X = np.asarray(delta_X, dtype=np.float32)

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
        reversion_rate = None
        state_reversion_available = False
        obs = data.obs
        cond_col = "condition" if "condition" in obs.columns else None

        # Both states must use the same fitted basis. Existing Harmony/PCA
        # coordinates cannot be added to a delta projected onto a new basis.
        centered_base = X - np.mean(X, axis=0)
        _, _, vt = np.linalg.svd(centered_base, full_matrices=False)
        pca_components = vt[:min(10, vt.shape[0])].T
        orig_pca = centered_base @ pca_components
        perturbed_pca = orig_pca + delta_X @ pca_components

        # If disease and control conditions exist, calculate reversion rate
        reversion_conditions = None
        if cond_col and len(obs[cond_col].dropna().unique()) >= 2:
            conditions = list(obs[cond_col].dropna().unique())
            cond_disease = params.get("condition_disease")
            cond_ctrl = params.get("condition_control")
            if cond_disease is None or cond_ctrl is None:
                # Infer only common, unambiguous labels.  Arbitrary A/B labels
                # do not establish which group is a disease baseline.
                normalized = {str(value).strip().lower(): value for value in conditions}
                disease_keys = ("ad", "disease", "case", "treated", "cko", "mutant")
                control_keys = ("control", "ctrl", "healthy", "wt", "wildtype", "con")
                cond_disease = next((normalized[key] for key in disease_keys if key in normalized), None)
                cond_ctrl = next((normalized[key] for key in control_keys if key in normalized), None)
            if cond_disease is not None and cond_ctrl is not None and cond_disease != cond_ctrl:
                reversion_conditions = {"condition_disease": str(cond_disease), "condition_control": str(cond_ctrl)}

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
                        state_reversion_available = True

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
            "reversion_conditions": reversion_conditions,
            "mean_shift": mean_shift,
            "network_source": network_source,
            "network_verified": network_verified,
            "network_semantics": "regulatory_network" if network_verified else "expression_association_counterfactual",
            "state_reversion_available": state_reversion_available,
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
                "network_source": network_source,
                "network_verified": network_verified,
                "state_reversion_available": state_reversion_available,
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
                "network_source": network_source,
                "network_verified": network_verified,
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
                "state_reversion_available": state_reversion_available,
                "mean_absolute_shift": float(np.mean(abs_shift)),
                "top_perturbed_genes": perturb_df["gene"].head(10).tolist(),
                "top_downstream_genes": top_downstream,
                "network_source": network_source,
                "network_verified": network_verified,
            }
        )


def generate_genetic_perturbation_evidence(
    contract: TaskContract,
    result: TaskResult,
    target_gene: str,
    reversion_rate: Optional[float],
    perturbation_type: str = "knockout",
) -> EvidenceNode:
    """
    Generates a calibrated EvidenceType.PERTURBATION node from genetic simulation results.
    Strictly caps in silico causal confidence score at <= 0.50.
    """
    out_uris = result.output_artifacts
    task_id = contract.task_id
    
    state_reversion_available = result.metrics.get(
        "state_reversion_available", reversion_rate is not None
    )
    network_source = result.metrics.get("network_source", "unspecified")
    network_verified = bool(result.metrics.get("network_verified", False))

    # A missing state annotation is a normal, explicitly represented outcome;
    # it must not be turned into a zero/positive reversion claim by a default.
    if reversion_rate is None or not state_reversion_available:
        return EvidenceNode(
            evidence_id=f"E_perturb_{target_gene}_{task_id}",
            type=EvidenceType.PERTURBATION,
            polarity=EvidencePolarity.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
            score=0.0,
            summary=(
                f"In silico {perturbation_type} simulation of {target_gene} produced a "
                "counterfactual expression matrix, but state reversion was not estimable "
                "because no validated disease/control state labels were available."
            ),
            source_task_id=task_id,
            source_artifact_uris=out_uris,
            metrics={
                "target_gene": target_gene,
                "perturbation_type": perturbation_type,
                "state_reversion_available": False,
                "network_source": network_source,
                "network_verified": network_verified,
                "in_silico_confidence_cap": 0.50,
            },
            biological_context={
                "target_gene": target_gene,
                "causal_status": "in_silico_perturbed",
                "state_claim": "not_estimable",
            },
            is_simulated=True,
            source_verified=False,
        )

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
            "network_source": network_source,
            "network_verified": network_verified,
            "state_reversion_available": True,
        },
        biological_context={
            "target_gene": target_gene,
            "causal_status": "in_silico_perturbed",
            "network_semantics": "validated_regulatory_network" if network_verified else "expression_association_counterfactual",
        },
        is_simulated=True,
        source_verified=False,
    )
