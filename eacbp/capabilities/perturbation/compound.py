"""
Compound Perturbation Capability: In silico small molecule / drug response simulation
and counterfactual state transition modeling using CMAP-style signature discordance scores
and transition probability matrices.

Mathematical Formulations:
    Discordance Score (Reversal) = -\\cos(s_{disease}, s_{drug}) = -\\frac{s_{disease} \\cdot s_{drug}}{\\|s_{disease}\\|_2 \\|s_{drug}\\|_2}
    P(State_v | Cell_i) = \\frac{\\exp(-d(x_{i, drug}, c_v)^2 / (2\\sigma^2))}{\\sum_w \\exp(-d(x_{i, drug}, c_w)^2 / (2\\sigma^2))}
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


# Built-in reference compound signatures (log2FC profiles for known mechanisms)
REFERENCE_COMPOUND_DATABASE: Dict[str, Dict[str, float]] = {
    "Bexarotene": {
        "Apoe": -1.8,
        "Trem2": -1.4,
        "Clec7a": -2.1,
        "Itgax": -1.9,
        "P2ry12": 1.6,
        "Cx3cr1": 1.5,
        "Tmem119": 1.4,
        "C1qa": -0.8,
        "C1qb": -0.7,
        "Cst7": -1.6,
    },
    "GW3965": {
        "Apoe": -1.5,
        "Trem2": -1.2,
        "Clec7a": -1.7,
        "P2ry12": 1.3,
        "Cx3cr1": 1.2,
        "Tmem119": 1.1,
        "Lrp1": 1.4,
        "Abca1": 2.0,
    },
    "Anti_Inflammatory_Small_Molecule": {
        "Apoe": -1.6,
        "Trem2": -1.3,
        "Clec7a": -1.9,
        "Itgax": -1.5,
        "P2ry12": 1.5,
        "Cx3cr1": 1.4,
        "Tmem119": 1.3,
        "Il1b": -2.2,
        "Tnf": -2.0,
    },
    "Mock_Exacerbator": {
        "Apoe": 2.0,
        "Trem2": 1.8,
        "Clec7a": 2.2,
        "Itgax": 1.9,
        "P2ry12": -1.8,
        "Cx3cr1": -1.5,
        "Tmem119": -1.4,
    },
}


def compute_cmap_cosine_discordance(
    disease_sig: np.ndarray,
    drug_sig: np.ndarray,
    n_permutations: int = 500,
    random_seed: int = 42,
) -> Tuple[float, float, float]:
    """
    Computes CMAP-style cosine discordance score and empirical permutation p-value.
    
    Discordance Score = -cos(s_disease, s_drug) = - (s_disease . s_drug) / (||s_disease|| * ||s_drug||)
    
    Returns:
        discordance_score: float in [-1.0, 1.0] (positive values indicate therapeutic reversal)
        cosine_similarity: float in [-1.0, 1.0]
        p_value: empirical significance of reversal score
    """
    norm_dis = np.linalg.norm(disease_sig)
    norm_drug = np.linalg.norm(drug_sig)

    if norm_dis < 1e-8 or norm_drug < 1e-8:
        return 0.0, 0.0, 1.0

    cosine_sim = float(np.dot(disease_sig, drug_sig) / (norm_dis * norm_drug))
    discordance_score = float(-cosine_sim)

    # Permutation test
    np.random.seed(random_seed)
    perm_scores = []
    
    for _ in range(n_permutations):
        shuffled_drug = np.random.permutation(drug_sig)
        p_sim = np.dot(disease_sig, shuffled_drug) / (norm_dis * norm_drug)
        perm_scores.append(-p_sim)

    perm_scores = np.array(perm_scores)
    # One-tailed p-value for reversal: how often null score >= observed discordance
    p_val = float((np.sum(perm_scores >= discordance_score) + 1.0) / (n_permutations + 1.0))

    return discordance_score, cosine_sim, p_val


class CompoundPerturbationCapability(BaseCapability):
    """
    Simulates small molecule / compound response and counterfactual cell state transitions.
    Evaluates transcriptomic reversal against disease signatures and calculates
    cell state transition probability matrices.
    """

    def __init__(self, implementation_id: str = "in_silico_compound_response_v1"):
        super().__init__(
            capability_name="compound_perturbation_simulation",
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["scRNA", "spatial"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
            output_types=[ArtifactType.TABLE, ArtifactType.ANNDATA],
            suitable_for=["drug_response", "counterfactual_transition", "cmap_matching"],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        # Handle input either as SCData or Table
        if isinstance(payload, SCData) or (isinstance(payload, dict) and "X" in payload):
            data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
            X = data.X.copy()
            obs = data.obs.copy()
            n_cells, n_genes = X.shape
            
            if "gene_name" in data.var.columns:
                gene_names = np.array(data.var["gene_name"].values, dtype=str)
            elif data.var.index is not None and len(data.var.index) == n_genes:
                gene_names = np.array(data.var.index.values, dtype=str)
            else:
                gene_names = np.array([f"Gene_{i}" for i in range(n_genes)], dtype=str)
        else:
            # If input is a DEG table
            deg_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            gene_names = np.array(deg_df["gene"].values, dtype=str)
            n_genes = len(gene_names)
            n_cells = 100
            X = np.ones((n_cells, n_genes), dtype=np.float32)
            obs = pd.DataFrame({"cell_id": [f"c_{i}" for i in range(n_cells)], "condition": ["AD"] * 50 + ["control"] * 50})
            data = SCData(X=X, obs=obs, var=pd.DataFrame({"gene_name": gene_names}))

        gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Contract Parameters
        params = contract.parameters
        compound_name = str(params.get("compound_name", "Anti_Inflammatory_Small_Molecule"))
        dosage = float(params.get("dosage", params.get("scale_factor", 1.0)))
        n_perms = int(params.get("n_permutations", 500))
        custom_drug_sig = params.get("drug_signature", None)
        custom_disease_sig = params.get("disease_signature", None)

        # 1. Determine Disease Signature s_disease (1 x N_genes)
        disease_sig = np.zeros(n_genes, dtype=np.float32)
        if custom_disease_sig is not None:
            if isinstance(custom_disease_sig, dict):
                for g, val in custom_disease_sig.items():
                    if g in gene_name_to_idx:
                        disease_sig[gene_name_to_idx[g]] = float(val)
            elif len(custom_disease_sig) == n_genes:
                disease_sig = np.asarray(custom_disease_sig, dtype=np.float32)
        else:
            # Compute from condition in obs: AD vs Control
            cond_col = "condition" if "condition" in obs.columns else None
            if cond_col and len(obs[cond_col].unique()) >= 2:
                conditions = obs[cond_col].unique()
                cond_ad = "AD" if "AD" in conditions else conditions[0]
                cond_ctrl = "control" if "control" in conditions else conditions[1]

                mean_ad = np.mean(X[obs[cond_col] == cond_ad], axis=0)
                mean_ctrl = np.mean(X[obs[cond_col] == cond_ctrl], axis=0)
                disease_sig = np.log2((mean_ad + 1e-3) / (mean_ctrl + 1e-3)).astype(np.float32)
            else:
                # Fallback: variance-weighted difference
                disease_sig = (X[0] - np.mean(X, axis=0)).astype(np.float32)

        # 2. Determine Drug Perturbation Signature s_drug (1 x N_genes)
        drug_sig = np.zeros(n_genes, dtype=np.float32)
        if custom_drug_sig is not None:
            if isinstance(custom_drug_sig, dict):
                for g, val in custom_drug_sig.items():
                    if g in gene_name_to_idx:
                        drug_sig[gene_name_to_idx[g]] = float(val)
            elif len(custom_drug_sig) == n_genes:
                drug_sig = np.asarray(custom_drug_sig, dtype=np.float32)
        elif compound_name in REFERENCE_COMPOUND_DATABASE:
            ref_dict = REFERENCE_COMPOUND_DATABASE[compound_name]
            for g, val in ref_dict.items():
                if g in gene_name_to_idx:
                    drug_sig[gene_name_to_idx[g]] = float(val)
        else:
            # Synthetic reversal signature matching inverse of top disease genes
            top_dis_idx = np.argsort(np.abs(disease_sig))[-min(10, n_genes):]
            for idx in top_dis_idx:
                drug_sig[idx] = -1.5 * np.sign(disease_sig[idx])

        # 3. Compute CMAP Cosine Discordance Score
        reversal_score, cosine_sim, p_val = compute_cmap_cosine_discordance(
            disease_sig=disease_sig,
            drug_sig=drug_sig,
            n_permutations=n_perms,
            random_seed=42,
        )

        # 4. Simulate Counterfactual Treatment & Cell State Transitions
        # X_drug = max(0, X + dosage * drug_sig)
        X_drug = np.maximum(0.0, X + dosage * drug_sig).astype(np.float32)

        # Identify Cell States (e.g. microglia_state, leiden, or condition)
        state_col = (
            "microglia_state" if "microglia_state" in obs.columns
            else ("leiden" if "leiden" in obs.columns else ("condition" if "condition" in obs.columns else None))
        )

        if state_col:
            unique_states = [str(s) for s in obs[state_col].unique()]
        else:
            unique_states = ["State_0", "State_1"]
            obs["inferred_state"] = ["State_0" if i < n_cells // 2 else "State_1" for i in range(n_cells)]
            state_col = "inferred_state"

        n_states = len(unique_states)
        state_to_idx = {s: i for i, s in enumerate(unique_states)}

        # Compute baseline centroids for each state
        centroids = []
        for s in unique_states:
            mask = (obs[state_col] == s).values
            if np.sum(mask) > 0:
                centroids.append(np.mean(X[mask], axis=0))
            else:
                centroids.append(np.mean(X, axis=0))
        centroids = np.array(centroids, dtype=np.float32)  # (K, G)

        # For each drug-treated cell, compute probability of transitioning to each state
        # Distance to state centroids
        diffs = X_drug[:, np.newaxis, :] - centroids[np.newaxis, :, :]  # (N, K, G)
        sq_dists = np.sum(diffs ** 2, axis=2)  # (N, K)

        # Softmax over negative distances
        sigma_sq = float(np.mean(sq_dists)) + 1e-4
        scaled_logits = -sq_dists / (2.0 * sigma_sq)
        scaled_logits -= np.max(scaled_logits, axis=1, keepdims=True)
        exp_logits = np.exp(scaled_logits)
        trans_probs = exp_logits / np.sum(exp_logits, axis=1, keepdims=True)  # (N, K)

        predicted_state_idx = np.argmax(trans_probs, axis=1)
        predicted_states = [unique_states[idx] for idx in predicted_state_idx]

        # Compute State-to-State Transition Probability Matrix T (K x K)
        T_mat = np.zeros((n_states, n_states), dtype=np.float32)
        for u_idx, u_state in enumerate(unique_states):
            u_mask = (obs[state_col] == u_state).values
            if np.sum(u_mask) > 0:
                T_mat[u_idx, :] = np.mean(trans_probs[u_mask], axis=0)
            else:
                T_mat[u_idx, u_idx] = 1.0

        # Compute disease-to-healthy transition rate
        disease_transition_rate = 0.0
        dam_states = [s for s in unique_states if "DAM" in s or "M3" in s or "AD" in s]
        homeo_states = [s for s in unique_states if "Homeo" in s or "M1" in s or "control" in s]

        if dam_states and homeo_states:
            dam_u = state_to_idx[dam_states[0]]
            homeo_v = state_to_idx[homeo_states[0]]
            disease_transition_rate = float(T_mat[dam_u, homeo_v])
        elif n_states >= 2:
            disease_transition_rate = float(T_mat[0, 1])

        # Transition Matrix DataFrame
        trans_df = pd.DataFrame(
            T_mat,
            index=[f"from_{s}" for s in unique_states],
            columns=[f"to_{s}" for s in unique_states],
        )

        # Gene-level signature comparison table
        top_reversed_genes = []
        gene_summary_records = []
        for i, g_name in enumerate(gene_names):
            d_val = float(disease_sig[i])
            c_val = float(drug_sig[i])
            is_reversed = bool(d_val * c_val < -1e-4)
            if is_reversed:
                top_reversed_genes.append(str(g_name))

            gene_summary_records.append({
                "gene": str(g_name),
                "disease_log2fc": d_val,
                "compound_log2fc": c_val,
                "reversal_concordance": - (d_val * c_val),
                "is_reversed": is_reversed,
            })

        gene_df = pd.DataFrame(gene_summary_records).sort_values("reversal_concordance", ascending=False).reset_index(drop=True)

        # Step 6: Create and Register Output Artifacts
        uri_obj = ArtifactURI.parse(in_uri)
        clean_name = compound_name.lower().replace(" ", "_")
        if contract.expected_outputs and len(contract.expected_outputs) >= 1:
            out_table_uri = contract.expected_outputs[0]
            out_adata_uri = contract.expected_outputs[1] if len(contract.expected_outputs) >= 2 else f"adata://{uri_obj.study_id}/perturbation_drug_{clean_name}/v1"
        else:
            out_table_uri = f"table://{uri_obj.study_id}/compound_perturbation_{clean_name}/v1"
            out_adata_uri = f"adata://{uri_obj.study_id}/perturbation_drug_{clean_name}/v1"

        # Increment versions if already exists (Invariant 2)
        table_uri_obj = ArtifactURI.parse(out_table_uri)
        while registry.exists(table_uri_obj.to_string()):
            table_uri_obj = table_uri_obj.next_version()
        out_table_uri = table_uri_obj.to_string()

        adata_uri_obj = ArtifactURI.parse(out_adata_uri)
        while registry.exists(adata_uri_obj.to_string()):
            adata_uri_obj = adata_uri_obj.next_version()
        out_adata_uri = adata_uri_obj.to_string()

        # Update SCData with drug treatment
        res_data = data.copy()
        res_data.X = X_drug
        for k_idx, s in enumerate(unique_states):
            res_data.obs[f"prob_transition_{s}"] = trans_probs[:, k_idx]
        res_data.obs["predicted_transition_state"] = predicted_states
        res_data.uns["compound_perturbation"] = {
            "compound_name": compound_name,
            "dosage": dosage,
            "reversal_score": reversal_score,
            "cosine_similarity": cosine_sim,
            "p_value": p_val,
            "transition_matrix": T_mat.tolist(),
            "disease_transition_rate": disease_transition_rate,
        }

        # Register artifacts
        registry.register(
            uri_str=out_table_uri,
            payload=trans_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="simulate_compound_response",
            parent_uris=[in_uri],
            parameters={
                "compound_name": compound_name,
                "dosage": dosage,
                "n_permutations": n_perms,
            },
            summary_metrics={
                "compound_name": compound_name,
                "reversal_score": reversal_score,
                "cosine_similarity": cosine_sim,
                "p_value": p_val,
                "disease_transition_rate": disease_transition_rate,
                "top_reversed_genes": top_reversed_genes[:5],
                "therapeutic_potential": bool(reversal_score > 0.30 and p_val < 0.05),
            }
        )

        registry.register(
            uri_str=out_adata_uri,
            payload=res_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="simulate_counterfactual_cells",
            parent_uris=[in_uri],
            parameters={"compound_name": compound_name, "dosage": dosage},
            summary_metrics={
                "compound_name": compound_name,
                "reversal_score": reversal_score,
                "cells_shifted": n_cells,
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_table_uri, out_adata_uri],
            executed_operations=[
                "compute_disease_signature",
                "calculate_cmap_discordance",
                "simulate_counterfactual_transitions",
                "compute_transition_matrix",
            ],
            metrics={
                "compound_name": compound_name,
                "reversal_score": reversal_score,
                "cosine_similarity": cosine_sim,
                "p_value": p_val,
                "disease_transition_rate": disease_transition_rate,
                "top_reversed_genes": top_reversed_genes[:5],
                "therapeutic_potential": bool(reversal_score > 0.30 and p_val < 0.05),
            }
        )


def generate_compound_perturbation_evidence(
    contract: TaskContract,
    result: TaskResult,
    compound_name: str,
    reversal_score: float,
    transition_rate: float,
    p_value: float = 0.01,
) -> EvidenceNode:
    """
    Generates a calibrated EvidenceType.PERTURBATION node from small molecule response simulation.
    Strictly caps in silico causal confidence score at <= 0.50.
    """
    out_uris = result.output_artifacts
    task_id = contract.task_id

    # Quantitative score calibrated and capped at 0.50
    normalized_score = max(0.10, min(0.50, float(max(0.0, reversal_score))))

    strength = (
        EvidenceStrength.STRONG if reversal_score >= 0.50
        else (EvidenceStrength.MODERATE if reversal_score >= 0.20 else EvidenceStrength.WEAK)
    )

    polarity = EvidencePolarity.SUPPORTING if reversal_score > 0 else EvidencePolarity.CONTRADICTING
    action = "therapeutic reversal" if reversal_score > 0 else "disease exacerbation"

    summary = (
        f"In silico compound simulation with {compound_name} predicts {action} of disease signature "
        f"(CMAP discordance score: {reversal_score:.2f}, p-val: {p_value:.3f}, "
        f"counterfactual homeostatic transition rate: {transition_rate*100:.1f}%)."
    )

    return EvidenceNode(
        evidence_id=f"E_compound_{compound_name}_{task_id}",
        type=EvidenceType.PERTURBATION,
        polarity=polarity,
        strength=strength,
        score=normalized_score,
        summary=summary,
        source_task_id=task_id,
        source_artifact_uris=out_uris,
        metrics={
            "compound_name": compound_name,
            "reversal_score": reversal_score,
            "transition_rate": transition_rate,
            "p_value": p_value,
            "in_silico_confidence_cap": 0.50,
        },
        biological_context={
            "compound": compound_name,
            "causal_status": "in_silico_perturbed",
        },
    )
