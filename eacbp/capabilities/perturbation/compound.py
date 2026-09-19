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


MAX_COMPOUND_ELEMENTS = 20_000_000


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
    disease_sig = np.asarray(disease_sig, dtype=np.float32)
    drug_sig = np.asarray(drug_sig, dtype=np.float32)
    if disease_sig.ndim != 1 or drug_sig.ndim != 1 or disease_sig.shape != drug_sig.shape:
        raise ValueError("Disease and drug signatures must be one-dimensional vectors with equal length")
    if not np.isfinite(disease_sig).all() or not np.isfinite(drug_sig).all():
        raise ValueError("Disease and drug signatures must contain only finite values")
    if int(n_permutations) < 1:
        raise ValueError("n_permutations must be at least 1")

    norm_dis = np.linalg.norm(disease_sig)
    norm_drug = np.linalg.norm(drug_sig)

    if norm_dis < 1e-8 or norm_drug < 1e-8:
        return 0.0, 0.0, 1.0

    cosine_sim = float(np.dot(disease_sig, drug_sig) / (norm_dis * norm_drug))
    discordance_score = float(-cosine_sim)

    # Permutation test
    rng = np.random.default_rng(random_seed)
    perm_scores = []
    
    for _ in range(int(n_permutations)):
        shuffled_drug = rng.permutation(drug_sig)
        p_sim = np.dot(disease_sig, shuffled_drug) / (norm_dis * norm_drug)
        perm_scores.append(-p_sim)

    perm_scores = np.array(perm_scores)
    # One-tailed p-value for reversal: how often null score >= observed discordance
    p_val = float((np.sum(perm_scores >= discordance_score) + 1.0) / (int(n_permutations) + 1.0))

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
        params = contract.parameters
        simulation_mode = str(params.get("mode", "real")).strip().lower()
        if simulation_mode not in {"real", "demo", "simulation", "simulated"}:
            raise ValueError("mode must be one of: real, demo, simulation")

        # Handle input either as SCData or Table
        if isinstance(payload, SCData) or (isinstance(payload, dict) and "X" in payload):
            data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
            n_cells, n_genes = data.X.shape
            if int(n_cells) * int(n_genes) > MAX_COMPOUND_ELEMENTS:
                raise ValueError(
                    f"Compound counterfactual simulation is limited to {MAX_COMPOUND_ELEMENTS:,} expression elements; subset the input"
                )
            X = data.X.toarray() if hasattr(data.X, "toarray") else np.asarray(data.X, dtype=np.float32).copy()
            obs = data.obs.copy()
            n_cells, n_genes = X.shape
            
            if "gene_name" in data.var.columns:
                gene_names = np.array(data.var["gene_name"].values, dtype=str)
            elif data.var.index is not None and len(data.var.index) == n_genes:
                gene_names = np.array(data.var.index.values, dtype=str)
            else:
                gene_names = np.array([f"Gene_{i}" for i in range(n_genes)], dtype=str)
        else:
            # A DEG table has gene-level summaries and cannot supply the cell
            # states needed for a counterfactual.  Never manufacture a fixed
            # 100-cell all-ones matrix.  An explicit demo may provide the
            # actual expression matrix and labels through the contract.
            deg_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            if "gene" not in deg_df.columns:
                raise ValueError("Compound perturbation table input requires a 'gene' column")
            expression_matrix = params.get("expression_matrix")
            obs_input = params.get("obs")
            if expression_matrix is None or obs_input is None:
                raise ValueError(
                    "A DEG table does not contain cell-level expression or state labels; "
                    "provide an AnnData/SCData input or explicitly provide expression_matrix and obs"
                )
            if simulation_mode not in {"demo", "simulation", "simulated"}:
                raise ValueError(
                    "DEG-table counterfactuals require explicit mode='demo' or mode='simulation'"
                )
            gene_names = np.array(deg_df["gene"].values, dtype=str)
            if isinstance(expression_matrix, pd.DataFrame):
                if set(gene_names).issubset(expression_matrix.columns):
                    X = expression_matrix.loc[:, list(gene_names)].to_numpy(dtype=np.float32)
                else:
                    X = expression_matrix.to_numpy(dtype=np.float32)
            else:
                X = np.asarray(expression_matrix, dtype=np.float32)
            if X.ndim != 2 or X.shape[1] != len(gene_names):
                raise ValueError(
                    "expression_matrix must be a 2D cells x DEG-genes matrix with one column per gene"
                )
            obs = obs_input.copy(deep=True) if isinstance(obs_input, pd.DataFrame) else pd.DataFrame(obs_input)
            if len(obs) != X.shape[0]:
                raise ValueError("obs length must match expression_matrix rows")
            n_cells, n_genes = X.shape
            data = SCData(X=X, obs=obs, var=pd.DataFrame({"gene_name": gene_names}))

        gene_name_to_idx = {g: i for i, g in enumerate(gene_names)}

        # Contract Parameters
        compound_name = str(params.get("compound_name", "Anti_Inflammatory_Small_Molecule"))
        dosage = float(params.get("dosage", params.get("scale_factor", 1.0)))
        n_perms = int(params.get("n_permutations", 500))
        custom_drug_sig = params.get("drug_signature", None)
        custom_disease_sig = params.get("disease_signature", None)
        signature_conditions = None

        # State labels are required for a state transition result.  They are
        # also a valid source for a disease signature when no separate
        # condition column is present (for example Homeostatic vs DAM states).
        state_col = (
            "microglia_state" if "microglia_state" in obs.columns
            else ("cluster" if "cluster" in obs.columns else ("condition" if "condition" in obs.columns else None))
        )
        if state_col is None and "state_labels" in params:
            state_labels = params["state_labels"]
            if len(state_labels) != n_cells:
                raise ValueError("state_labels length must match expression matrix rows")
            obs = obs.copy(deep=True)
            obs["state_labels"] = list(state_labels)
            state_col = "state_labels"

        # 1. Determine Disease Signature s_disease (1 x N_genes)
        disease_sig = np.zeros(n_genes, dtype=np.float32)
        disease_signature_source = None
        disease_signature_measured = False
        if custom_disease_sig is not None:
            if isinstance(custom_disease_sig, dict):
                unknown = set(custom_disease_sig).difference(gene_name_to_idx)
                if len(unknown) == len(custom_disease_sig):
                    raise ValueError("disease_signature contains no genes present in the input")
                for g, val in custom_disease_sig.items():
                    if g in gene_name_to_idx:
                        disease_sig[gene_name_to_idx[g]] = float(val)
            else:
                disease_sig_array = np.asarray(custom_disease_sig, dtype=np.float32)
                if disease_sig_array.ndim != 1 or len(disease_sig_array) != n_genes:
                    raise ValueError("disease_signature must have one value per input gene")
                disease_sig = disease_sig_array
            if not np.isfinite(disease_sig).all():
                raise ValueError("disease_signature contains non-finite values")
            disease_signature_source = str(params.get("disease_signature_source", "user_provided_unverified"))
            disease_signature_measured = bool(params.get("disease_signature_measured", False))
            signature_conditions = {"source": disease_signature_source}
        else:
            cond_col = "condition" if "condition" in obs.columns else None
            cond_a = params.get("condition_a")
            cond_b = params.get("condition_b")
            if cond_col is not None:
                conditions = list(obs[cond_col].dropna().unique())
                if cond_a is None or cond_b is None:
                    normalized = {str(value).strip().lower(): value for value in conditions}
                    disease_keys = ("ad", "disease", "case", "treated", "cko", "mutant")
                    control_keys = ("control", "ctrl", "healthy", "wt", "wildtype", "con")
                    cond_a = next((normalized[key] for key in disease_keys if key in normalized), None)
                    cond_b = next((normalized[key] for key in control_keys if key in normalized), None)
                if cond_a is not None and cond_b is not None and cond_a != cond_b:
                    mask_a = (obs[cond_col] == cond_a).values
                    mask_b = (obs[cond_col] == cond_b).values
                    if mask_a.any() and mask_b.any():
                        mean_a = np.mean(X[mask_a], axis=0)
                        mean_b = np.mean(X[mask_b], axis=0)
                        disease_sig = np.log2((mean_a + 1e-3) / (mean_b + 1e-3)).astype(np.float32)
                        disease_signature_source = "observed_condition_contrast"
                        disease_signature_measured = True
                        signature_conditions = {"condition_a": str(cond_a), "condition_b": str(cond_b)}

            if disease_signature_source is None and state_col is not None:
                labels = obs[state_col].astype(str)
                disease_mask = labels.str.lower().map(
                    lambda value: any(token in value for token in ("dam", "disease", "case", "ad", "m3", "cko", "mutant"))
                ).to_numpy()
                control_mask = labels.str.lower().map(
                    lambda value: any(token in value for token in ("homeo", "control", "healthy", "wt", "m1", "con"))
                ).to_numpy()
                if disease_mask.any() and control_mask.any():
                    disease_sig = np.log2(
                        (np.mean(X[disease_mask], axis=0) + 1e-3)
                        / (np.mean(X[control_mask], axis=0) + 1e-3)
                    ).astype(np.float32)
                    disease_signature_source = "observed_state_contrast"
                    disease_signature_measured = True
                    signature_conditions = {
                        "disease_states": sorted(labels[disease_mask].unique().tolist()),
                        "control_states": sorted(labels[control_mask].unique().tolist()),
                    }

            if disease_signature_source is None:
                raise ValueError(
                    "Disease signature requires an explicit disease_signature or two observed, semantically labelled groups; "
                    "a single unlabeled cohort cannot define disease"
                )

        # 2. Determine Drug Perturbation Signature s_drug (1 x N_genes)
        drug_sig = np.zeros(n_genes, dtype=np.float32)
        drug_signature_source = None
        drug_signature_verified = False
        if custom_drug_sig is not None:
            if isinstance(custom_drug_sig, dict):
                unknown = set(custom_drug_sig).difference(gene_name_to_idx)
                if len(unknown) == len(custom_drug_sig):
                    raise ValueError("drug_signature contains no genes present in the input")
                for g, val in custom_drug_sig.items():
                    if g in gene_name_to_idx:
                        drug_sig[gene_name_to_idx[g]] = float(val)
            else:
                drug_sig_array = np.asarray(custom_drug_sig, dtype=np.float32)
                if drug_sig_array.ndim != 1 or len(drug_sig_array) != n_genes:
                    raise ValueError("drug_signature must have one value per input gene")
                drug_sig = drug_sig_array
            if not np.isfinite(drug_sig).all():
                raise ValueError("drug_signature contains non-finite values")
            drug_signature_source = str(params.get("drug_signature_source", "user_provided_unverified"))
            drug_signature_verified = bool(params.get("drug_signature_verified", False))
        elif compound_name in REFERENCE_COMPOUND_DATABASE:
            ref_dict = REFERENCE_COMPOUND_DATABASE[compound_name]
            for g, val in ref_dict.items():
                if g in gene_name_to_idx:
                    drug_sig[gene_name_to_idx[g]] = float(val)
            drug_signature_source = "local_reference_unverified"
            drug_signature_verified = False
        else:
            raise ValueError("Unknown compound requires an explicit drug_signature; no reversal signature is synthesized")

        if not np.isfinite(drug_sig).all() or np.linalg.norm(drug_sig) < 1e-8:
            raise ValueError("Drug signature has no finite values for genes present in the input")

        # 3. Compute CMAP Cosine Discordance Score
        reversal_score, cosine_sim, p_val = compute_cmap_cosine_discordance(
            disease_sig=disease_sig,
            drug_sig=drug_sig,
            n_permutations=n_perms,
            random_seed=int(params.get("random_seed", 42)),
        )
        signature_provenance_verified = bool(
            disease_signature_measured and drug_signature_verified
        )
        therapeutic_potential = bool(
            signature_provenance_verified and reversal_score > 0.30 and p_val < 0.05
        )

        # 4. Simulate Counterfactual Treatment & Cell State Transitions
        # X_drug = max(0, X + dosage * drug_sig)
        # Signatures are log2 fold-change profiles, so apply them as fold
        # changes to the measured baseline rather than adding them as raw
        # counts.  The operation remains explicitly counterfactual.
        fold_change = np.power(2.0, dosage * drug_sig, dtype=np.float32)
        X_drug = np.maximum(0.0, X * fold_change[np.newaxis, :]).astype(np.float32)

        # Identify Cell States (e.g. microglia_state, cluster, or condition)
        if state_col is None:
            raise ValueError(
                "Compound counterfactual transitions require an observed state column "
                "(microglia_state, cluster, condition) or explicit state_labels"
            )
        state_values = obs[state_col].astype(str)
        unique_states = list(dict.fromkeys(state_values.dropna().tolist()))
        if len(unique_states) < 2:
            raise ValueError("At least two observed states are required for transition inference")

        if int(n_cells) * int(max(1, len(unique_states))) * int(n_genes) > MAX_COMPOUND_ELEMENTS:
            raise ValueError(
                f"Compound transition probabilities exceed the {MAX_COMPOUND_ELEMENTS:,}-element resource limit; reduce cells/states/genes"
            )

        n_states = len(unique_states)
        state_to_idx = {s: i for i, s in enumerate(unique_states)}

        # Compute baseline centroids for each state
        centroids = []
        for s in unique_states:
            mask = (state_values == s).values
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
            u_mask = (state_values == u_state).values
            if np.sum(u_mask) > 0:
                T_mat[u_idx, :] = np.mean(trans_probs[u_mask], axis=0)
            else:
                T_mat[u_idx, u_idx] = 1.0

        # Compute disease-to-healthy transition rate
        disease_transition_rate = None
        disease_transition_available = False
        dam_states = [str(params["source_state"])] if params.get("source_state") is not None else []
        homeo_states = [str(params["target_state"])] if params.get("target_state") is not None else []
        if any(s not in state_to_idx for s in dam_states + homeo_states):
            raise ValueError("Requested transition states must occur in the observed state column")

        if dam_states and homeo_states:
            dam_u = state_to_idx[dam_states[0]]
            homeo_v = state_to_idx[homeo_states[0]]
            disease_transition_rate = float(T_mat[dam_u, homeo_v])
            disease_transition_available = True

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
            "disease_transition_available": disease_transition_available,
            "signature_conditions": signature_conditions,
            "disease_signature_source": disease_signature_source,
            "disease_signature_measured": disease_signature_measured,
            "drug_signature_source": drug_signature_source,
            "drug_signature_verified": drug_signature_verified,
            "simulation_mode": simulation_mode,
            "input_artifact_uri": in_uri,
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
                "therapeutic_potential": therapeutic_potential,
                "therapeutic_signal_unverified": bool(reversal_score > 0.30 and p_val < 0.05 and not signature_provenance_verified),
                "disease_signature_source": disease_signature_source,
                "disease_signature_measured": disease_signature_measured,
                "drug_signature_source": drug_signature_source,
                "drug_signature_verified": drug_signature_verified,
                "signature_provenance_verified": signature_provenance_verified,
                "disease_transition_available": disease_transition_available,
                "simulation_mode": simulation_mode,
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
                "disease_transition_available": disease_transition_available,
                "disease_signature_source": disease_signature_source,
                "disease_signature_measured": disease_signature_measured,
                "drug_signature_source": drug_signature_source,
                "drug_signature_verified": drug_signature_verified,
                "signature_provenance_verified": signature_provenance_verified,
                "simulation_mode": simulation_mode,
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
                "therapeutic_potential": therapeutic_potential,
                "therapeutic_signal_unverified": bool(reversal_score > 0.30 and p_val < 0.05 and not signature_provenance_verified),
                "disease_signature_source": disease_signature_source,
                "disease_signature_measured": disease_signature_measured,
                "drug_signature_source": drug_signature_source,
                "drug_signature_verified": drug_signature_verified,
                "disease_transition_available": disease_transition_available,
                "simulation_mode": simulation_mode,
            }
        )


def generate_compound_perturbation_evidence(
    contract: TaskContract,
    result: TaskResult,
    compound_name: str,
    reversal_score: float,
    transition_rate: Optional[float],
    p_value: float = 0.01,
) -> EvidenceNode:
    """
    Generates a calibrated EvidenceType.PERTURBATION node from small molecule response simulation.
    Strictly caps in silico causal confidence score at <= 0.50.
    """
    out_uris = result.output_artifacts
    task_id = contract.task_id

    transition_available = bool(
        result.metrics.get("disease_transition_available", transition_rate is not None)
    )
    disease_source = result.metrics.get("disease_signature_source", "unspecified")
    drug_source = result.metrics.get("drug_signature_source", "unspecified")
    source_verified = bool(result.metrics.get("signature_provenance_verified", False))
    if transition_rate is None or not transition_available:
        return EvidenceNode(
            evidence_id=f"E_compound_{compound_name}_{task_id}",
            type=EvidenceType.PERTURBATION,
            polarity=EvidencePolarity.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
            score=0.0,
            summary=(
                f"In silico compound simulation with {compound_name} produced a transcriptomic "
                "counterfactual, but disease-to-homeostatic transition was not estimable "
                "from the available state labels."
            ),
            source_task_id=task_id,
            source_artifact_uris=out_uris,
            metrics={
                "compound_name": compound_name,
                "reversal_score": reversal_score,
                "transition_rate": None,
                "disease_transition_available": False,
                "disease_signature_source": disease_source,
                "drug_signature_source": drug_source,
                "source_verified": source_verified,
                "in_silico_confidence_cap": 0.50,
            },
            biological_context={
                "compound": compound_name,
                "causal_status": "in_silico_perturbed",
                "state_claim": "not_estimable",
            },
            is_simulated=True,
            source_verified=source_verified,
        )

    # Quantitative score calibrated and capped at 0.50
    normalized_score = max(0.10, min(0.50, float(max(0.0, reversal_score))))

    strength = (
        EvidenceStrength.STRONG if reversal_score >= 0.50
        else (EvidenceStrength.MODERATE if reversal_score >= 0.20 else EvidenceStrength.WEAK)
    )

    polarity = EvidencePolarity.SUPPORTING if reversal_score > 0 else EvidencePolarity.CONTRADICTING
    action = "therapeutic reversal" if reversal_score > 0 else "disease exacerbation"

    verification_note = "verified signature sources" if source_verified else "unverified local/provided signature sources"
    summary = (
        f"In silico compound simulation with {compound_name} predicts {action} of disease signature "
        f"(CMAP discordance score: {reversal_score:.2f}, p-val: {p_value:.3f}, "
        f"counterfactual homeostatic transition rate: {transition_rate*100:.1f}%; {verification_note})."
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
            "disease_signature_source": disease_source,
            "drug_signature_source": drug_source,
            "source_verified": source_verified,
        },
        biological_context={
            "compound": compound_name,
            "causal_status": "in_silico_perturbed",
            "signature_source_status": "verified" if source_verified else "unverified",
        },
        is_simulated=True,
        source_verified=source_verified,
    )
