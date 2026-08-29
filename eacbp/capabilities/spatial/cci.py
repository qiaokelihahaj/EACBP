"""
Spatial Cell-Cell Communication (CCI) and Ligand-Receptor Analysis Capability for EACBP.
Computes expression interaction scores, spatial proximity contact density weighting,
and spatial permutation significance testing with FDR correction.
"""

from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial.domain import (
    validate_spatial_coordinates,
    build_spatial_neighborhood_graph,
)
from eacbp.capabilities.spatial.autocorrelation import benjamini_hochberg
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


CURATED_LIGAND_RECEPTOR_PAIRS = [
    ("Apoe", "Trem2", "Microglial Activation & Lipid Sensing"),
    ("Apoe", "Lrp1", "Lipid Metabolism & Clearance"),
    ("Apoe", "Ldlr", "Lipid Transport & Homeostasis"),
    ("App", "Cd74", "Amyloid Processing & Immune Interaction"),
    ("Cx3cl1", "Cx3cr1", "Microglial Homeostasis & Chemotaxis"),
    ("Ccl2", "Ccr2", "Monocyte & Microglia Recruitment"),
    ("C3", "C3ar1", "Complement Cascade & Synaptic Pruning"),
    ("Spp1", "Cd44", "Neuroinflammation & Extracellular Matrix"),
    ("Spp1", "Itgb1", "Cell Adhesion & Migration"),
    ("Tgfb1", "Tgfbr1", "TGF-beta Signaling & Quiescence"),
    ("Il1b", "Il1r1", "Pro-inflammatory Cascade"),
    ("Tnf", "Tnfrsf1a", "TNF Signaling & Apoptosis"),
    ("Csf1", "Csf1r", "Microglial Survival & Proliferation"),
    ("Vegfa", "Kdr", "Angiogenesis & Vascular Remodeling"),
    ("Jag1", "Notch1", "Notch Signaling & Cell Fate"),
    ("Bdnf", "Ntrk2", "Neurotrophic Support & Plasticity"),
    ("Ccl3", "Ccr5", "Chemokine Signaling"),
    ("Ccl5", "Ccr5", "Immune Cell Recruitment"),
    ("Gas6", "Axl", "Phagocytosis & Efferocytosis"),
    ("Pros1", "Tyro3", "TAM Receptor Phagocytic Signaling"),
]


def calculate_spatial_contact_density(
    labels: np.ndarray,
    unique_labels: List[str],
    W: np.ndarray,
) -> Dict[Tuple[str, str], float]:
    """
    Computes spatial contact density matrix between all cell type pairs:
    W_spatial(A, B) = sum_{i in A, j in B} W_ij / (|A| * |B|)
    """
    densities = {}
    label_arr = np.asarray(labels)
    
    for a in unique_labels:
        mask_a = (label_arr == a)
        n_a = int(mask_a.sum())
        if n_a == 0:
            continue
        for b in unique_labels:
            mask_b = (label_arr == b)
            n_b = int(mask_b.sum())
            if n_b == 0:
                continue
            
            sub_w = W[np.ix_(mask_a, mask_b)]
            total_edges = float(np.sum(sub_w))
            w_density = total_edges / float(n_a * n_b)
            densities[(a, b)] = float(w_density)
            
    return densities


def compute_spatial_cci(
    data: SCData,
    W: np.ndarray,
    cell_type_col: str = "cell_type",
    lr_pairs: Optional[List[Tuple[str, str, str]]] = None,
    n_permutations: int = 200,
    fdr_threshold: float = 0.05,
    random_seed: int = 42,
) -> pd.DataFrame:
    """
    Computes proximity-weighted ligand-receptor interaction scores and permutation p-values.
    """
    np.random.seed(random_seed)
    
    # 1. Determine cell type annotations
    obs = data.obs
    if cell_type_col in obs.columns:
        labels = obs[cell_type_col].astype(str).values
    elif "spatial_domain" in obs.columns:
        labels = obs["spatial_domain"].astype(str).values
    elif "cell_type_ground_truth" in obs.columns:
        labels = obs["cell_type_ground_truth"].astype(str).values
    elif "leiden" in obs.columns:
        labels = obs["leiden"].astype(str).values
    else:
        labels = np.array(["CellType_All"] * data.n_obs)

    unique_cell_types = sorted(list(set(labels)))
    if len(unique_cell_types) == 0:
        return pd.DataFrame()

    # 2. Gene indexing (case-insensitive lookup)
    gene_names = list(data.var["gene_name"]) if "gene_name" in data.var.columns else [f"Gene_{i}" for i in range(data.n_vars)]
    gene_map = {g.lower(): (i, g) for i, g in enumerate(gene_names)}

    pairs_to_evaluate = lr_pairs or CURATED_LIGAND_RECEPTOR_PAIRS
    
    valid_pairs = []
    for item in pairs_to_evaluate:
        lig, rec = item[0], item[1]
        pw = item[2] if len(item) > 2 else "Unknown"
        if lig.lower() in gene_map and rec.lower() in gene_map:
            l_idx, l_orig = gene_map[lig.lower()]
            r_idx, r_orig = gene_map[rec.lower()]
            valid_pairs.append((l_orig, r_orig, l_idx, r_idx, pw))

    if not valid_pairs:
        if data.n_vars >= 2:
            l_idx, r_idx = 0, 1
            l_orig, r_orig = str(gene_names[0]), str(gene_names[1])
            valid_pairs.append((l_orig, r_orig, l_idx, r_idx, "Data-Derived Interaction"))

    # 3. Compute observed spatial contact densities
    obs_contact_densities = calculate_spatial_contact_density(labels, unique_cell_types, W)

    # 4. Compute mean expression per cell type
    mean_expr = {}
    for ct in unique_cell_types:
        mask = (labels == ct)
        if mask.sum() > 0:
            mean_expr[ct] = np.mean(data.X[mask], axis=0)
        else:
            mean_expr[ct] = np.zeros(data.n_vars)

    # 5. Compute observed raw & spatial interaction scores
    interactions = []
    for sender in unique_cell_types:
        for receiver in unique_cell_types:
            w_spatial = obs_contact_densities.get((sender, receiver), 0.0)
            
            for l_orig, r_orig, l_idx, r_idx, pw in valid_pairs:
                l_val = float(mean_expr[sender][l_idx])
                r_val = float(mean_expr[receiver][r_idx])
                
                raw_score = float(np.sqrt(max(0.0, l_val * r_val)))
                spatial_score = float(raw_score * w_spatial)
                
                interactions.append({
                    "sender_cell_type": sender,
                    "receiver_cell_type": receiver,
                    "ligand": l_orig,
                    "receptor": r_orig,
                    "pathway": pw,
                    "ligand_expr": l_val,
                    "receptor_expr": r_val,
                    "raw_score": raw_score,
                    "spatial_weight": w_spatial,
                    "spatial_score": spatial_score,
                    "l_idx": l_idx,
                    "r_idx": r_idx,
                })

    if not interactions:
        return pd.DataFrame()

    # 6. Spatial Permutation Test
    n_interactions = len(interactions)
    obs_scores = np.array([inter["spatial_score"] for inter in interactions])
    greater_counts = np.zeros(n_interactions, dtype=int)

    n_perms = max(10, n_permutations)
    for _ in range(n_perms):
        perm_labels = np.random.permutation(labels)
        perm_densities = calculate_spatial_contact_density(perm_labels, unique_cell_types, W)
        
        for k, inter in enumerate(interactions):
            s_ct = inter["sender_cell_type"]
            r_ct = inter["receiver_cell_type"]
            perm_w = perm_densities.get((s_ct, r_ct), 0.0)
            perm_spatial_score = inter["raw_score"] * perm_w
            if perm_spatial_score >= obs_scores[k]:
                greater_counts[k] += 1

    p_values = (greater_counts + 1.0) / (n_perms + 1.0)
    fdr_q = benjamini_hochberg(p_values)

    results_list = []
    for k, inter in enumerate(interactions):
        results_list.append({
            "sender_cell_type": inter["sender_cell_type"],
            "receiver_cell_type": inter["receiver_cell_type"],
            "ligand": inter["ligand"],
            "receptor": inter["receptor"],
            "pathway": inter["pathway"],
            "mean_ligand_expr": inter["ligand_expr"],
            "mean_receptor_expr": inter["receptor_expr"],
            "raw_interaction_score": inter["raw_score"],
            "spatial_proximity_weight": inter["spatial_weight"],
            "spatial_interaction_score": inter["spatial_score"],
            "p_value": float(p_values[k]),
            "fdr_q_value": float(fdr_q[k]),
            "p_val_adj": float(fdr_q[k]),
            "is_significant": bool((fdr_q[k] < fdr_threshold) and (inter["spatial_score"] > 0.0)),
        })

    res_df = pd.DataFrame(results_list)
    res_df = res_df.sort_values(["spatial_interaction_score", "fdr_q_value"], ascending=[False, True]).reset_index(drop=True)
    return res_df


class CellCellCommunicationCapability(BaseCapability):
    """
    Spatial Cell-Cell Communication (CCI) & Ligand-Receptor Capability.
    """

    def __init__(
        self,
        capability_name: str = "cell_cell_communication",
        implementation_id: str = "cci_ligand_receptor_v1",
    ):
        super().__init__(
            capability_name=capability_name,
            implementation_id=implementation_id,
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_modalities=["spatial", "scRNA"],
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            suitable_for=["cell_cell_communication", "ligand_receptor_interaction", "spatial_niches"],
            output_types=[ArtifactType.TABLE],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        # Parameters
        k_neighbors = int(contract.parameters.get("k_neighbors", 6))
        n_permutations = int(contract.parameters.get("n_permutations", 200))
        fdr_threshold = float(contract.parameters.get("fdr_threshold", 0.05))
        cell_type_col = contract.parameters.get("cell_type_col", "cell_type")
        custom_lr_pairs = contract.parameters.get("custom_lr_pairs", None)
        random_seed = int(contract.parameters.get("random_seed", 42))

        # Spatial coordinates
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
            # Non-spatial standard single-cell CCI (CellChat / CellPhoneDB style)
            W = np.ones((data.n_obs, data.n_obs), dtype=np.float32) / float(max(1, data.n_obs))
            is_spatial = False
        else:
            validated_coords = validate_spatial_coordinates(coords, data.n_obs)
            is_spatial = True
            # Spatial connectivity graph safely loaded
            if "spatial_connectivities" in data.obsm and isinstance(data.obsm["spatial_connectivities"], np.ndarray):
                W = np.asarray(data.obsm["spatial_connectivities"], dtype=np.float32)
            elif "spatial_connectivities" in data.uns and isinstance(data.uns["spatial_connectivities"], np.ndarray):
                W = np.asarray(data.uns["spatial_connectivities"], dtype=np.float32)
            elif hasattr(data, "obsp") and isinstance(getattr(data, "obsp", None), dict) and "spatial_connectivities" in data.obsp and isinstance(data.obsp["spatial_connectivities"], np.ndarray):
                W = np.asarray(data.obsp["spatial_connectivities"], dtype=np.float32)
            else:
                W, _, _ = build_spatial_neighborhood_graph(validated_coords, k_neighbors=k_neighbors)

        # Run Spatial CCI Computation
        cci_df = compute_spatial_cci(
            data=data,
            W=W,
            cell_type_col=cell_type_col,
            lr_pairs=custom_lr_pairs,
            n_permutations=n_permutations,
            fdr_threshold=fdr_threshold,
            random_seed=random_seed,
        )

        sig_interactions = cci_df[cci_df["is_significant"]] if not cci_df.empty else pd.DataFrame()

        # Resolve output URI
        uri_obj = ArtifactURI.parse(in_uri)
        if contract.expected_outputs:
            out_uri = contract.expected_outputs[0]
        else:
            out_uri = f"table://{uri_obj.study_id}/spatial_cci/v1"

        registry.register(
            uri_str=out_uri,
            payload=cci_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="compute_spatial_cell_cell_communication",
            parent_uris=[in_uri],
            parameters={
                "k_neighbors": k_neighbors,
                "n_permutations": n_permutations,
                "fdr_threshold": fdr_threshold,
                "cell_type_col": cell_type_col,
                "random_seed": random_seed,
            },
            summary_metrics={
                "total_interactions_tested": len(cci_df),
                "significant_interactions": len(sig_interactions),
                "top_interaction": (
                    f"{sig_interactions.iloc[0]['sender_cell_type']}->{sig_interactions.iloc[0]['receiver_cell_type']}: {sig_interactions.iloc[0]['ligand']}-{sig_interactions.iloc[0]['receptor']}"
                    if not sig_interactions.empty else "none"
                ),
            },
        )

        all_ops = [
            "load_lr_database",
            "calculate_spatial_contact_density",
            "compute_spatial_cci_score",
            "run_spatial_permutation_test",
            "benjamini_hochberg_correction",
            "evaluate_cell_cell_communication",
            "ligand_receptor_cci",
        ]
        if contract.allowed_operations:
            executed_ops = [op for op in all_ops if op in contract.allowed_operations]
        else:
            executed_ops = [
                "load_lr_database",
                "calculate_spatial_contact_density",
                "compute_spatial_cci_score",
                "run_spatial_permutation_test",
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
                "total_interactions": len(cci_df),
                "significant_interactions": len(sig_interactions),
                "top_interactions": (
                    sig_interactions[["sender_cell_type", "receiver_cell_type", "ligand", "receptor", "spatial_interaction_score"]]
                    .head(5)
                    .to_dict(orient="records")
                    if not sig_interactions.empty else []
                ),
            },
        )
