"""
Capability Router: explicit overrides, local policies and registry lookup.
Routes spatial capabilities, agent adapters, knowledge engine, and perturbation simulations.
"""

from typing import Dict, Any, Optional
from eacbp.schemas.task import TaskContract
from eacbp.schemas.study import StudyManifest
from eacbp.capabilities.registry import CapabilityRegistry


class CapabilityRouter:
    """Resolves target capabilities and specific methods using 3-tier hierarchical routing."""

    def __init__(self, capability_registry: CapabilityRegistry):
        self.registry = capability_registry

    def resolve_method(
        self,
        capability_name: str,
        manifest: StudyManifest,
        current_state: Dict[str, Any],
    ) -> str:
        override = current_state.get("method_overrides", {}).get(capability_name)
        if override:
            return self.registry.get(capability_name, override).implementation_id
        if capability_name == "deg" and current_state.get("advanced_analysis"):
            return self.registry.get("deg", "pydeseq2_pseudobulk_v1").implementation_id
        profile = current_state.get("method_profile", "baseline")
        if profile not in {"baseline", "standard"}:
            raise ValueError("method_profile must be 'baseline' or 'standard'")
        if profile == "standard":
            if capability_name == "integration" and len(current_state.get("batches", manifest.experimental_design.batches)) < 2:
                return "no_correction_v1"
            standard = {"integration": "harmonypy_v1", "clustering": "scanpy_leiden_umap_v1",
                        "trajectory_inference": "scanpy_dpt_v1"}
            if capability_name in standard:
                return self.registry.get(capability_name, standard[capability_name]).implementation_id
        # Layer 1: Hard constraints & data specifications
        if capability_name in ("quantification", "fastq_quantification"):
            return current_state.get("quant_tool", "kb_python_v1")

        if capability_name == "trajectory_inference":
            return "root_distance_pseudotime_v1"

        # Layer 2: Method policies based on data audit & metadata
        if capability_name == "integration":
            batch_possible = current_state.get("batch_effect_possible", True)
            batches = current_state.get("batches", manifest.experimental_design.batches)
            if batch_possible and len(batches) > 1:
                return "batch_mean_centering_v1"
            return "no_correction_v1"

        if capability_name == "deg":
            # The implementation determines eligibility from actual donor metadata.
            return "donor_pseudobulk_welch_v1"

        if capability_name == "qc":
            return "sc_qc_v1"

        if capability_name == "clustering":
            return "kmeans_marker_embedding_v1"

        if capability_name == "subset_cells":
            return "subset_cells_v1"

        if capability_name == "dataset_audit":
            return "sc_audit_v1"

        if capability_name == "differential_abundance":
            return "state_abundance_donor_welch_v1"

        if capability_name == "normalization":
            return "sc_normalize_log1p_v1"

        # Spatial Plane
        if capability_name == "spatial_domain":
            return "spatial_domain_knn_v1"

        if capability_name == "spatial_deg":
            return "spatial_deg_morans_i_v1"

        if capability_name in ("cell_cell_communication", "spatial_cci"):
            return "cci_ligand_receptor_v1"

        # External Agent Adapters
        if capability_name in ("spacell_microenvironment_analysis", "spacell_analysis", "adapter_spacell"):
            return "spacell_agent_v1"

        if capability_name in ("chatcell_dialogue_prediction", "chatcell_reasoning", "adapter_chatcell"):
            return "chatcell_agent_v1"

        if capability_name in ("gene_function_reasoning", "gene_agent_reasoning", "adapter_geneagent"):
            return "gene_agent_v1"

        # Multi-Source Knowledge Engine
        if capability_name == "knowledge_retrieval":
            is_prior_guided = (
                manifest.analysis_policy.prior_guided_analysis
                or bool(manifest.hypotheses.user_provided)
                or current_state.get("prior_guided", False)
            )
            return "knowledge_engine_prior_v1" if is_prior_guided else "knowledge_engine_discovery_v1"

        # In Silico Perturbation Simulation Plane
        if capability_name in ("genetic_perturbation_simulation", "genetic_perturbation"):
            ptype = current_state.get("perturbation_type", "knockout")
            if ptype == "overexpression":
                return "in_silico_overexpression_v1"
            return "in_silico_crispr_ko_v1"

        if capability_name in ("compound_perturbation_simulation", "compound_perturbation"):
            return "in_silico_compound_response_v1"

        # Layer 3: Registry lookup fallback
        if capability_name in self.registry._capabilities:
            cap = self.registry.get(capability_name)
            return cap.implementation_id

        # If not registered, raise informative error
        raise KeyError(
            f"Unregistered capability '{capability_name}' in CapabilityRegistry. "
            f"Available capabilities: {list(self.registry._capabilities.keys())}"
        )
