"""
Intent Parser translating natural language research queries into structured StudyManifest and scientific objectives.
"""

from typing import Dict, Any, List, Optional
import re
from eacbp.schemas.study import (
    StudyManifest,
    BiologicalDesign,
    ExperimentalDesign,
    DataSpec,
    Hypotheses,
    AnalysisPolicy,
    Constraints,
    ReproducibilityConfig,
)


class IntentParser:
    """Parses scientific questions into structured study objectives without premature algorithm locking."""

    @staticmethod
    def parse_prompt_to_manifest(
        user_prompt: str,
        study_id: str = "AD_mouse_001",
        raw_artifact_uri: Optional[str] = None,
    ) -> StudyManifest:
        prompt_lower = user_prompt.lower()

        # Species detection
        species = "mus_musculus" if any(w in prompt_lower for w in ["mouse", "小鼠", "murine"]) else "homo_sapiens"

        # Tissue detection
        tissue = "brain" if any(w in prompt_lower for w in ["brain", "脑", "cortex", "hippocampus"]) else "tissue_sample"

        # Disease / Intervention detection
        if any(w in prompt_lower for w in ["kat8", "mof", "myst1", "cko", "敲除"]):
            disease = "Kat8_cKO_Developmental_Disruption"
            conditions = ["cKO", "con"]
        elif any(w in prompt_lower for w in ["ad", "alzheimer", "阿尔茨海默", "痴呆"]):
            disease = "Alzheimer"
            conditions = ["AD", "control"]
        else:
            disease = "healthy"
            conditions = ["condition_1", "condition_2"]

        # Target cell types
        target_cells = []
        if any(w in prompt_lower for w in ["microglia", "小胶质"]):
            target_cells.append("Microglia")
        if any(w in prompt_lower for w in ["astrocyte", "星形胶质"]):
            target_cells.append("Astrocytes")
        if any(w in prompt_lower for w in ["neuron", "神经元"]):
            target_cells.append("Neurons")
        if any(w in prompt_lower for w in ["progenitor", "前体", "干细胞", "stem"]):
            target_cells.append("Progenitors")

        # Modalities
        modalities = ["scRNA"]
        has_spatial = any(w in prompt_lower for w in ["spatial", "空间", "visium", "stereoseq", "merfish"])
        if has_spatial:
            modalities.append("spatial")
        if any(w in prompt_lower for w in ["perturb", "crispr", "ko", "knockout", "敲除", "cko"]):
            modalities.append("perturbation")

        # Prior-guided mode detection
        is_kat8 = any(w in prompt_lower for w in ["kat8", "mof", "myst1", "h4k16ac", "组蛋白乙酰化"])
        is_dam = any(w in prompt_lower for w in ["dam假说", "prior", "基于", "trem2-apoe", "trem2", "apoe"])
        is_prior_guided = is_kat8 or is_dam or ("prior" in prompt_lower)
        user_hypotheses = []
        if is_kat8:
            user_hypotheses.append("Kat8 (Mof) H4K16ac epigenetic loss impairs cell cycle progression and lineage differentiation")
        elif is_dam:
            user_hypotheses.append("DAM subpopulation regulation via TREM2-APOE axis")

        manifest = StudyManifest(
            study_id=study_id,
            title=f"Single-Cell Study: {disease} in {species} {tissue}",
            biological_design=BiologicalDesign(
                species=species,
                tissue=tissue,
                disease=disease,
                conditions=conditions,
                target_cell_types=target_cells or ["Microglia"],
            ),
            experimental_design=ExperimentalDesign(
                biological_unit="mouse" if species == "mus_musculus" else "donor",
                batches=["batch_1", "batch_2"],
                total_samples=12,
                donor_replicates_per_condition={"AD": 6, "control": 6},
            ),
            data=DataSpec(
                modalities=modalities,
                raw_artifact_uri=raw_artifact_uri or f"adata://{study_id}/raw/v1",
                has_spatial_coordinates=has_spatial,
                has_rna_velocity=False,
            ),
            hypotheses=Hypotheses(
                user_provided=user_hypotheses,
            ),
            analysis_policy=AnalysisPolicy(
                discovery_mode=not is_prior_guided,
                prior_guided_analysis=is_prior_guided,
                strict_reproducibility=True,
                prefer_pseudobulk=True,
            ),
            reproducibility=ReproducibilityConfig(
                random_seed=42,
            ),
        )
        return manifest
