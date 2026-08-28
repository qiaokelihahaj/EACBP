"""
Evidence and Claim schema defining the 4-tier language protocol and evidence-backed claims.
"""

from typing import List, Dict, Any, Optional
from enum import Enum
from datetime import datetime, timezone
from pydantic import BaseModel, Field


class EvidenceType(str, Enum):
    DATASET_AUDIT = "dataset_audit"
    QC_METRICS = "qc_metrics"
    INTEGRATION_EVALUATION = "integration_evaluation"
    CLUSTERING_STABILITY = "clustering_stability"
    DIFFERENTIAL_ABUNDANCE = "differential_abundance"
    PSEUDOBULK_DEG = "pseudobulk_deg"
    CELL_LEVEL_DEG = "cell_level_deg"
    TRAJECTORY_STABILITY = "trajectory_stability"
    ROOT_SENSITIVITY = "root_sensitivity"
    SPATIAL_LOCALIZATION = "spatial_localization"
    PATHWAY_ENRICHMENT = "pathway_enrichment"
    LITERATURE_SUPPORT = "literature_support"
    PERTURBATION = "perturbation"


class EvidencePolarity(str, Enum):
    SUPPORTING = "supporting"
    CONTRADICTING = "contradicting"
    NEUTRAL = "neutral"


class EvidenceStrength(str, Enum):
    VERY_STRONG = "very_strong"
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"
    INSUFFICIENT = "insufficient"


class LanguageTier(str, Enum):
    LEVEL_1_OBSERVATION = "Level 1 - Observation"
    LEVEL_2_STATISTICAL_INFERENCE = "Level 2 - Statistical inference"
    LEVEL_3_SUPPORTED_INTERPRETATION = "Level 3 - Supported interpretation"
    LEVEL_4_HYPOTHESIS = "Level 4 - Hypothesis"


class ClaimType(str, Enum):
    DESCRIPTIVE = "descriptive"
    COMPARATIVE = "comparative"
    STATE_TRANSITION = "state_transition"
    REGULATORY = "regulatory"
    MECHANISTIC_HYPOTHESIS = "mechanistic_hypothesis"


class ConfidenceScore(BaseModel):
    association: float = Field(0.0, ge=0.0, le=1.0, description="Observational & statistical strength")
    mechanistic: float = Field(0.0, ge=0.0, le=1.0, description="Pathway & regulatory consistency")
    causal: float = Field(0.0, ge=0.0, le=1.0, description="Perturbation or counterfactual confirmation")
    overall: float = Field(0.0, ge=0.0, le=1.0, description="Aggregated confidence metric")


class EvidenceNode(BaseModel):
    evidence_id: str = Field(..., description="Unique evidence ID, e.g. E203")
    type: EvidenceType
    polarity: EvidencePolarity = Field(EvidencePolarity.SUPPORTING)
    strength: EvidenceStrength = Field(EvidenceStrength.MODERATE)
    score: float = Field(0.5, ge=0.0, le=1.0, description="Normalized quantitative score (0-1)")
    summary: str = Field(..., description="Human-readable summary of the evidence")
    
    # Provenance
    source_task_id: str = Field(..., description="Task ID where this evidence was computed")
    source_artifact_uris: List[str] = Field(default_factory=list, description="Artifact URIs backing this evidence")
    
    # Quantitative support details
    metrics: Dict[str, Any] = Field(default_factory=dict, description="e.g. {'p_val_adj': 1e-6, 'log2fc': 1.8, 'n_replicates': 6}")
    biological_context: Dict[str, Any] = Field(default_factory=dict, description="e.g. {'cell_type': 'microglia', 'gene': 'Apoe'}")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))


class ClaimNode(BaseModel):
    claim_id: str = Field(..., description="Unique claim ID, e.g. C103")
    statement: str = Field(..., description="The scientific statement")
    language_tier: LanguageTier = Field(LanguageTier.LEVEL_2_STATISTICAL_INFERENCE)
    claim_type: ClaimType = Field(ClaimType.MECHANISTIC_HYPOTHESIS)
    causal_status: str = Field("observational", description="'observational' vs 'experimental_perturbed'")
    
    # Evidence linkages
    support_evidence_ids: List[str] = Field(default_factory=list)
    contradiction_evidence_ids: List[str] = Field(default_factory=list)
    
    confidence: ConfidenceScore = Field(default_factory=ConfidenceScore)
    provenance_summary: str = Field("", description="Traceable summary of evidence path")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
