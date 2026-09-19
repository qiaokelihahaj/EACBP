"""
Confidence score computation for evidence and multi-modal claims across 5 pillars.
"""

from typing import List
from eacbp.schemas.evidence import (
    EvidenceNode,
    ConfidenceScore,
    EvidenceStrength,
    EvidencePolarity,
    EvidenceType,
)


class ConfidenceCalculator:
    """
    Calculates multidimensional confidence scores across 5 pillars:
    1. Spatial localization (++++, alpha = 1.0)
    2. Pseudobulk DEG (++++, alpha = 1.0)
    3. Literature support (++++, alpha = 0.95)
    4. Trajectory inference (+++, alpha = 0.75)
    5. In silico perturbation (+, alpha = 0.35, causal capped at 0.50)
    """

    STRENGTH_WEIGHTS = {
        EvidenceStrength.VERY_STRONG: 1.0,
        EvidenceStrength.STRONG: 0.8,
        EvidenceStrength.MODERATE: 0.6,
        EvidenceStrength.WEAK: 0.3,
        EvidenceStrength.INSUFFICIENT: 0.0,
    }

    PILLAR_WEIGHTS = {
        EvidenceType.SPATIAL_LOCALIZATION: 1.00,
        EvidenceType.PSEUDOBULK_DEG: 1.00,
        EvidenceType.CELL_LEVEL_DEG: 0.60,
        EvidenceType.DIFFERENTIAL_ABUNDANCE: 0.85,
        EvidenceType.DATASET_AUDIT: 0.80,
        EvidenceType.QC_METRICS: 0.80,
        EvidenceType.CLUSTERING_STABILITY: 0.80,
        EvidenceType.LITERATURE_SUPPORT: 0.95,
        EvidenceType.PATHWAY_ENRICHMENT: 0.90,
        EvidenceType.TRAJECTORY_STABILITY: 0.75,
        EvidenceType.ROOT_SENSITIVITY: 0.75,
        EvidenceType.PERTURBATION: 0.35,
        EvidenceType.FUNCTIONAL_ACTIVITY: 0.90,
        EvidenceType.SENSITIVITY_ANALYSIS: 0.80,
        EvidenceType.CELL_COMMUNICATION: 0.60,
        EvidenceType.CELL_ANNOTATION: 0.60,
    }

    @classmethod
    def calculate(
        cls,
        supporting_evidence: List[EvidenceNode],
        contradicting_evidence: List[EvidenceNode],
    ) -> ConfidenceScore:
        if not supporting_evidence and not contradicting_evidence:
            return ConfidenceScore(association=0.0, mechanistic=0.0, causal=0.0, overall=0.0)

        assoc_scores = []
        mech_scores = []
        causal_scores = []

        # Collapse overlapping input origins within each dimension. A conservative
        # minimum prevents repeated analyses of one matrix from increasing a
        # dimension's score. This remains a heuristic, not a probability model.
        grouped = {"association": [], "mechanistic": [], "causal": []}
        for ev in supporting_evidence:
            # Strength weight * score
            base_score = cls.STRENGTH_WEIGHTS.get(ev.strength, 0.6) * ev.score
            
            # Categorize by evidence nature
            if ev.type in (
                EvidenceType.DIFFERENTIAL_ABUNDANCE,
                EvidenceType.PSEUDOBULK_DEG,
                EvidenceType.CELL_LEVEL_DEG,
                EvidenceType.QC_METRICS,
                EvidenceType.DATASET_AUDIT,
                EvidenceType.CLUSTERING_STABILITY,
                EvidenceType.CELL_ANNOTATION,
            ):
                dimension = "association"
            elif ev.type in (
                EvidenceType.TRAJECTORY_STABILITY,
                EvidenceType.ROOT_SENSITIVITY,
                EvidenceType.PATHWAY_ENRICHMENT,
                EvidenceType.LITERATURE_SUPPORT,
                EvidenceType.SPATIAL_LOCALIZATION,
                EvidenceType.FUNCTIONAL_ACTIVITY,
                EvidenceType.CELL_COMMUNICATION,
            ):
                dimension = "mechanistic"
            elif ev.type == EvidenceType.PERTURBATION:
                dimension = "causal"
            else:
                # Sensitivity describes reliability, and supplies no extra
                # independent support or causal credit.
                continue
            origins = set(ev.data_origin_uris) or {"task:" + ev.source_task_id}
            score = base_score
            remaining = []
            for existing_origins, existing_score in grouped[dimension]:
                if origins & existing_origins:
                    origins |= existing_origins
                    score = min(score, existing_score)
                else:
                    remaining.append((existing_origins, existing_score))
            # A newly merged group can bridge previously disjoint roots.
            changed = True
            while changed:
                changed = False
                kept = []
                for existing_origins, existing_score in remaining:
                    if origins & existing_origins:
                        origins |= existing_origins
                        score = min(score, existing_score)
                        changed = True
                    else:
                        kept.append((existing_origins, existing_score))
                remaining = kept
            grouped[dimension] = remaining + [(origins, score)]

        assoc_scores = [score for _, score in grouped["association"]]
        mech_scores = [score for _, score in grouped["mechanistic"]]
        causal_scores = [score for _, score in grouped["causal"]]

        # Average dimensions
        assoc_conf = float(sum(assoc_scores) / len(assoc_scores)) if assoc_scores else 0.0
        mech_conf = float(sum(mech_scores) / len(mech_scores)) if mech_scores else 0.0
        
        # In silico causal confidence is capped at 0.50
        raw_causal = float(sum(causal_scores) / len(causal_scores)) if causal_scores else 0.0
        causal_conf = min(0.50, raw_causal)

        # Contradiction penalty
        contra_penalty = 0.0
        for cev in contradicting_evidence:
            contra_penalty += cls.STRENGTH_WEIGHTS.get(cev.strength, 0.4) * 0.25

        # Weighted overall composite
        if causal_conf > 0:
            raw_overall = (0.35 * assoc_conf) + (0.35 * mech_conf) + (0.30 * causal_conf)
        else:
            raw_overall = (0.55 * assoc_conf) + (0.45 * mech_conf)

        final_overall = max(0.0, min(1.0, raw_overall - contra_penalty))

        return ConfidenceScore(
            association=round(min(1.0, max(0.0, assoc_conf)), 3),
            mechanistic=round(min(1.0, max(0.0, mech_conf)), 3),
            causal=round(min(1.0, max(0.0, causal_conf)), 3),
            overall=round(final_overall, 3),
        )
