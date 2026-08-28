"""
Unit tests for Evidence DAG, 4-Tier Language protocol, and Claim Engine.
"""

import pytest
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
    LanguageTier,
    ClaimType,
)
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.evidence.language import LanguageEnforcer
from eacbp.evidence.claim import ClaimEngine


def test_confidence_calculator_with_contradictions():
    sup1 = EvidenceNode(
        evidence_id="E1",
        type=EvidenceType.PSEUDOBULK_DEG,
        strength=EvidenceStrength.STRONG,
        score=0.9,
        summary="Apoe is upregulated in AD.",
        source_task_id="task_008",
    )
    sup2 = EvidenceNode(
        evidence_id="E2",
        type=EvidenceType.TRAJECTORY_STABILITY,
        strength=EvidenceStrength.STRONG,
        score=0.85,
        summary="Trajectory consistently points from homeostatic to DAM.",
        source_task_id="task_009",
    )

    conf_clean = ConfidenceCalculator.calculate([sup1, sup2], [])
    assert conf_clean.association > 0.6
    assert conf_clean.mechanistic > 0.6
    assert conf_clean.overall > 0.6

    # Add contradicting evidence
    contra = EvidenceNode(
        evidence_id="E3",
        type=EvidenceType.TRAJECTORY_STABILITY,
        polarity=EvidencePolarity.CONTRADICTING,
        strength=EvidenceStrength.STRONG,
        score=0.8,
        summary="Alternative trajectory suggests reverse progression.",
        source_task_id="task_010",
    )
    conf_contradicted = ConfidenceCalculator.calculate([sup1, sup2], [contra])
    assert conf_contradicted.overall < conf_clean.overall


def test_language_enforcer_catches_unsupported_causality():
    # Observational study claiming causal verb 'drives' at Level 1 or 2
    valid, msg = LanguageEnforcer.audit_statement(
        statement="APOE expression drives microglial neuroinflammation.",
        tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE,
        causal_status="observational"
    )
    assert not valid
    assert "EPISTEMIC VIOLATION" in msg

    # Qualified statement passes
    valid_qual, _ = LanguageEnforcer.audit_statement(
        statement="APOE expression is significantly associated with microglial state transition.",
        tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE,
        causal_status="observational"
    )
    assert valid_qual


def test_claim_engine_creates_traceable_claims():
    eg = EvidenceGraph()
    engine = ClaimEngine(evidence_graph=eg)

    ev1 = EvidenceNode(
        evidence_id="E201",
        type=EvidenceType.DIFFERENTIAL_ABUNDANCE,
        strength=EvidenceStrength.STRONG,
        score=0.88,
        summary="M3 state is enriched in AD.",
        source_task_id="task_007",
    )
    eg.add_evidence(ev1)

    claim = engine.create_claim(
        claim_id="C103",
        statement="APOE-high microglia may represent an AD-associated transitional state.",
        language_tier=LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION,
        claim_type=ClaimType.MECHANISTIC_HYPOTHESIS,
        causal_status="observational",
        support_evidence_ids=["E201"],
    )

    assert claim.claim_id == "C103"
    assert claim.confidence.overall > 0.0
    assert "E201" in claim.support_evidence_ids
    
    prov = eg.get_claim_provenance("C103")
    assert len(prov["supporting_evidence"]) == 1
    assert prov["supporting_evidence"][0]["evidence_id"] == "E201"
