"""
Claim Engine orchestrating scientific claim synthesis, validation, and multi-evidence linking.
"""

from typing import List, Optional
from eacbp.schemas.evidence import (
    ClaimNode,
    ClaimType,
    LanguageTier,
    EvidenceNode,
)
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.evidence.language import LanguageEnforcer


class ClaimEngine:
    """Creates, validates, and links scientific claims with multi-hop evidence traceability."""

    def __init__(self, evidence_graph: Optional[EvidenceGraph] = None):
        self.evidence_graph = evidence_graph or EvidenceGraph()

    def create_claim(
        self,
        claim_id: str,
        statement: str,
        language_tier: LanguageTier,
        claim_type: ClaimType = ClaimType.MECHANISTIC_HYPOTHESIS,
        causal_status: str = "observational",
        support_evidence_ids: Optional[List[str]] = None,
        contradiction_evidence_ids: Optional[List[str]] = None,
    ) -> ClaimNode:
        support_evidence_ids = support_evidence_ids or []
        contradiction_evidence_ids = contradiction_evidence_ids or []

        # 1. Epistemic language validation
        valid_lang, lang_error = LanguageEnforcer.audit_statement(
            statement=statement,
            tier=language_tier,
            causal_status=causal_status,
        )
        if not valid_lang:
            raise ValueError(lang_error)

        # 2. Gather EvidenceNodes
        sup_nodes = [
            self.evidence_graph.evidence_nodes[eid]
            for eid in support_evidence_ids
            if eid in self.evidence_graph.evidence_nodes
        ]
        contra_nodes = [
            self.evidence_graph.evidence_nodes[eid]
            for eid in contradiction_evidence_ids
            if eid in self.evidence_graph.evidence_nodes
        ]

        # 3. Calculate multidimensional confidence
        confidence = ConfidenceCalculator.calculate(sup_nodes, contra_nodes)

        # 4. Generate provenance summary
        ev_summaries = [f"[{e.evidence_id}: {e.summary}]" for e in sup_nodes]
        prov_summary = (
            f"Claim {claim_id} supported by {len(sup_nodes)} evidence items: "
            + "; ".join(ev_summaries)
        )

        claim = ClaimNode(
            claim_id=claim_id,
            statement=statement,
            language_tier=language_tier,
            claim_type=claim_type,
            causal_status=causal_status,
            support_evidence_ids=support_evidence_ids,
            contradiction_evidence_ids=contradiction_evidence_ids,
            confidence=confidence,
            provenance_summary=prov_summary,
        )

        self.evidence_graph.add_claim(claim)
        return claim
