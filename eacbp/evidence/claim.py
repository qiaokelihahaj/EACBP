"""
Claim Engine orchestrating scientific claim synthesis, validation, and multi-evidence linking.
"""

from typing import List, Optional
import math
from eacbp.schemas.evidence import (
    ClaimNode,
    ClaimType,
    LanguageTier,
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    ConfidenceScore,
)
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.evidence.language import LanguageEnforcer


class ClaimEngine:
    """Creates, validates, and links scientific claims with multi-hop evidence traceability."""

    def __init__(self, evidence_graph: Optional[EvidenceGraph] = None):
        self.evidence_graph = evidence_graph or EvidenceGraph()

    @staticmethod
    def has_valid_statistics(node: EvidenceNode) -> bool:
        """A score or a raw p value is not a multiple-testing-adjusted result."""
        if not node.audit_passed or not node.source_artifact_uris:
            return False
        if node.type not in (EvidenceType.PSEUDOBULK_DEG, EvidenceType.DIFFERENTIAL_ABUNDANCE,
                             EvidenceType.SPATIAL_LOCALIZATION, EvidenceType.PATHWAY_ENRICHMENT, EvidenceType.FUNCTIONAL_ACTIVITY):
            return False
        alpha = 0.05
        threshold_key = next(
            (key for key in ("alpha", "significance_alpha") if key in node.metrics),
            None,
        )
        if threshold_key is not None:
            try:
                alpha = float(node.metrics[threshold_key])
            except (TypeError, ValueError):
                return False
            if not math.isfinite(alpha) or not 0.0 < alpha < 1.0:
                return False
        for key in ("fdr_q_value", "fdr", "p_val_adj"):
            try:
                value = float(node.metrics[key])
            except (KeyError, TypeError, ValueError):
                continue
            if math.isfinite(value) and 0 <= value < alpha:
                return True
        return False

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
        if claim_type == ClaimType.RESULT_SUMMARY:
            raise ValueError("Use create_result_summary for an audited neutral result")
        support_evidence_ids = support_evidence_ids or []
        contradiction_evidence_ids = contradiction_evidence_ids or []
        missing = set(support_evidence_ids + contradiction_evidence_ids) - self.evidence_graph.evidence_nodes.keys()
        if missing:
            raise ValueError(f"Unknown evidence IDs: {sorted(missing)}")
        if not support_evidence_ids:
            raise ValueError("A claim requires at least one supporting evidence item.")

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
        if any(e.polarity != EvidencePolarity.SUPPORTING for e in sup_nodes):
            raise ValueError("Neutral or contradicting evidence cannot be registered as support.")
        if any(not e.audit_passed or not e.source_artifact_uris for e in sup_nodes):
            raise ValueError("Supporting evidence must have passed audit and reference an artifact.")
        if language_tier in (LanguageTier.LEVEL_2_STATISTICAL_INFERENCE, LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION):
            statistical_support = [e for e in sup_nodes if self.has_valid_statistics(e)]
            if not statistical_support:
                raise ValueError("Statistical conclusions require audited significant FDR evidence.")
            # Arbitrary prose cannot inherit significance from an unrelated result.
            # Until a typed interpretation schema exists, admit exact audited
            # result statements only, rather than guessing semantic equivalence.
            if statement.removeprefix("[SIMULATED DATA] ") not in [e.summary for e in statistical_support]:
                raise ValueError("Statistical claim must restate its supporting evidence summary exactly")
            statistical_support = [e for e in statistical_support
                                   if e.summary == statement.removeprefix("[SIMULATED DATA] ")]
        if language_tier == LanguageTier.LEVEL_3_SUPPORTED_INTERPRETATION:
            if not any(e.source_verified and e.type in (EvidenceType.LITERATURE_SUPPORT, EvidenceType.PATHWAY_ENRICHMENT) for e in sup_nodes):
                raise ValueError("Supported interpretation requires verified biological/knowledge support.")
            knowledge = [e for e in sup_nodes if e.source_verified and e.type in (EvidenceType.LITERATURE_SUPPORT, EvidenceType.PATHWAY_ENRICHMENT)]
            matched = any(
                any(key in stat.biological_context and key in source.biological_context
                    and stat.biological_context[key] == source.biological_context[key]
                    for key in ("gene", "state", "pathway"))
                and all(stat.biological_context[key] == source.biological_context[key]
                        for key in stat.biological_context.keys() & source.biological_context.keys())
                for stat in statistical_support for source in knowledge)
            if not matched:
                raise ValueError("Verified knowledge must match the statistical evidence biological context")
        simulated = any(e.is_simulated for e in sup_nodes)
        if simulated and "[SIMULATED DATA]" not in statement:
            statement = "[SIMULATED DATA] " + statement

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
            is_simulated=simulated,
            inference_contract_ids=sorted({e.inference_contract_id for e in sup_nodes + contra_nodes if e.inference_contract_id}),
        )

        self.evidence_graph.add_claim(claim)
        return claim

    def create_result_summary(self, claim_id: str, evidence_id: str) -> ClaimNode:
        """Describe an audited result without promoting it to a positive discovery."""
        node = self.evidence_graph.evidence_nodes.get(evidence_id)
        if (node is None or node.type != EvidenceType.STATISTICAL_RESULT
                or node.polarity != EvidencePolarity.NEUTRAL or not node.audit_passed
                or not node.source_artifact_uris or not node.inference_contract_id):
            raise ValueError("A result summary requires audited neutral statistical evidence and its inference contract")
        if node.metrics.get("scientific_status") not in {"estimated_supported", "estimated_inconclusive", "not_estimable"}:
            raise ValueError("An unsupported scientific state cannot be admitted as a result summary")
        statement = node.summary
        if node.is_simulated and not statement.startswith("[SIMULATED DATA] "):
            statement = "[SIMULATED DATA] " + statement
        claim = ClaimNode(
            claim_id=claim_id, statement=statement,
            language_tier=LanguageTier.LEVEL_1_OBSERVATION,
            claim_type=ClaimType.RESULT_SUMMARY,
            support_evidence_ids=[evidence_id],
            inference_contract_ids=[node.inference_contract_id],
            confidence=ConfidenceScore(), is_simulated=node.is_simulated,
            provenance_summary=f"Audited result description from {evidence_id}; no positive or equivalence conclusion is implied.",
        )
        self.evidence_graph.add_claim(claim)
        return claim
