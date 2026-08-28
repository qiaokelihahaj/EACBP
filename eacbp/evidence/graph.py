"""
Evidence DAG (Directed Acyclic Graph) managing scientific evidence, claims, and provenance traces.
"""

from typing import Dict, List, Any, Optional
import networkx as nx
from eacbp.schemas.evidence import EvidenceNode, ClaimNode, EvidencePolarity


class EvidenceGraph:
    """The Evidence DAG maintaining all evidence units, claims, and their multi-hop provenance."""

    def __init__(self):
        self.graph = nx.DiGraph()
        self.evidence_nodes: Dict[str, EvidenceNode] = {}
        self.claim_nodes: Dict[str, ClaimNode] = {}

    def add_evidence(self, evidence: EvidenceNode):
        eid = evidence.evidence_id
        self.evidence_nodes[eid] = evidence
        self.graph.add_node(
            eid,
            node_type="evidence",
            evidence_type=evidence.type.value,
            polarity=evidence.polarity.value,
            strength=evidence.strength.value,
            score=evidence.score,
            summary=evidence.summary,
            source_task=evidence.source_task_id,
            source_artifacts=evidence.source_artifact_uris,
        )

    def add_claim(self, claim: ClaimNode):
        cid = claim.claim_id
        self.claim_nodes[cid] = claim
        self.graph.add_node(
            cid,
            node_type="claim",
            statement=claim.statement,
            language_tier=claim.language_tier.value,
            claim_type=claim.claim_type.value,
            causal_status=claim.causal_status,
            confidence=claim.confidence.overall,
        )

        for eid in claim.support_evidence_ids:
            if self.graph.has_node(eid):
                self.graph.add_edge(eid, cid, relationship="supports")

        for eid in claim.contradiction_evidence_ids:
            if self.graph.has_node(eid):
                self.graph.add_edge(eid, cid, relationship="contradicts")

    def link_evidence_to_claim(self, evidence_id: str, claim_id: str, is_support: bool = True):
        if not self.graph.has_node(evidence_id) or not self.graph.has_node(claim_id):
            raise KeyError(f"Nodes '{evidence_id}' or '{claim_id}' not found in EvidenceGraph")
        rel = "supports" if is_support else "contradicts"
        self.graph.add_edge(evidence_id, claim_id, relationship=rel)
        
        claim = self.claim_nodes[claim_id]
        if is_support and evidence_id not in claim.support_evidence_ids:
            claim.support_evidence_ids.append(evidence_id)
        elif not is_support and evidence_id not in claim.contradiction_evidence_ids:
            claim.contradiction_evidence_ids.append(evidence_id)

    def get_claim_provenance(self, claim_id: str) -> Dict[str, Any]:
        """Traces a claim back to all supporting and contradicting evidence nodes and their source artifacts."""
        if claim_id not in self.claim_nodes:
            raise KeyError(f"Claim ID '{claim_id}' not found in EvidenceGraph")
        
        claim = self.claim_nodes[claim_id]
        support_items = []
        for eid in claim.support_evidence_ids:
            ev = self.evidence_nodes.get(eid)
            if ev:
                support_items.append(ev.model_dump())

        contra_items = []
        for eid in claim.contradiction_evidence_ids:
            ev = self.evidence_nodes.get(eid)
            if ev:
                contra_items.append(ev.model_dump())

        return {
            "claim_id": claim.claim_id,
            "statement": claim.statement,
            "language_tier": claim.language_tier.value,
            "causal_status": claim.causal_status,
            "confidence": claim.confidence.model_dump(),
            "supporting_evidence": support_items,
            "contradicting_evidence": contra_items,
        }

    def to_mermaid(self) -> str:
        """Generates Mermaid diagram representing the Evidence-to-Claim DAG."""
        lines = ["graph TD"]
        for cid, claim in self.claim_nodes.items():
            conf = claim.confidence.overall
            lines.append(f'    "{cid}"["Claim {cid}: {claim.statement}<br><b>Tier: {claim.language_tier.value} | Conf: {conf:.2f}</b>"]')

        for eid, ev in self.evidence_nodes.items():
            lines.append(f'    "{eid}"["Evidence {eid}: {ev.summary}<br><i>{ev.type.value} ({ev.strength.value})</i>"]')

        for src, dst, data in self.graph.edges(data=True):
            rel = data.get("relationship", "linked")
            lines.append(f'    "{src}" -->|{rel}| "{dst}"')

        return "\n".join(lines)
