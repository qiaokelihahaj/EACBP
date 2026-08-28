"""
Provenance Engine: Maps every scientific statement back to Claims, Evidence nodes, Tasks, and Artifacts.
"""

from typing import Dict, Any, List, Optional
from eacbp.schemas.evidence import ClaimNode, EvidenceNode
from eacbp.schemas.artifact import ArtifactMetadata
from eacbp.schemas.task import TaskResult
from eacbp.evidence.graph import EvidenceGraph
from eacbp.artifact.registry import ArtifactRegistry


class SentenceProvenanceTracker:
    """Provides complete multi-hop sentence-to-raw-data provenance resolution."""

    def __init__(self, evidence_graph: EvidenceGraph, artifact_registry: ArtifactRegistry, task_history: List[TaskResult]):
        self.evidence_graph = evidence_graph
        self.artifact_registry = artifact_registry
        self.task_map = {t.task_id: t for t in task_history}

    def resolve_claim_provenance_card(self, claim_id: str) -> Dict[str, Any]:
        """Resolves the complete trace for a specific claim: Claim -> Evidence -> Task -> Artifact -> Raw."""
        if claim_id not in self.evidence_graph.claim_nodes:
            raise KeyError(f"Claim ID '{claim_id}' not found.")

        claim = self.evidence_graph.claim_nodes[claim_id]
        
        evidence_chain = []
        for eid in claim.support_evidence_ids:
            ev = self.evidence_graph.evidence_nodes.get(eid)
            if not ev:
                continue

            task = self.task_map.get(ev.source_task_id)
            artifact_details = []
            for uri in ev.source_artifact_uris:
                if self.artifact_registry.exists(uri):
                    meta = self.artifact_registry.get_metadata(uri)
                    lineage_path = self.artifact_registry.lineage.get_lineage_path_from_root(uri)
                    artifact_details.append({
                        "uri": uri,
                        "type": meta.type.value,
                        "operation": meta.operation,
                        "sha256": meta.sha256_hash,
                        "parameters": meta.parameters,
                        "software": meta.software_versions,
                        "random_seed": meta.random_seed,
                        "lineage_path": lineage_path,
                    })

            evidence_chain.append({
                "evidence_id": ev.evidence_id,
                "evidence_type": ev.type.value,
                "strength": ev.strength.value,
                "score": ev.score,
                "summary": ev.summary,
                "source_task": {
                    "task_id": task.task_id if task else ev.source_task_id,
                    "capability": task.capability if task else "unknown",
                    "method": task.method_used if task else "unknown",
                    "executed_operations": task.executed_operations if task else [],
                },
                "artifacts": artifact_details,
            })

        return {
            "claim_id": claim.claim_id,
            "statement": claim.statement,
            "language_tier": claim.language_tier.value,
            "causal_status": claim.causal_status,
            "confidence": claim.confidence.model_dump(),
            "evidence_count": len(evidence_chain),
            "evidence_chain": evidence_chain,
        }
