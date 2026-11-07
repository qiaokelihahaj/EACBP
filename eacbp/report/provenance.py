"""
Provenance Engine: Maps every scientific statement back to Claims, Evidence nodes, Tasks, and Artifacts.
"""

from typing import Dict, Any, List
from eacbp.schemas.task import TaskResult
from eacbp.evidence.graph import EvidenceGraph
from eacbp.evidence.extraction import validate_evidence_row_binding
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
        
        def resolve_evidence(eid):
            ev = self.evidence_graph.evidence_nodes.get(eid)
            if not ev:
                raise ValueError(f"Claim {claim_id} references missing evidence {eid}")

            task = self.task_map.get(ev.source_task_id)
            if not ev.audit_passed or task is None or getattr(task.status, "value", task.status) != "success":
                raise ValueError(f"Evidence {eid} lacks a successful audited source task")
            if not ev.source_artifact_uris:
                raise ValueError(f"Evidence {eid} has no source artifacts")
            if ev.inference_contract_id != getattr(task, "inference_contract_id", None):
                raise ValueError(f"Evidence {eid} inference contract differs from its source task")
            row_status = validate_evidence_row_binding(ev, self.artifact_registry, task)
            artifact_details = []
            for uri in ev.source_artifact_uris:
                if not self.artifact_registry.exists(uri) or uri not in task.output_artifacts:
                    raise ValueError(f"Evidence {eid} references missing or foreign artifact {uri}")
                else:
                    meta = self.artifact_registry.get_metadata(uri)
                    if meta.created_by_task != ev.source_task_id:
                        raise ValueError(f"Evidence {eid} artifact was created by a different task")
                    lineage_path = self.artifact_registry.lineage.get_lineage_path_from_root(uri)
                    ancestors = set(self.artifact_registry.lineage.get_ancestors(uri)) | {uri}
                    lineage_edges = [(a, b) for a, b in self.artifact_registry.lineage.graph.edges
                                     if a in ancestors and b in ancestors]
                    artifact_details.append({
                        "uri": uri,
                        "type": meta.type.value,
                        "operation": meta.operation,
                        "sha256": meta.sha256_hash,
                        "parameters": meta.parameters,
                        "software": meta.software_versions,
                        "random_seed": meta.random_seed,
                        "lineage_path": lineage_path,
                        "ancestor_artifacts": sorted(ancestors),
                        "lineage_edges": lineage_edges,
                    })

            return {
                "evidence_id": ev.evidence_id,
                "evidence_type": ev.type.value,
                "strength": ev.strength.value,
                "score": ev.score,
                "summary": ev.summary,
                "inference_contract_id": ev.inference_contract_id,
                "metrics": ev.metrics,
                "biological_context": ev.biological_context,
                "row_selector": ev.row_selector.model_dump(mode="json") if ev.row_selector else None,
                "row_verification": row_status,
                "source_task": {
                    "task_id": task.task_id if task else ev.source_task_id,
                    "capability": task.capability if task else "unknown",
                    "method": task.method_used if task else "unknown",
                    "executed_operations": task.executed_operations if task else [],
                },
                "artifacts": artifact_details,
            }

        evidence_chain = [resolve_evidence(eid) for eid in claim.support_evidence_ids]
        contradiction_chain = [resolve_evidence(eid) for eid in claim.contradiction_evidence_ids]

        return {
            "claim_id": claim.claim_id,
            "statement": claim.statement,
            "language_tier": claim.language_tier.value,
            "causal_status": claim.causal_status,
            "confidence": claim.confidence.model_dump(),
            "evidence_count": len(evidence_chain),
            "inference_contract_ids": list(claim.inference_contract_ids),
            "evidence_chain": evidence_chain,
            "contradiction_evidence_count": len(contradiction_chain),
            "contradiction_evidence_chain": contradiction_chain,
        }
