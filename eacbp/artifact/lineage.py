"""
Lineage graph tracking data provenance, parent-child transitions, and branch diffs.
"""

from typing import Dict, List, Any, Optional, Set
import networkx as nx
from eacbp.schemas.artifact import ArtifactMetadata, LineageNode


class LineageGraph:
    """Directed Acyclic Graph (DAG) maintaining data provenance across all artifact versions."""

    def __init__(self):
        self.graph = nx.DiGraph()
        self.metadata_store: Dict[str, ArtifactMetadata] = {}

    def add_artifact(self, metadata: ArtifactMetadata):
        uri = metadata.uri
        self.metadata_store[uri] = metadata
        self.graph.add_node(
            uri,
            artifact_id=metadata.artifact_id,
            type=metadata.type.value,
            operation=metadata.operation,
            task_id=metadata.created_by_task,
            hash=metadata.sha256_hash,
            created_at=metadata.created_at.isoformat(),
        )

        for parent_uri in metadata.parent_uris:
            if not self.graph.has_node(parent_uri):
                # Add placeholder node if parent hasn't been explicitly registered yet
                self.graph.add_node(parent_uri)
            self.graph.add_edge(parent_uri, uri, task_id=metadata.created_by_task, operation=metadata.operation)

    def get_parents(self, uri: str) -> List[str]:
        if not self.graph.has_node(uri):
            return []
        return list(self.graph.predecessors(uri))

    def get_children(self, uri: str) -> List[str]:
        if not self.graph.has_node(uri):
            return []
        return list(self.graph.successors(uri))

    def get_ancestors(self, uri: str) -> List[str]:
        if not self.graph.has_node(uri):
            return []
        return list(nx.ancestors(self.graph, uri))

    def get_descendants(self, uri: str) -> List[str]:
        if not self.graph.has_node(uri):
            return []
        return list(nx.descendants(self.graph, uri))

    def get_lineage_path_from_root(self, uri: str) -> List[str]:
        """Returns the shortest path from any root artifact to the given uri."""
        if not self.graph.has_node(uri):
            return [uri]
        roots = [n for n in self.graph.nodes if self.graph.in_degree(n) == 0]
        for root in roots:
            if nx.has_path(self.graph, root, uri):
                return list(nx.shortest_path(self.graph, root, uri))
        return [uri]

    def compare_branches(self, uri_a: str, uri_b: str) -> Dict[str, Any]:
        """Compares two branched artifacts, finding their common ancestor and operation diffs."""
        ancestors_a = set(self.get_ancestors(uri_a)).union({uri_a})
        ancestors_b = set(self.get_ancestors(uri_b)).union({uri_b})
        common_ancestors = ancestors_a.intersection(ancestors_b)

        meta_a = self.metadata_store.get(uri_a)
        meta_b = self.metadata_store.get(uri_b)

        return {
            "branch_a": uri_a,
            "branch_b": uri_b,
            "common_ancestors": list(common_ancestors),
            "operation_a": meta_a.operation if meta_a else "unknown",
            "operation_b": meta_b.operation if meta_b else "unknown",
            "parameters_a": meta_a.parameters if meta_a else {},
            "parameters_b": meta_b.parameters if meta_b else {},
            "summary_metrics_a": meta_a.summary_metrics if meta_a else {},
            "summary_metrics_b": meta_b.summary_metrics if meta_b else {},
        }

    def to_mermaid(self) -> str:
        """Generates a Mermaid graph representation of the lineage DAG."""
        lines = ["graph TD"]
        for node in self.graph.nodes:
            meta = self.metadata_store.get(node)
            label = f"{node}<br><i>{meta.operation}</i>" if meta else node
            lines.append(f'    "{node}"["{label}"]')
        for src, dst, data in self.graph.edges(data=True):
            op = data.get("operation", "")
            lines.append(f'    "{src}" -->|{op}| "{dst}"')
        return "\n".join(lines)
