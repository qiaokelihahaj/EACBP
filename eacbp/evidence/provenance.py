"""Artifact-owned provenance rules shared by all evidence extractors.

Extractors propose scientific observations. They cannot choose the data roots,
remove a simulation marker, or certify their own literature sources.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

from eacbp.artifact.uri import ArtifactURI
from eacbp.auditor.base import ValidationReport, ValidationSeverity
from eacbp.schemas.artifact import ArtifactMetadata
from eacbp.schemas.evidence import EvidenceNode
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


@dataclass(frozen=True)
class ArtifactProvenance:
    root_uris: tuple[str, ...]
    is_simulated: bool


def resolve_artifact_provenance(
    source_uris: Iterable[str],
    metadata_lookup: Callable[[str], ArtifactMetadata],
) -> ArtifactProvenance:
    """Resolve every ancestor, rejecting missing metadata and lineage cycles.

    A metadata lookup can come from a live registry or a validated snapshot.
    No payload is loaded and no registry state is modified.
    """
    roots: set[str] = set()
    state: dict[str, int] = {}
    simulated = False
    pending = [(ArtifactURI.parse(uri).to_string(), False) for uri in source_uris]
    while pending:
        uri, complete = pending.pop()
        if complete:
            state[uri] = 2
            continue
        if state.get(uri) == 2:
            continue
        if state.get(uri) == 1:
            raise ValueError(f"Evidence artifact lineage contains a cycle at {uri}")
        metadata = metadata_lookup(uri)
        if metadata.uri != uri:
            raise ValueError(f"Evidence artifact metadata URI does not match {uri}")
        simulated = simulated or bool(metadata.summary_metrics.get("is_simulated"))
        state[uri] = 1
        pending.append((uri, True))
        parents = [ArtifactURI.parse(parent).to_string() for parent in metadata.parent_uris]
        if parents:
            pending.extend((parent, False) for parent in parents)
        else:
            roots.add(uri)
    return ArtifactProvenance(tuple(sorted(roots)), simulated)


def normalize_evidence_candidates(
    candidates: Iterable[EvidenceNode],
    contract: TaskContract,
    result: TaskResult,
    report: ValidationReport,
    registry,
) -> list[EvidenceNode]:
    """Admit candidate observations with provenance derived from artifacts.

    Identity mismatches are rejected rather than repaired. A candidate may
    conservatively add a simulation marker, but cannot erase one recorded by
    the computation or any artifact ancestor. Source verification is rejected
    until an independent, verifiable source-receipt mechanism exists; a DOI or
    an extractor-provided boolean is not such a receipt.
    """
    if report.target_task_id != contract.task_id:
        raise ValueError("Evidence audit must refer to the current task")
    if result.task_id != contract.task_id or result.capability != contract.capability:
        raise ValueError("Evidence result must refer to the current task and capability")
    if result.status != TaskStatus.SUCCESS:
        return []
    if (not report.overall_passed or report.stop_rule_triggered
            or any(not check.passed and check.severity in (
                ValidationSeverity.ERROR, ValidationSeverity.STOP_RULE,
            ) for check in report.checks)):
        return []

    output_uris = {ArtifactURI.parse(uri).to_string() for uri in result.output_artifacts}
    seen_ids: set[str] = set()
    admitted = []
    # Reuse metadata across candidates from the same task, including large
    # DEG tables whose evidence nodes all cite the same output.
    metadata_cache: dict[str, ArtifactMetadata] = {}
    provenance_cache: dict[tuple[str, ...], ArtifactProvenance] = {}

    def metadata_lookup(uri: str) -> ArtifactMetadata:
        if uri not in metadata_cache:
            metadata_cache[uri] = registry.get_metadata(uri)
        return metadata_cache[uri]

    for candidate in candidates:
        if not isinstance(candidate, EvidenceNode):
            raise TypeError("Evidence extractors must return EvidenceNode objects")
        node = EvidenceNode.model_validate(candidate.model_dump())
        if not node.evidence_id or node.evidence_id in seen_ids:
            raise ValueError("Evidence IDs must be non-empty and unique within a task")
        if node.source_task_id != contract.task_id:
            raise ValueError("Evidence must cite the current task as source_task_id")
        sources = tuple(sorted({ArtifactURI.parse(uri).to_string()
                                for uri in node.source_artifact_uris}))
        if not sources or not set(sources).issubset(output_uris):
            raise ValueError("Evidence must cite only current task output artifacts")
        for uri in sources:
            if metadata_lookup(uri).created_by_task != contract.task_id:
                raise ValueError("Evidence output artifact was created by a different task")
        if node.source_verified:
            raise ValueError("source_verified requires independent source verification; no trusted source receipt is available")
        if sources not in provenance_cache:
            provenance_cache[sources] = resolve_artifact_provenance(sources, metadata_lookup)
        provenance = provenance_cache[sources]
        node.source_task_id = contract.task_id
        node.source_artifact_uris = list(sources)
        node.data_origin_uris = list(provenance.root_uris)
        node.is_simulated = bool(node.is_simulated or provenance.is_simulated
                                 or result.metrics.get("is_simulated"))
        node.audit_passed = True
        seen_ids.add(node.evidence_id)
        admitted.append(node)
    return admitted
