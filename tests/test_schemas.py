"""
Unit tests for EACBP Schemas.
"""

import pytest
from eacbp.schemas.study import StudyManifest, BiologicalDesign, ExperimentalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, RetryPolicy
from eacbp.schemas.artifact import ArtifactMetadata, ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType, ClaimNode, LanguageTier, ConfidenceScore


def test_study_manifest_schema():
    manifest = StudyManifest(
        study_id="AD_test_001",
        biological_design=BiologicalDesign(
            species="mus_musculus",
            tissue="brain",
            disease="Alzheimer",
            conditions=["AD", "control"],
            target_cell_types=["Microglia"],
        ),
        experimental_design=ExperimentalDesign(
            biological_unit="mouse",
            batches=["b1", "b2"],
            total_samples=12,
            donor_replicates_per_condition={"AD": 6, "control": 6},
        )
    )
    assert manifest.study_id == "AD_test_001"
    assert manifest.biological_design.species == "mus_musculus"
    assert manifest.experimental_design.donor_replicates_per_condition["AD"] == 6

    # Test serialization
    data = manifest.model_dump()
    reconstructed = StudyManifest.model_validate(data)
    assert reconstructed.study_id == manifest.study_id


def test_task_contract_schema():
    contract = TaskContract(
        task_id="task_018",
        capability="trajectory_inference",
        input_artifacts=["adata://AD/microglia/v4"],
        allowed_operations=["build_neighbor_graph", "infer_trajectory"],
        forbidden_operations=["filter_cells", "normalize", "recluster"],
    )
    assert contract.task_id == "task_018"
    assert "recluster" in contract.forbidden_operations
    assert contract.retry_policy.max_execution_retry == 2


def test_artifact_metadata_schema():
    meta = ArtifactMetadata(
        artifact_id="adata://AD/microglia/v4",
        uri="adata://AD/microglia/v4",
        type=ArtifactType.ANNDATA,
        study_id="AD_001",
        created_by_task="task_017",
        operation="subset_cells",
        sha256_hash="sha256:abcd1234efgh5678",
        storage_path="/path/to/v4.h5ad",
    )
    assert meta.type == ArtifactType.ANNDATA
    assert meta.sha256_hash.startswith("sha256:")


def test_evidence_and_claim_schema():
    ev = EvidenceNode(
        evidence_id="E201",
        type=EvidenceType.PSEUDOBULK_DEG,
        score=0.92,
        summary="Apoe is upregulated in AD microglia.",
        source_task_id="task_008",
        source_artifact_uris=["table://AD/deg/v1"],
    )
    assert ev.evidence_id == "E201"

    claim = ClaimNode(
        claim_id="C101",
        statement="Apoe is significantly upregulated in AD microglia.",
        language_tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE,
        support_evidence_ids=["E201"],
        confidence=ConfidenceScore(association=0.92, mechanistic=0.85, causal=0.0, overall=0.89),
    )
    assert claim.claim_id == "C101"
    assert claim.confidence.overall == 0.89
