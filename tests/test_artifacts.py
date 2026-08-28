"""
Unit tests for Artifact storage, immutability guarantees, and lineage DAG.
"""

import pytest
import os
import shutil
import pandas as pd
import numpy as np

from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.uri import ArtifactURI
from eacbp.artifact.storage import ArtifactStorageBackend, ArtifactAlreadyExistsError
from eacbp.artifact.lineage import LineageGraph
from eacbp.artifact.registry import ArtifactRegistry


@pytest.fixture
def temp_storage(tmp_path):
    storage_dir = tmp_path / "artifacts_test"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    yield registry
    if storage_dir.exists():
        shutil.rmtree(storage_dir)


def test_artifact_uri_parsing():
    uri_str = "adata://AD_mouse_001/microglia_subset/v4"
    uri = ArtifactURI.parse(uri_str)
    assert uri.scheme == "adata"
    assert uri.study_id == "AD_mouse_001"
    assert uri.name == "microglia_subset"
    assert uri.version == "v4"

    next_uri = uri.next_version()
    assert next_uri.to_string() == "adata://AD_mouse_001/microglia_subset/v5"

    branch_uri = uri.branch("harmony")
    assert branch_uri.to_string() == "adata://AD_mouse_001/microglia_subset/v4_harmony"


def test_artifact_immutability_and_overwrites(temp_storage):
    registry = temp_storage
    df = pd.DataFrame({"gene": ["Apoe", "Trem2"], "log2fc": [2.1, 1.8]})

    # Register initial v1
    meta1 = registry.register(
        uri_str="table://AD_001/deg_results/v1",
        payload=df,
        artifact_type=ArtifactType.TABLE,
        study_id="AD_001",
        created_by_task="task_001",
        operation="run_deg",
    )
    assert meta1.sha256_hash.startswith("sha256:")
    assert registry.exists("table://AD_001/deg_results/v1")

    # Invariant 2 test: Attempting to overwrite existing v1 must raise ArtifactAlreadyExistsError
    with pytest.raises(ArtifactAlreadyExistsError):
        registry.register(
            uri_str="table://AD_001/deg_results/v1",
            payload=df,
            artifact_type=ArtifactType.TABLE,
            study_id="AD_001",
            created_by_task="task_002",
            operation="re_run_deg",
        )

    # Saving as v2 succeeds
    meta2 = registry.register(
        uri_str="table://AD_001/deg_results/v2",
        payload=df,
        artifact_type=ArtifactType.TABLE,
        study_id="AD_001",
        created_by_task="task_002",
        operation="re_run_deg",
        parent_uris=["table://AD_001/deg_results/v1"],
    )
    assert meta2.uri == "table://AD_001/deg_results/v2"


def test_lineage_graph_and_branch_comparison(temp_storage):
    registry = temp_storage
    
    # Simulate branching: v3 -> (v4a_harmony, v4b_scvi)
    v3_meta = registry.register(
        uri_str="adata://AD/normalized/v3",
        payload={"X": np.ones((10, 10))},
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_003",
        operation="normalize",
    )

    v4a_meta = registry.register(
        uri_str="adata://AD/integrated/v4a_harmony",
        payload={"X": np.ones((10, 10)) * 2},
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_004a",
        operation="harmony_integration",
        parent_uris=["adata://AD/normalized/v3"],
        parameters={"method": "harmony"},
    )

    v4b_meta = registry.register(
        uri_str="adata://AD/integrated/v4b_scvi",
        payload={"X": np.ones((10, 10)) * 3},
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_004b",
        operation="scvi_integration",
        parent_uris=["adata://AD/normalized/v3"],
        parameters={"method": "scvi"},
    )

    # Check ancestry
    ancestors = registry.lineage.get_ancestors("adata://AD/integrated/v4a_harmony")
    assert "adata://AD/normalized/v3" in ancestors

    # Compare branches
    comparison = registry.lineage.compare_branches(
        "adata://AD/integrated/v4a_harmony",
        "adata://AD/integrated/v4b_scvi"
    )
    assert "adata://AD/normalized/v3" in comparison["common_ancestors"]
    assert comparison["operation_a"] == "harmony_integration"
    assert comparison["operation_b"] == "scvi_integration"
