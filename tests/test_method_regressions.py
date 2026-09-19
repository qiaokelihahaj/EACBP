"""Regression tests for honest method provenance and statistical units."""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.clustering import ClusteringCapability
from eacbp.capabilities.deg import DifferentialExpressionCapability
from eacbp.capabilities.integration import IntegrationCapability
from eacbp.capabilities.normalization import NormalizationCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial.domain import build_spatial_neighborhood_graph
from eacbp.capabilities.trajectory import TrajectoryCapability
from eacbp.capabilities.perturbation.compound import CompoundPerturbationCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _register(registry, uri, data, artifact_type=ArtifactType.ANNDATA):
    registry.register(
        uri_str=uri,
        payload=data.to_dict(),
        artifact_type=artifact_type,
        study_id=uri.split("/")[2],
        created_by_task="test_input",
        operation="input",
    )


def _small_data(with_counts=False):
    rng = np.random.default_rng(3)
    X = rng.poisson(2, size=(24, 8)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cell_id": [f"c{i}" for i in range(24)],
            "condition": ["treated"] * 12 + ["baseline"] * 12,
            "donor": [f"t{i}" for i in range(4) for _ in range(3)] + [f"b{i}" for i in range(4) for _ in range(3)],
            "cell_type_ground_truth": ["LEAKED_TRUTH"] * 24,
        }
    )
    var = pd.DataFrame({"gene_name": ["Cx3cr1", "P2ry12", "Tmem119", "Apoe", "Trem2", "G1", "G2", "G3"]})
    return SCData(X=X, obs=obs, var=var, layers={"counts": X.copy()} if with_counts else {})


def test_clustering_does_not_use_ground_truth_and_reports_kmeans(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data()
    uri = "adata://method_test/raw/v1"
    _register(reg, uri, data)
    contract = TaskContract(
        task_id="cluster_test",
        capability="clustering",
        method="kmeans_marker_embedding_v1",
        input_artifacts=[uri],
        expected_outputs=["adata://method_test/annotated/v4"],
        parameters={"k_clusters": 2, "random_seed": 7},
    )
    result = ClusteringCapability().execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert result.method_used == "kmeans_marker_embedding_v1"
    _, payload = reg.get(result.output_artifacts[0])
    output = SCData.from_dict(payload)
    assert "LEAKED_TRUTH" not in set(output.obs["cell_type"])
    assert "cluster" in output.obs and "X_embedding_2d" in output.obsm
    assert "leiden" not in output.obs and "X_umap" not in output.obsm


def test_aliases_return_honest_method_ids():
    assert IntegrationCapability("harmony").implementation_id == "batch_mean_centering_v1"
    assert TrajectoryCapability("paga_dpt").implementation_id == "root_distance_pseudotime_v1"


def test_unknown_compound_cannot_invent_reversal_signature(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    uri = "adata://method_test/raw/v1"
    _register(reg, uri, _small_data())
    contract = TaskContract(task_id="unknown_drug", capability="compound_perturbation_simulation",
                            method="in_silico_compound_response_v1", input_artifacts=[uri],
                            parameters={"compound_name": "unmeasured_compound", "condition_a": "treated", "condition_b": "baseline"})
    with pytest.raises(ValueError, match="explicit drug_signature"):
        CompoundPerturbationCapability().execute(contract, reg)


def test_normalization_preserves_counts_layer(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data(with_counts=True)
    uri = "adata://method_test/raw/v1"
    _register(reg, uri, data)
    contract = TaskContract(
        task_id="norm_test",
        capability="normalization",
        method="library_size_log1p_v1",
        input_artifacts=[uri],
        expected_outputs=["adata://method_test/normalized/v2"],
    )
    result = NormalizationCapability().execute(contract, reg)
    _, payload = reg.get(result.output_artifacts[0])
    output = SCData.from_dict(payload)
    assert "counts" in output.layers
    assert np.array_equal(output.layers["counts"], data.layers["counts"])
    assert output.uns["normalization"]["counts_layer"] == "counts"


def test_trajectory_constant_input_returns_schema_stable_empty_table(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data()
    data.X[:] = 1.0
    uri = "adata://method_test/raw/v1"
    _register(reg, uri, data)
    contract = TaskContract(
        task_id="trajectory_test",
        capability="trajectory_inference",
        method="root_distance_pseudotime_v1",
        input_artifacts=[uri],
        expected_outputs=["table://method_test/trajectory_results/v1"],
        parameters={"random_seed": 11},
    )
    result = TrajectoryCapability().execute(contract, reg)
    assert result.method_used == "root_distance_pseudotime_v1"
    assert result.metrics["stability_score"] != 0.85
    _, table = reg.get(result.output_artifacts[0])
    assert list(table.columns) == ["gene", "spearman_rho", "p_value", "fdr_q_value", "trend"]


def test_deg_uses_counts_donor_pseudobulk_and_observed_condition_names(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data(with_counts=True)
    uri = "adata://method_test/raw/v1"
    _register(reg, uri, data)
    contract = TaskContract(
        task_id="deg_test",
        capability="deg",
        method="donor_pseudobulk_welch_v1",
        input_artifacts=[uri],
        expected_outputs=["table://method_test/deg_results/v1"],
        parameters={"condition_a": "treated", "condition_b": "baseline"},
    )
    result = DifferentialExpressionCapability().execute(contract, reg)
    assert result.method_used == "donor_pseudobulk_welch_v1"
    _, table = reg.get(result.output_artifacts[0])
    assert set(table["statistical_unit"]) == {"donor_pseudobulk"}
    assert set(table["condition_a"]) == {"treated"}
    assert set(table["condition_b"]) == {"baseline"}
    assert table["effect_definition"].str.contains("donor").all()


def test_large_spatial_graph_uses_sparse_knn_storage():
    coords = np.random.default_rng(9).random((1601, 2))
    W, distances, W_norm = build_spatial_neighborhood_graph(coords, k_neighbors=4)
    assert sparse.issparse(W)
    assert sparse.issparse(distances)
    assert sparse.issparse(W_norm)
    assert W.nnz <= 2 * 1601 * 4
