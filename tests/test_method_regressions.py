"""Regression tests for honest method provenance and statistical units."""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from scipy.sparse.linalg import svds as scipy_svds

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.clustering import ClusteringCapability
from eacbp.capabilities.deg import DifferentialExpressionCapability, DifferentialAbundanceCapability
from eacbp.auditor.statistical import StatisticalValidator
from eacbp.capabilities.integration import IntegrationCapability, compute_pca
from eacbp.capabilities.normalization import NormalizationCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.subset import SubsetCapability
from eacbp.capabilities.spatial.domain import build_spatial_neighborhood_graph
from eacbp.capabilities.trajectory import TrajectoryCapability
from eacbp.capabilities.perturbation.compound import CompoundPerturbationCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus, ScientificResultStatus


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
    output = reg.load_scdata(result.output_artifacts[0])
    assert "LEAKED_TRUTH" not in set(output.obs["cell_type"])
    assert "cluster" in output.obs and "X_embedding_2d" in output.obsm
    assert "leiden" not in output.obs and "X_umap" not in output.obsm


def test_aliases_return_honest_method_ids():
    assert IntegrationCapability("harmony").implementation_id == "batch_mean_centering_v1"
    assert TrajectoryCapability("paga_dpt").implementation_id == "root_distance_pseudotime_v1"


def test_integration_uses_the_explicit_batch_column(tmp_path, monkeypatch):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data()
    data.X = sparse.csr_matrix(data.X)
    data.obs["study_batch"] = ["site_a"] * 12 + ["site_b"] * 12
    captured = {}

    def capture_svds_start_vector(matrix, *, k, v0):
        captured["v0"] = v0.copy()
        return scipy_svds(matrix, k=k, v0=v0)

    import scipy.sparse.linalg
    monkeypatch.setattr(scipy.sparse.linalg, "svds", capture_svds_start_vector)
    uri = "adata://batch_mapping/raw/v1"
    _register(reg, uri, data)
    contract = TaskContract(
        task_id="integration_batch_mapping",
        capability="integration",
        method="batch_mean_centering_v1",
        input_artifacts=[uri],
        expected_outputs=["adata://batch_mapping/integrated/v3"],
        parameters={"batch_col": "study_batch", "n_components": 2, "random_seed": 7},
    )

    result = IntegrationCapability().execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert "batch_mean_centering_pca" in result.executed_operations
    expected_v0 = np.random.default_rng(7).standard_normal(min(data.X.shape))
    np.testing.assert_array_equal(captured["v0"], expected_v0)


def test_sparse_pca_is_stable_for_seed_after_global_rng_changes():
    matrix = sparse.diags(np.arange(1, 21, dtype=np.float32), shape=(30, 20), format="csr")
    global_rng_state = np.random.get_state()
    try:
        np.random.seed(1)
        first = compute_pca(matrix, n_components=4, random_seed=23)
        np.random.seed(987)
        second = compute_pca(matrix, n_components=4, random_seed=23)
    finally:
        np.random.set_state(global_rng_state)

    np.testing.assert_allclose(first, second, rtol=1e-6, atol=1e-6)


def test_subset_passes_contract_seed_into_sparse_pca(tmp_path, monkeypatch):
    import scipy.sparse.linalg

    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    data = _small_data()
    data.X = sparse.csr_matrix(data.X)
    data.obs["cell_type"] = "target"
    captured = {}

    def capture_svds_start_vector(matrix, *, k, v0):
        captured["v0"] = v0.copy()
        return scipy_svds(matrix, k=k, v0=v0)

    monkeypatch.setattr(scipy.sparse.linalg, "svds", capture_svds_start_vector)
    uri = "adata://subset_seed/raw/v1"
    _register(reg, uri, data)
    seed = 31
    contract = TaskContract(
        task_id="subset_seed",
        capability="subset_cells",
        method="subset_cells_v1",
        input_artifacts=[uri],
        expected_outputs=["adata://subset_seed/target_subset/v5"],
        parameters={"cell_type": "target", "random_seed": seed},
    )

    result = SubsetCapability().execute(contract, reg)

    assert result.status == TaskStatus.SUCCESS
    expected_v0 = np.random.default_rng(seed).standard_normal(min(data.X.shape))
    np.testing.assert_array_equal(captured["v0"], expected_v0)


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
    output = reg.load_scdata(result.output_artifacts[0])
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


def _baseline_science_case(tmp_path, capability, supported):
    rng = np.random.default_rng(304)
    first = rng.poisson(5, (40, 8))
    counts = np.vstack([first, first]).astype(np.int64)
    if supported:
        counts[:40, 0] += 100
    proportions = [8, 9, 7, 8, 2, 1, 3, 2] if supported else [8, 9, 7, 8] * 2
    obs = pd.DataFrame({"condition": ["treated"] * 40 + ["baseline"] * 40,
                        "donor": [f"d{i // 10}" for i in range(80)],
                        "state": [state for n in proportions for state in ["T"] * n + ["B"] * (10 - n)]},
                       index=[f"c{i}" for i in range(80)])
    data = SCData(counts, obs, pd.DataFrame({"gene_name": [f"g{i}" for i in range(8)]}), layers={"counts": counts})
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://method_test/raw/v1"
    _register(registry, uri, data)
    impl = DifferentialExpressionCapability() if capability == "deg" else DifferentialAbundanceCapability()
    task = TaskContract(task_id="baseline", capability=capability, method=impl.implementation_id,
                        input_artifacts=[uri], expected_outputs=["table://method_test/results/v1"],
                        parameters={"condition_a": "treated", "condition_b": "baseline", "state_col": "state"})
    return registry, task, impl.execute(task, registry)


@pytest.mark.parametrize("capability", ["deg", "differential_abundance"])
@pytest.mark.parametrize("supported", [False, True])
def test_baseline_science_recounts_supported_and_inconclusive_estimates(tmp_path, capability, supported):
    registry, task, result = _baseline_science_case(tmp_path, capability, supported)
    science = result.scientific_result
    assert science.status == (ScientificResultStatus.ESTIMATED_SUPPORTED if supported else ScientificResultStatus.ESTIMATED_INCONCLUSIVE)
    assert science.n_features_tested == (8 if capability == "deg" else 2)
    assert bool(science.n_features_supported) is supported
    report = StatisticalValidator().audit(task, result, registry)
    assert report.overall_passed, report
    assert next(check for check in report.checks if check.check_name == "baseline_scientific_result").passed


@pytest.mark.parametrize("capability", ["deg", "differential_abundance"])
@pytest.mark.parametrize("tamper", ["status", "tested", "supported", "missing"])
def test_baseline_science_rejects_forged_executor_completion(tmp_path, capability, tamper):
    registry, task, result = _baseline_science_case(tmp_path, capability, False)
    if tamper == "status":
        result.scientific_result.status = ScientificResultStatus.ESTIMATED_SUPPORTED
    elif tamper == "tested":
        result.scientific_result.n_features_tested += 10
    elif tamper == "supported":
        result.scientific_result.n_features_supported = 1
    else:
        result.scientific_result = None
    report = StatisticalValidator().audit(task, result, registry)
    check = next(check for check in report.checks if check.check_name == "baseline_scientific_result")
    assert not check.passed and not report.overall_passed


def test_baseline_descriptive_fallback_never_counts_as_an_estimate(tmp_path):
    data = _small_data(with_counts=True)
    data.obs = data.obs.drop(columns="donor")
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://method_test/raw/v1"
    _register(registry, uri, data)
    task = TaskContract(task_id="descriptive", capability="deg", method="donor_pseudobulk_welch_v1",
                        input_artifacts=[uri], expected_outputs=["table://method_test/results/v1"],
                        parameters={"condition_a": "treated", "condition_b": "baseline"})
    result = DifferentialExpressionCapability().execute(task, registry)
    assert result.scientific_result.status == ScientificResultStatus.NOT_ESTIMABLE
    assert result.scientific_result.n_features_tested == result.scientific_result.n_features_supported == 0
    report = StatisticalValidator().audit(task, result, registry)
    assert next(check for check in report.checks if check.check_name == "baseline_scientific_result").passed
    result.scientific_result.status = ScientificResultStatus.ESTIMATED_INCONCLUSIVE
    result.scientific_result.n_features_tested = 8
    report = StatisticalValidator().audit(task, result, registry)
    assert not next(check for check in report.checks if check.check_name == "baseline_scientific_result").passed


def test_large_spatial_graph_uses_sparse_knn_storage():
    coords = np.random.default_rng(9).random((1601, 2))
    W, distances, W_norm = build_spatial_neighborhood_graph(coords, k_neighbors=4)
    assert sparse.issparse(W)
    assert sparse.issparse(distances)
    assert sparse.issparse(W_norm)
    assert W.nnz <= 2 * 1601 * 4
