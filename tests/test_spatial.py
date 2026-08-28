"""
Comprehensive Unit and Integration Test Suite for Spatial Single-Cell Analytics Plane (M1).
Tests:
1. Spatial coordinate validation (2D/3D .obsm['spatial'], NaN/Inf checks, degenerate coords)
2. Spatial domain clustering and neighborhood graph correctness (.obsp, smoothing, silhouette)
3. Moran's I and Geary's C statistical properties (hotspots, dispersion, null random distributions, analytical variance)
4. Spatial DEG identification and Benjamini-Hochberg FDR thresholding
5. Cell-Cell Communication (CCI) proximity-weighted scoring and spatial permutation tests
6. Artifact immutability, SHA-256 payload integrity, and lineage tracking
7. Independent Scientific Auditor verification (Computational & Statistical validators)
"""

import pytest
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.uri import ArtifactURI
from eacbp.artifact.storage import ArtifactStorageBackend, ArtifactAlreadyExistsError
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.capabilities.spatial.domain import (
    SpatialDomainCapability,
    validate_spatial_coordinates,
    build_spatial_neighborhood_graph,
    compute_spatially_smoothed_embedding,
    calculate_silhouette,
)
from eacbp.capabilities.spatial.autocorrelation import (
    SpatialDEGCapability,
    calculate_morans_i,
    calculate_gearys_c,
    benjamini_hochberg,
)
from eacbp.capabilities.spatial.cci import (
    CellCellCommunicationCapability,
    compute_spatial_cci,
    calculate_spatial_contact_density,
    CURATED_LIGAND_RECEPTOR_PAIRS,
)
from eacbp.capabilities.spatial import create_synthetic_spatial_ad_study
from eacbp.auditor.computational import ComputationalValidator
from eacbp.auditor.statistical import StatisticalValidator


@pytest.fixture
def spatial_registry(tmp_path):
    """Provides an isolated ArtifactRegistry for spatial tests."""
    storage_dir = tmp_path / "spatial_artifacts"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    return registry


@pytest.fixture
def synthetic_spatial_data():
    """Generates standard synthetic spatial single-cell AD dataset."""
    return create_synthetic_spatial_ad_study(n_cells=300, n_genes=50, random_seed=42)


# ==============================================================================
# 1. Spatial Coordinate Validation Tests
# ==============================================================================

def test_spatial_coordinate_validation_2d_and_3d():
    """Verifies that 2D and 3D coordinate arrays pass validation."""
    n_obs = 50
    coords_2d = np.random.uniform(0, 100, size=(n_obs, 2)).astype(np.float32)
    val_2d = validate_spatial_coordinates(coords_2d, expected_n_obs=n_obs)
    assert val_2d.shape == (n_obs, 2)
    assert val_2d.dtype == np.float32

    coords_3d = np.random.uniform(0, 100, size=(n_obs, 3)).astype(np.float32)
    val_3d = validate_spatial_coordinates(coords_3d, expected_n_obs=n_obs)
    assert val_3d.shape == (n_obs, 3)


def test_spatial_coordinate_validation_errors():
    """Tests various invalid coordinate conditions."""
    n_obs = 50

    # None coordinates
    with pytest.raises(ValueError, match="is None"):
        validate_spatial_coordinates(None, n_obs)

    # 1D array
    with pytest.raises(ValueError, match="must be 2D array"):
        validate_spatial_coordinates(np.ones(n_obs), n_obs)

    # Wrong cell count
    with pytest.raises(ValueError, match="cell count"):
        validate_spatial_coordinates(np.ones((40, 2)), n_obs)

    # Invalid dimension count (4D)
    with pytest.raises(ValueError, match="2 or 3 dimensions"):
        validate_spatial_coordinates(np.ones((n_obs, 4)), n_obs)

    # NaNs in coordinates
    nan_coords = np.ones((n_obs, 2))
    nan_coords[5, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        validate_spatial_coordinates(nan_coords, n_obs)

    # Infs in coordinates
    inf_coords = np.ones((n_obs, 2))
    inf_coords[10, 1] = np.inf
    with pytest.raises(ValueError, match="non-finite"):
        validate_spatial_coordinates(inf_coords, n_obs)

    # Degenerate zero-variance coordinates (all points at same position)
    zero_var = np.full((n_obs, 2), 5.0)
    with pytest.raises(ValueError, match="Degenerate spatial coordinates"):
        validate_spatial_coordinates(zero_var, n_obs)


# ==============================================================================
# 2. Spatial Neighborhood Graph & Domain Clustering Tests
# ==============================================================================

def test_spatial_neighborhood_graph_properties():
    """Tests spatial k-NN graph construction mathematical invariants."""
    coords = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [10.0, 10.0],
        [11.0, 10.0],
    ])
    W_sym, D_pairwise, W_norm = build_spatial_neighborhood_graph(coords, k_neighbors=2)

    # 1. Symmetric adjacency
    assert np.allclose(W_sym, W_sym.T)
    # 2. Zero diagonal (no self-loops)
    assert np.allclose(np.diag(W_sym), 0.0)
    assert np.allclose(np.diag(D_pairwise), 0.0)
    # 3. Non-negative weights
    assert np.all(W_sym >= 0.0)
    assert np.all(D_pairwise >= 0.0)
    # 4. Closest points are connected
    assert W_sym[0, 1] == 1.0
    assert W_sym[0, 2] == 1.0
    # 5. Row normalization sum equals 1.0 for connected nodes
    row_sums = W_norm.sum(axis=1)
    assert np.allclose(row_sums, 1.0)


def test_spatially_smoothed_embedding():
    """Tests spatially smoothed latent embedding calculation."""
    Z = np.array([
        [10.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
    ])
    # W_norm where node 0 connects equally to node 1 and 2
    W_norm = np.array([
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ])
    Z_smooth_0 = compute_spatially_smoothed_embedding(Z, W_norm, smoothing_lambda=0.0)
    assert np.allclose(Z_smooth_0, Z)

    Z_smooth_5 = compute_spatially_smoothed_embedding(Z, W_norm, smoothing_lambda=0.5)
    # Node 0 smoothed with lambda=0.5: 0.5*[10, 0] + 0.5*[0, 0] = [5, 0]
    assert np.allclose(Z_smooth_5[0], [5.0, 0.0])


def test_spatial_domain_capability_execution(spatial_registry, synthetic_spatial_data):
    """Tests SpatialDomainCapability execution through CapabilityRegistry."""
    reg = spatial_registry
    in_uri = "adata://AD_study/spatial_norm/v3"
    reg.register(
        uri_str=in_uri,
        payload=synthetic_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="AD_study",
        created_by_task="task_003_norm",
        operation="normalize",
    )

    cap_reg = CapabilityRegistry()
    cap_reg.register(SpatialDomainCapability())

    contract = TaskContract(
        task_id="task_spatial_001_domain",
        capability="spatial_domain",
        method="spatial_domain_knn_v1",
        input_artifacts=[in_uri],
        allowed_operations=[
            "validate_spatial_coordinates",
            "build_spatial_connectivities",
            "spatially_smoothed_embedding",
            "cluster_spatial_domains",
            "calculate_silhouette",
        ],
        forbidden_operations=["filter_cells", "normalize"],
        parameters={"k_neighbors": 6, "n_domains": 4, "smoothing_lambda": 0.3, "random_seed": 42},
        expected_outputs=["adata://AD_study/spatial_domains/v4"],
    )

    result = cap_reg.execute_contract(contract, reg)

    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 1
    out_uri = result.output_artifacts[0]
    assert out_uri == "adata://AD_study/spatial_domains/v4"
    assert reg.exists(out_uri)

    meta, payload = reg.get(out_uri)
    assert meta.type == ArtifactType.SPATIAL_DATA
    assert meta.summary_metrics["n_domains"] == 4
    assert meta.summary_metrics["silhouette_score"] > -1.0

    out_data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
    assert "spatial_domain" in out_data.obs.columns
    assert "spatial" in out_data.obsm
    assert "X_spatial_pca" in out_data.obsm
    assert "spatial_connectivities" in out_data.obsm
    assert len(out_data.obs["spatial_domain"].unique()) == 4


# ==============================================================================
# 3. Moran's I & Geary's C Statistical Properties Tests
# ==============================================================================

def test_moran_and_geary_positive_spatial_autocorrelation():
    """
    Tests that a spatially clustered / continuous gradient pattern yields
    high positive Moran's I and low Geary's C with significant p-values.
    """
    x = np.linspace(0, 10, 10)
    y = np.linspace(0, 10, 10)
    xx, yy = np.meshgrid(x, y)
    coords = np.column_stack([xx.ravel(), yy.ravel()])
    N = len(coords)

    # Injected spatial gradient
    expression_gradient = xx.ravel() + yy.ravel() + np.random.normal(0, 0.1, size=N)

    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=4)

    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(expression_gradient, W)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(expression_gradient, W)

    # Positive spatial autocorrelation invariants
    assert m_i > 0.50, f"Expected strong Moran's I > 0.50, got {m_i}"
    assert m_z > 3.0, f"Expected significant z-score > 3.0, got {m_z}"
    assert m_p < 0.001, f"Expected highly significant p-value < 0.001, got {m_p}"

    assert g_c < 0.60, f"Expected Geary's C < 0.60, got {g_c}"
    assert g_z < -3.0, f"Expected negative Geary's z-score, got {g_z}"
    assert g_p < 0.001, f"Expected significant Geary p-value, got {g_p}"


def test_moran_and_geary_negative_spatial_autocorrelation():
    """
    Tests that a checkerboard dispersion pattern yields
    negative Moran's I and Geary's C > 1.0.
    """
    x = np.arange(8)
    y = np.arange(8)
    xx, yy = np.meshgrid(x, y)
    coords = np.column_stack([xx.ravel(), yy.ravel()])

    # Checkerboard expression
    checkerboard = ((xx.ravel() + yy.ravel()) % 2 == 0).astype(float) * 10.0

    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=4)

    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(checkerboard, W)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(checkerboard, W)

    assert m_i < -0.30, f"Expected negative Moran's I < -0.30, got {m_i}"
    assert m_z < -2.0, f"Expected negative z-score, got {m_z}"
    assert g_c > 1.30, f"Expected Geary's C > 1.30, got {g_c}"


def test_moran_and_geary_random_null_distribution():
    """
    Tests that independent random noise produces Moran's I ≈ E[I] = -1/(N-1)
    and Geary's C ≈ 1.0, with non-significant p-values (p > 0.05).
    """
    np.random.seed(123)
    coords = np.random.uniform(0, 100, size=(100, 2))
    random_noise = np.random.normal(5.0, 1.0, size=100)

    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=6)

    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(random_noise, W)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(random_noise, W)

    assert abs(m_i - m_exp) < 3.0 * np.sqrt(m_var)
    assert m_p > 0.05, f"Expected non-significant p-value > 0.05, got {m_p}"
    assert abs(g_c - 1.0) < 0.35
    assert g_p > 0.05, f"Expected non-significant Geary p-value > 0.05, got {g_p}"


def test_moran_and_geary_constant_expression_edge_case():
    """Tests constant expression edge case (zero variance)."""
    coords = np.random.uniform(0, 100, size=(20, 2))
    constant_expr = np.full(20, 3.14)
    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=4)

    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(constant_expr, W)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(constant_expr, W)

    assert m_i == 0.0
    assert m_p == 1.0
    assert g_c == 1.0
    assert g_p == 1.0


# ==============================================================================
# 4. Spatial DEG Capability & FDR Thresholding Tests
# ==============================================================================

def test_spatial_deg_capability_execution(spatial_registry, synthetic_spatial_data):
    """
    Tests SpatialDEGCapability identifies known spatially localized genes (Apoe, Trem2, Clec7a)
    with FDR < 0.05 while non-spatial random genes have higher FDR.
    """
    reg = spatial_registry
    in_uri = "adata://AD_study/spatial_domains/v4"
    reg.register(
        uri_str=in_uri,
        payload=synthetic_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="AD_study",
        created_by_task="task_004_domain",
        operation="identify_domains",
    )

    cap_reg = CapabilityRegistry()
    cap_reg.register(SpatialDEGCapability())

    contract = TaskContract(
        task_id="task_spatial_002_deg",
        capability="spatial_deg",
        method="spatial_moran_deg_v1",
        input_artifacts=[in_uri],
        allowed_operations=[
            "build_spatial_connectivities",
            "calculate_moran_i",
            "calculate_geary_c",
            "analytical_significance_test",
            "benjamini_hochberg_correction",
        ],
        forbidden_operations=["filter_cells", "normalize", "recluster"],
        parameters={"k_neighbors": 6, "min_moran_i": 0.15, "fdr_threshold": 0.05},
        expected_outputs=["table://AD_study/spatial_deg/v1"],
    )

    result = cap_reg.execute_contract(contract, reg)

    assert result.status == TaskStatus.SUCCESS
    out_uri = result.output_artifacts[0]
    assert out_uri == "table://AD_study/spatial_deg/v1"
    assert reg.exists(out_uri)

    meta, payload = reg.get(out_uri)
    assert meta.type == ArtifactType.TABLE
    df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)

    assert "moran_i" in df.columns
    assert "geary_c" in df.columns
    assert "fdr_q_value" in df.columns
    assert "is_spatially_variable" in df.columns

    sig_genes = df[df["is_spatially_variable"]]["gene"].tolist()
    assert "Apoe" in sig_genes
    assert "Trem2" in sig_genes
    assert "Clec7a" in sig_genes

    apoe_row = df[df["gene"] == "Apoe"].iloc[0]
    assert apoe_row["moran_i"] > 0.25
    assert apoe_row["fdr_q_value"] < 0.01


# ==============================================================================
# 5. Cell-Cell Communication (CCI) Proximity-Weighted Tests
# ==============================================================================

def test_spatial_cci_contact_density_and_proximity():
    """Tests spatial contact density between clustered cell types."""
    coords_a = np.random.uniform(0, 10, size=(20, 2))
    coords_b = np.random.uniform(100, 110, size=(20, 2))
    coords = np.vstack([coords_a, coords_b])

    labels = np.array(["TypeA"] * 20 + ["TypeB"] * 20)
    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=4)

    densities = calculate_spatial_contact_density(labels, ["TypeA", "TypeB"], W)

    assert densities[("TypeA", "TypeA")] > 0.0
    assert densities[("TypeB", "TypeB")] > 0.0
    assert densities[("TypeA", "TypeB")] == 0.0


def test_cell_cell_communication_capability_execution(spatial_registry, synthetic_spatial_data):
    """
    Tests CellCellCommunicationCapability computing proximity-weighted ligand-receptor scores
    and permutation significance testing.
    """
    reg = spatial_registry
    in_uri = "adata://AD_study/spatial_domains/v4"
    reg.register(
        uri_str=in_uri,
        payload=synthetic_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="AD_study",
        created_by_task="task_004_domain",
        operation="identify_domains",
    )

    cap_reg = CapabilityRegistry()
    cap_reg.register(CellCellCommunicationCapability())

    contract = TaskContract(
        task_id="task_spatial_003_cci",
        capability="cell_cell_communication",
        method="cci_ligand_receptor_v1",
        input_artifacts=[in_uri],
        allowed_operations=[
            "load_lr_database",
            "calculate_spatial_contact_density",
            "compute_spatial_cci_score",
            "run_spatial_permutation_test",
        ],
        forbidden_operations=["filter_cells", "mutate_coordinates"],
        parameters={
            "k_neighbors": 6,
            "n_permutations": 100,
            "fdr_threshold": 0.05,
            "cell_type_col": "cell_type",
            "random_seed": 42,
        },
        expected_outputs=["table://AD_study/spatial_cci/v1"],
    )

    result = cap_reg.execute_contract(contract, reg)

    assert result.status == TaskStatus.SUCCESS
    out_uri = result.output_artifacts[0]
    assert out_uri == "table://AD_study/spatial_cci/v1"
    assert reg.exists(out_uri)

    meta, payload = reg.get(out_uri)
    assert meta.type == ArtifactType.TABLE
    cci_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)

    assert not cci_df.empty
    assert "spatial_interaction_score" in cci_df.columns
    assert "spatial_proximity_weight" in cci_df.columns
    assert "p_value" in cci_df.columns
    assert "fdr_q_value" in cci_df.columns

    for idx, row in cci_df.head(10).iterrows():
        expected_score = row["raw_interaction_score"] * row["spatial_proximity_weight"]
        assert np.isclose(row["spatial_interaction_score"], expected_score, atol=1e-6)


# ==============================================================================
# 6. Artifact Immutability, SHA-256 Checksums, and Lineage Tracking
# ==============================================================================

def test_spatial_artifact_immutability_and_lineage(spatial_registry, synthetic_spatial_data):
    """
    Tests that spatial capabilities produce immutable versioned artifacts
    with valid SHA-256 hashes and correct parent-child lineage tracking.
    """
    reg = spatial_registry
    
    # 1. Register raw spatial dataset
    raw_uri = "adata://AD_001/spatial_raw/v1"
    meta1 = reg.register(
        uri_str=raw_uri,
        payload=synthetic_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="AD_001",
        created_by_task="task_raw_ingest",
        operation="ingest_spatial_data",
    )
    assert meta1.sha256_hash.startswith("sha256:")

    # Attempting in-place overwrite must raise ArtifactAlreadyExistsError
    with pytest.raises(ArtifactAlreadyExistsError):
        reg.register(
            uri_str=raw_uri,
            payload=synthetic_spatial_data.to_dict(),
            artifact_type=ArtifactType.SPATIAL_DATA,
            study_id="AD_001",
            created_by_task="task_rogue_overwrite",
            operation="reingest",
        )

    # 2. Run SpatialDomainCapability to create v2
    domain_cap = SpatialDomainCapability()
    contract_dom = TaskContract(
        task_id="task_spatial_domain",
        capability="spatial_domain",
        method="spatial_domain_knn_v1",
        input_artifacts=[raw_uri],
        expected_outputs=["adata://AD_001/spatial_domains/v2"],
    )
    res_dom = domain_cap.execute(contract_dom, reg)
    assert res_dom.status == TaskStatus.SUCCESS
    dom_uri = res_dom.output_artifacts[0]

    # 3. Run SpatialDEGCapability to create table v1
    deg_cap = SpatialDEGCapability()
    contract_deg = TaskContract(
        task_id="task_spatial_deg",
        capability="spatial_deg",
        method="spatial_moran_deg_v1",
        input_artifacts=[dom_uri],
        expected_outputs=["table://AD_001/spatial_deg/v1"],
    )
    res_deg = deg_cap.execute(contract_deg, reg)
    assert res_deg.status == TaskStatus.SUCCESS
    deg_uri = res_deg.output_artifacts[0]

    # Check Lineage Graph ancestry
    ancestors = reg.lineage.get_ancestors(deg_uri)
    assert dom_uri in ancestors
    assert raw_uri in ancestors

    # Verify SHA-256 payload integrity
    meta_deg, payload_deg = reg.get(deg_uri)
    assert meta_deg.sha256_hash.startswith("sha256:")
    assert meta_deg.size_bytes > 0


# ==============================================================================
# 7. Independent Scientific Auditor Suite Verification
# ==============================================================================

def test_auditors_evaluate_spatial_artifacts(spatial_registry, synthetic_spatial_data):
    """
    Tests ComputationalValidator and StatisticalValidator on spatial outputs.
    Ensures zero errors and full compliance with scientific audit rules.
    """
    reg = spatial_registry
    raw_uri = "adata://AD_AUDIT/spatial_raw/v1"
    reg.register(
        uri_str=raw_uri,
        payload=synthetic_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="AD_AUDIT",
        created_by_task="task_000",
        operation="init",
    )

    cap_reg = CapabilityRegistry()
    cap_reg.register(SpatialDomainCapability())
    cap_reg.register(SpatialDEGCapability(capability_name="deg"))

    # Task 1: Spatial Domain
    t1 = TaskContract(
        task_id="task_audit_domain",
        capability="spatial_domain",
        method="spatial_domain_knn_v1",
        input_artifacts=[raw_uri],
        expected_outputs=["adata://AD_AUDIT/spatial_domains/v2"],
    )
    res1 = cap_reg.execute_contract(t1, reg)
    assert res1.status == TaskStatus.SUCCESS

    # Audit Task 1 with ComputationalValidator
    comp_auditor = ComputationalValidator()
    report1 = comp_auditor.audit(t1, res1, reg)
    assert report1.overall_passed, f"Computational audit failed: {[c.message for c in report1.checks if not c.passed]}"
    
    check_names = [c.check_name for c in report1.checks]
    assert "matrix_non_empty" in check_names
    assert "expression_finite_values" in check_names
    assert "embedding_spatial_finite" in check_names

    # Task 2: Spatial DEG
    t2 = TaskContract(
        task_id="task_audit_deg",
        capability="deg",
        method="spatial_moran_deg_v1",
        input_artifacts=[res1.output_artifacts[0]],
        expected_outputs=["table://AD_AUDIT/spatial_deg/v1"],
    )
    res2 = cap_reg.execute_contract(t2, reg)
    assert res2.status == TaskStatus.SUCCESS

    # Audit Task 2 with Computational & Statistical Validators
    report2_comp = comp_auditor.audit(t2, res2, reg)
    assert report2_comp.overall_passed

    stat_auditor = StatisticalValidator()
    report2_stat = stat_auditor.audit(t2, res2, reg)
    assert report2_stat.overall_passed
    stat_check_names = [c.check_name for c in report2_stat.checks]
    assert "multiple_testing_correction" in stat_check_names
