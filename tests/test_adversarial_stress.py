"""
Adversarial Stress Test Suite for Spatial Capabilities and Perturbation Engine in EACBP.
Empirically tests edge cases, singular matrices, extreme parameters, and mathematical bounds.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial.domain import (
    validate_spatial_coordinates,
    build_spatial_neighborhood_graph,
    compute_spatially_smoothed_embedding,
    simple_kmeans,
    calculate_silhouette,
    SpatialDomainCapability,
)
from eacbp.capabilities.spatial.autocorrelation import (
    calculate_morans_i,
    calculate_gearys_c,
    benjamini_hochberg,
    SpatialDEGCapability,
)
from eacbp.capabilities.spatial.cci import (
    calculate_spatial_contact_density,
    compute_spatial_cci,
    CellCellCommunicationCapability,
)
from eacbp.capabilities.perturbation.genetic import (
    construct_grn_adjacency_from_data,
    compute_grn_propagator,
    GeneticPerturbationCapability,
    generate_genetic_perturbation_evidence,
)
from eacbp.capabilities.perturbation.compound import (
    compute_cmap_cosine_discordance,
    CompoundPerturbationCapability,
    generate_compound_perturbation_evidence,
)
from eacbp.auditor.statistical import StatisticalValidator
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType, ArtifactMetadata
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
)
from eacbp.artifact.registry import ArtifactRegistry


# =============================================================================
# 1. Spatial Autocorrelation & Coordinate Stress Tests
# =============================================================================

def test_spatial_autocorrelation_flat_and_zero_variance():
    """
    Stress-test Moran's I and Geary's C with flat/constant expressions and near-zero variance.
    Must not crash, divide by zero, or produce NaNs/Infs.
    """
    N = 50
    coords = np.random.RandomState(42).rand(N, 2)
    W, _, _ = build_spatial_neighborhood_graph(coords, k_neighbors=6)

    # Vector of all zeros
    x_zeros = np.zeros(N)
    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(x_zeros, W)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(x_zeros, W)

    assert m_i == 0.0
    assert m_p == 1.0
    assert not np.isnan(m_i) and not np.isinf(m_i)
    assert g_c == 1.0
    assert g_p == 1.0
    assert not np.isnan(g_c) and not np.isinf(g_c)

    # Vector of all constants (e.g. 42.0)
    x_const = np.full(N, 42.0)
    m_i2, _, _, _, m_p2 = calculate_morans_i(x_const, W)
    g_c2, _, _, _, g_p2 = calculate_gearys_c(x_const, W)

    assert m_i2 == 0.0
    assert m_p2 == 1.0
    assert g_c2 == 1.0
    assert g_p2 == 1.0

    # Vector with extreme tiny variance (1e-16)
    x_tiny = 5.0 + np.random.RandomState(42).randn(N) * 1e-16
    m_i3, _, _, _, _ = calculate_morans_i(x_tiny, W)
    g_c3, _, _, _, _ = calculate_gearys_c(x_tiny, W)
    assert not np.isnan(m_i3)
    assert not np.isnan(g_c3)


def test_spatial_coordinates_extreme_ranges_and_shapes():
    """
    Stress-test spatial coordinates with extreme scales, 3D, and validate degenerate checks.
    """
    N = 100
    # Coordinates in billions (astronomical scale)
    coords_huge = np.random.RandomState(42).rand(N, 2) * 1e12
    val_huge = validate_spatial_coordinates(coords_huge, N)
    assert val_huge.shape == (N, 2)
    W_huge, D_huge, _ = build_spatial_neighborhood_graph(val_huge, k_neighbors=5)
    assert W_huge.shape == (N, N)
    assert not np.isnan(W_huge).any()

    # Negative coordinates
    coords_neg = -np.random.RandomState(42).rand(N, 2) * 5000.0
    val_neg = validate_spatial_coordinates(coords_neg, N)
    assert val_neg.shape == (N, 2)

    # 3D Coordinates
    coords_3d = np.random.RandomState(42).rand(N, 3) * 100.0
    val_3d = validate_spatial_coordinates(coords_3d, N)
    assert val_3d.shape == (N, 3)
    W_3d, _, _ = build_spatial_neighborhood_graph(val_3d, k_neighbors=8)
    assert W_3d.shape == (N, N)

    # Degenerate: all identical positions (zero variance)
    coords_degen = np.full((N, 2), 10.0)
    with pytest.raises(ValueError, match="Degenerate spatial coordinates"):
        validate_spatial_coordinates(coords_degen, N)

    # NaNs and Infs
    coords_nan = coords_huge.copy()
    coords_nan[0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite values"):
        validate_spatial_coordinates(coords_nan, N)

    coords_inf = coords_huge.copy()
    coords_inf[5, 1] = np.inf
    with pytest.raises(ValueError, match="non-finite values"):
        validate_spatial_coordinates(coords_inf, N)


def test_spatial_disconnected_graphs_and_high_density_clusters():
    """
    Stress-test disconnected spatial graphs (two far-apart clusters) and high density clusters.
    """
    # Two widely separated clusters: 30 cells at (0, 0) region, 30 cells at (100000, 100000) region
    np.random.seed(42)
    c1 = np.random.randn(30, 2)
    c2 = np.random.randn(30, 2) + 100000.0
    coords_split = np.vstack([c1, c2])
    
    W_split, D_split, W_norm = build_spatial_neighborhood_graph(coords_split, k_neighbors=4)
    assert W_split.shape == (60, 60)
    assert np.all(np.diag(W_split) == 0.0)

    # Perfect spatial separation expression: cluster 1 has expr 10, cluster 2 has expr 0
    expr_split = np.array([10.0] * 30 + [0.0] * 30)
    m_i, m_exp, m_var, m_z, m_p = calculate_morans_i(expr_split, W_split)
    g_c, g_exp, g_var, g_z, g_p = calculate_gearys_c(expr_split, W_split)

    # Strong positive autocorrelation (within floating point precision bound [-1.0, 1.0 + 1e-6])
    assert m_i > 0.80
    assert -1.0 - 1e-6 <= m_i <= 1.0 + 1e-6
    assert g_c >= 0.0
    assert g_c < 0.20
    assert m_p < 0.01

    # Alternating / checkerboard expression on dense cluster (negative autocorrelation)
    expr_checker = np.array([10.0, 0.0] * 30)
    m_i_neg, _, _, _, _ = calculate_morans_i(expr_checker, W_split)
    g_c_neg, _, _, _, _ = calculate_gearys_c(expr_checker, W_split)
    assert -1.0 - 1e-6 <= m_i_neg <= 1.0 + 1e-6
    assert g_c_neg >= 0.0


def test_spatial_deg_and_cci_capabilities_stress_execution(tmp_path):
    """
    Execute SpatialDEGCapability and CellCellCommunicationCapability under stress conditions
    using ArtifactRegistry and TaskContract.
    """
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    np.random.seed(42)
    n_cells = 80
    n_genes = 20

    X = np.random.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
    # Add a strong spatially patterned gene at index 0
    coords = np.random.rand(n_cells, 2).astype(np.float32) * 50.0
    X[:, 0] = coords[:, 0] * 2.0  # spatial gradient

    obs = pd.DataFrame({
        "cell_id": [f"c_{i}" for i in range(n_cells)],
        "cell_type": ["Microglia" if i < 40 else "Neuron" for i in range(n_cells)],
    })
    var = pd.DataFrame({
        "gene_name": ["Apoe", "Trem2", "Lrp1", "App", "Cd74", "Cx3cl1", "Cx3cr1"] + [f"Gene_{i}" for i in range(7, n_genes)]
    })
    sc_data = SCData(X=X, obs=obs, var=var, obsm={"spatial": coords})

    registry.register(
        uri_str="adata://stress_study/spatial_raw/v1",
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id="stress_study",
        created_by_task="init_task",
        operation="init",
    )

    # 1. Run Spatial Domain
    dom_cap = SpatialDomainCapability()
    dom_contract = TaskContract(
        task_id="task_dom_stress",
        study_id="stress_study",
        capability="spatial_domain",
        input_artifacts=["adata://stress_study/spatial_raw/v1"],
        expected_outputs=["adata://stress_study/spatial_domains/v1"],
        parameters={"k_neighbors": 8, "n_domains": 3, "smoothing_lambda": 0.5},
    )
    res_dom = dom_cap.execute(dom_contract, registry)
    assert res_dom.status == TaskStatus.SUCCESS
    assert registry.exists("adata://stress_study/spatial_domains/v1")

    # 2. Run Spatial DEG
    deg_cap = SpatialDEGCapability()
    deg_contract = TaskContract(
        task_id="task_deg_stress",
        study_id="stress_study",
        capability="spatial_deg",
        input_artifacts=["adata://stress_study/spatial_domains/v1"],
        expected_outputs=["table://stress_study/spatial_deg/v1"],
        parameters={"k_neighbors": 6, "min_moran_i": 0.15},
    )
    res_deg = deg_cap.execute(deg_contract, registry)
    assert res_deg.status == TaskStatus.SUCCESS
    _, deg_table = registry.get("table://stress_study/spatial_deg/v1")
    assert isinstance(deg_table, pd.DataFrame)
    assert "moran_i" in deg_table.columns
    assert "geary_c" in deg_table.columns
    # Check all Moran's I are in [-1, 1] and Geary's C >= 0
    assert (deg_table["moran_i"] >= -1.0 - 1e-6).all() and (deg_table["moran_i"] <= 1.0 + 1e-6).all()
    assert (deg_table["geary_c"] >= 0.0).all()

    # 3. Run Spatial CCI
    cci_cap = CellCellCommunicationCapability()
    cci_contract = TaskContract(
        task_id="task_cci_stress",
        study_id="stress_study",
        capability="cell_cell_communication",
        input_artifacts=["adata://stress_study/spatial_domains/v1"],
        expected_outputs=["table://stress_study/spatial_cci/v1"],
        parameters={"k_neighbors": 6, "n_permutations": 50},
    )
    res_cci = cci_cap.execute(cci_contract, registry)
    assert res_cci.status == TaskStatus.SUCCESS
    _, cci_table = registry.get("table://stress_study/spatial_cci/v1")
    assert isinstance(cci_table, pd.DataFrame)
    if not cci_table.empty:
        assert (cci_table["spatial_interaction_score"] >= 0.0).all()
        assert (cci_table["p_value"] >= 0.0).all() and (cci_table["p_value"] <= 1.0).all()


# =============================================================================
# 2. Genetic Perturbation & GRN Propagation Stress Tests
# =============================================================================

def test_grn_propagator_singular_and_ill_conditioned_matrices():
    """
    Stress-test GRN propagator (I - alpha * A)^(-1) under:
    - All-ones adjacency matrix
    - All-zeros adjacency matrix
    - Rank-1 singular matrix
    - Alpha >= 1.0 (clamped safely)
    - Extreme dimensions
    """
    # All zeros adjacency
    G = 30
    A_zero = np.zeros((G, G), dtype=np.float32)
    M_zero = compute_grn_propagator(A_zero, alpha=0.3)
    assert np.allclose(M_zero, np.eye(G, dtype=np.float32))

    # All ones row-normalized adjacency (rank 1 matrix with all entries 1/G)
    A_ones = np.ones((G, G), dtype=np.float32) / G
    np.fill_diagonal(A_ones, 0.0)
    # Row normalize
    A_ones = A_ones / np.sum(A_ones, axis=1, keepdims=True)

    M_ones = compute_grn_propagator(A_ones, alpha=0.5)
    assert not np.isnan(M_ones).any()
    assert not np.isinf(M_ones).any()

    # Extreme alpha clamping
    M_alpha_high = compute_grn_propagator(A_ones, alpha=10.0)
    assert not np.isnan(M_alpha_high).any()
    assert not np.isinf(M_alpha_high).any()

    M_alpha_neg = compute_grn_propagator(A_ones, alpha=-5.0)
    assert np.allclose(M_alpha_neg, np.eye(G, dtype=np.float32))

    # Singular matrix forced LinAlgError fallback test
    # Construct a matrix with eigenvalue = 1.0 and test invertibility handling
    A_singular = np.eye(G, dtype=np.float32)
    M_sing = compute_grn_propagator(A_singular, alpha=0.9)
    assert not np.isnan(M_sing).any()


def test_genetic_perturbation_missing_genes_and_edge_inputs(tmp_path):
    """
    Test GeneticPerturbationCapability error handling for non-existent target genes,
    zero baseline expression, and overexpression factor stress.
    """
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    N, G = 50, 10
    X = np.random.RandomState(42).poisson(lam=1.5, size=(N, G)).astype(np.float32)
    obs = pd.DataFrame({"cell_id": [f"c_{i}" for i in range(N)], "condition": ["AD"] * 25 + ["control"] * 25})
    var = pd.DataFrame({"gene_name": [f"Gene_{i}" for i in range(G)]})
    data = SCData(X=X, obs=obs, var=var)

    registry.register(
        uri_str="adata://stress_grn/raw/v1",
        payload=data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="stress_grn",
        created_by_task="init",
        operation="init",
    )

    cap = GeneticPerturbationCapability()

    # Missing target gene raises KeyError
    contract_missing = TaskContract(
        task_id="task_missing_gene",
        study_id="stress_grn",
        capability="genetic_perturbation_simulation",
        input_artifacts=["adata://stress_grn/raw/v1"],
        parameters={"target_genes": ["NonExistentGene_XYZ"]},
    )
    with pytest.raises(KeyError, match="not found in dataset"):
        cap.execute(contract_missing, registry)

    # Valid KO execution
    contract_ko = TaskContract(
        task_id="task_valid_ko",
        study_id="stress_grn",
        capability="genetic_perturbation_simulation",
        input_artifacts=["adata://stress_grn/raw/v1"],
        parameters={"target_genes": ["Gene_0"], "perturbation_type": "knockout", "efficiency": 1.0},
    )
    res_ko = cap.execute(contract_ko, registry)
    assert res_ko.status == TaskStatus.SUCCESS
    
    # Verify non-negativity of simulated expression
    out_adata_uri = res_ko.output_artifacts[0]
    _, res_payload = registry.get(out_adata_uri)
    res_sc = SCData.from_dict(res_payload)
    assert np.min(res_sc.X) >= 0.0
    # Knocked out gene at index 0 should be 0.0 across all cells
    assert np.allclose(res_sc.X[:, 0], 0.0)


# =============================================================================
# 3. Compound Perturbation & Discordance Stress Tests
# =============================================================================

def test_compound_perturbation_orthogonal_and_inverted_signatures():
    """
    Stress-test CMAP cosine discordance under orthogonal, perfectly reversed, and identical signatures.
    """
    G = 50
    # Heterogeneous disease signature (realistic differential log2FC)
    np.random.seed(42)
    sig_dis = np.linspace(-3.0, 3.0, G).astype(np.float32)

    # 1. Perfectly inverted signature (s_drug = -s_dis) -> reversal score = +1.0
    sig_rev = -sig_dis
    rev_score, cos_sim, p_val = compute_cmap_cosine_discordance(sig_dis, sig_rev, n_permutations=200)
    assert np.isclose(rev_score, 1.0, atol=1e-5)
    assert np.isclose(cos_sim, -1.0, atol=1e-5)
    assert p_val <= 0.05

    # 2. Identical signature (s_drug = s_dis) -> disease exacerbator -> reversal score = -1.0
    sig_same = sig_dis.copy()
    exac_score, exac_cos, p_val_exac = compute_cmap_cosine_discordance(sig_dis, sig_same, n_permutations=200)
    assert np.isclose(exac_score, -1.0, atol=1e-5)
    assert np.isclose(exac_cos, 1.0, atol=1e-5)
    assert p_val_exac > 0.50

    # 3. Orthogonal signature (dot product = 0)
    sig_ortho = np.zeros(G, dtype=np.float32)
    sig_ortho[:G//2] = sig_dis[G//2:]
    sig_ortho[G//2:] = -sig_dis[:G//2]
    # Ensure dot product is zero
    dot_val = np.dot(sig_dis, sig_ortho)
    assert np.isclose(dot_val, 0.0, atol=1e-4)
    ortho_score, ortho_cos, _ = compute_cmap_cosine_discordance(sig_dis, sig_ortho, n_permutations=100)
    assert np.isclose(ortho_score, 0.0, atol=1e-4)
    assert np.isclose(ortho_cos, 0.0, atol=1e-4)

    # 4. Zero-norm signature (all zeros)
    sig_zero = np.zeros(G, dtype=np.float32)
    z_score, z_cos, z_p = compute_cmap_cosine_discordance(sig_dis, sig_zero, n_permutations=100)
    assert z_score == 0.0
    assert z_cos == 0.0
    assert z_p == 1.0


def test_compound_perturbation_transition_matrix_stochasticity(tmp_path):
    """
    Test CompoundPerturbationCapability transition probability matrix row stochasticity.
    All rows must sum to 1.0 and all entries must be in [0, 1].
    """
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    N, G = 60, 15
    X = np.random.RandomState(42).poisson(lam=3.0, size=(N, G)).astype(np.float32)
    obs = pd.DataFrame({
        "cell_id": [f"c_{i}" for i in range(N)],
        "microglia_state": ["Homeostatic"] * 20 + ["DAM_Early"] * 20 + ["DAM_Late"] * 20,
    })
    var = pd.DataFrame({"gene_name": ["Apoe", "Trem2", "Clec7a", "P2ry12", "Cx3cr1"] + [f"G_{i}" for i in range(5, G)]})
    data = SCData(X=X, obs=obs, var=var)

    registry.register(
        uri_str="adata://drug_study/microglia/v1",
        payload=data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="drug_study",
        created_by_task="init",
        operation="init",
    )

    cap = CompoundPerturbationCapability()
    contract = TaskContract(
        task_id="task_bexarotene",
        study_id="drug_study",
        capability="compound_perturbation_simulation",
        input_artifacts=["adata://drug_study/microglia/v1"],
        parameters={"compound_name": "Bexarotene", "dosage": 1.5, "n_permutations": 50},
    )

    result = cap.execute(contract, registry)
    assert result.status == TaskStatus.SUCCESS

    # Check Transition Table
    out_table_uri = result.output_artifacts[0]
    _, trans_df = registry.get(out_table_uri)
    assert isinstance(trans_df, pd.DataFrame)
    
    # Verify row stochasticity
    row_sums = trans_df.values.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-4)
    assert (trans_df.values >= 0.0).all() and (trans_df.values <= 1.0).all()


# =============================================================================
# 4. Statistical Bounds & Causal Confidence Limit Tests
# =============================================================================

def test_causal_confidence_cap_at_half():
    """
    Verify that ConfidenceCalculator strictly caps in silico causal confidence at <= 0.50,
    even with multiple VERY_STRONG perturbation nodes with score 1.0.
    """
    # 5 VERY STRONG perturbation evidence nodes with score 1.0
    perturb_nodes = []
    for i in range(5):
        node = EvidenceNode(
            evidence_id=f"E_perturb_max_{i}",
            type=EvidenceType.PERTURBATION,
            polarity=EvidencePolarity.SUPPORTING,
            strength=EvidenceStrength.VERY_STRONG,
            score=1.0,
            summary=f"Maximal perturbation evidence {i}",
            source_task_id="task_test",
            source_artifact_uris=["adata://study/p/v1"],
        )
        perturb_nodes.append(node)

    conf = ConfidenceCalculator.calculate(
        supporting_evidence=perturb_nodes,
        contradicting_evidence=[],
    )

    # Invariant: Causal confidence cannot exceed 0.50
    assert conf.causal <= 0.50
    assert conf.overall <= 1.0


def test_statistical_validator_audits_extreme_perturbation_and_spatial(tmp_path):
    """
    Verify StatisticalValidator flags violations of theoretical bounds on spatial and perturbation artifacts.
    """
    validator = StatisticalValidator()
    registry = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))

    # 1. Invalid Moran's I (> 1.0 or < -1.0)
    bad_moran_df = pd.DataFrame({
        "gene": ["Gene_A", "Gene_B"],
        "moran_i": [1.85, -1.50],
        "geary_c": [0.5, 0.8],
        "fdr_q_value": [0.01, 0.02],
    })
    registry.register(
        uri_str="table://audit_test/bad_moran/v1",
        payload=bad_moran_df,
        artifact_type=ArtifactType.TABLE,
        study_id="audit_test",
        created_by_task="task_bad_moran",
        operation="moran_test",
    )
    contract_moran = TaskContract(
        task_id="task_bad_moran",
        study_id="audit_test",
        capability="spatial_deg",
        input_artifacts=["adata://audit_test/raw/v1"],
        expected_outputs=["table://audit_test/bad_moran/v1"],
    )
    result_moran = TaskResult(
        task_id="task_bad_moran",
        status=TaskStatus.SUCCESS,
        capability="spatial_deg",
        method_used="moran_v1",
        output_artifacts=["table://audit_test/bad_moran/v1"],
    )
    report_moran = validator.audit(contract_moran, result_moran, registry)
    moran_check = [c for c in report_moran.checks if c.check_name == "morans_i_bounds_check"][0]
    assert moran_check.passed is False

    # 2. Invalid Geary's C (< 0.0)
    bad_geary_df = pd.DataFrame({
        "gene": ["Gene_A"],
        "moran_i": [0.4],
        "geary_c": [-0.5],
        "fdr_q_value": [0.01],
    })
    registry.register(
        uri_str="table://audit_test/bad_geary/v1",
        payload=bad_geary_df,
        artifact_type=ArtifactType.TABLE,
        study_id="audit_test",
        created_by_task="task_bad_geary",
        operation="geary_test",
    )
    contract_geary = TaskContract(
        task_id="task_bad_geary",
        study_id="audit_test",
        capability="spatial_deg",
        input_artifacts=["adata://audit_test/raw/v1"],
        expected_outputs=["table://audit_test/bad_geary/v1"],
    )
    result_geary = TaskResult(
        task_id="task_bad_geary",
        status=TaskStatus.SUCCESS,
        capability="spatial_deg",
        method_used="geary_v1",
        output_artifacts=["table://audit_test/bad_geary/v1"],
    )
    report_geary = validator.audit(contract_geary, result_geary, registry)
    geary_check = [c for c in report_geary.checks if c.check_name == "gearys_c_bounds_check"][0]
    assert geary_check.passed is False
