"""
Unit and integration tests for In Silico Perturbation Simulation Plane.
Tests GeneticPerturbationCapability (CRISPR KO / Overexpression) and CompoundPerturbationCapability (CMAP discordance).
"""

import gc
import pytest
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import (
    EvidenceNode,
    EvidenceType,
    EvidencePolarity,
    EvidenceStrength,
    LanguageTier,
    ClaimType,
)
from eacbp.capabilities.base import BaseCapability
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.storage import ArtifactAlreadyExistsError
from eacbp.evidence.confidence import ConfidenceCalculator
from eacbp.evidence.language import LanguageEnforcer
from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.graph import EvidenceGraph

from eacbp.capabilities.perturbation.genetic import (
    GeneticPerturbationCapability,
    construct_grn_adjacency_from_data,
    compute_grn_propagator,
    generate_genetic_perturbation_evidence,
)
from eacbp.capabilities.perturbation.compound import (
    CompoundPerturbationCapability,
    compute_cmap_cosine_discordance,
    generate_compound_perturbation_evidence,
    REFERENCE_COMPOUND_DATABASE,
)


@pytest.fixture(autouse=True)
def run_gc():
    """Ensures file handles are closed after each test on Windows."""
    yield
    gc.collect()


@pytest.fixture
def test_dataset():
    """Generates synthetic single-cell dataset with distinct AD vs control states."""
    return SCData.create_synthetic_ad_study(
        n_cells=200,
        n_genes=50,
        n_ad_mice=3,
        n_ctrl_mice=3,
        random_seed=42,
    )


@pytest.fixture
def populated_registry(tmp_path, test_dataset):
    """Sets up an ArtifactRegistry containing the initial dataset artifact."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    reg.register(
        uri_str="adata://AD_TEST/microglia_subset/v1",
        payload=test_dataset.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD_TEST",
        created_by_task="task_init",
        operation="initial_subset",
    )
    return reg


def test_genetic_perturbation_crispr_ko_propagation(populated_registry):
    """Verifies in silico CRISPR KO shifts target gene towards 0 and propagates shifts across GRN."""
    reg = populated_registry
    cap = GeneticPerturbationCapability()

    contract = TaskContract(
        task_id="task_perturb_ko_trem2",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        allowed_operations=[
            "construct_grn_adjacency",
            "simulate_genetic_perturbation",
            "propagate_network_shift",
            "compute_state_reversion",
        ],
        parameters={
            "target_gene": "Trem2",
            "perturbation_type": "knockout",
            "efficiency": 0.98,
            "network_attenuation": 0.30,
        },
    )

    result = cap.execute(contract, reg)

    # 1. Execution status check
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 2
    adata_uri, table_uri = result.output_artifacts[0], result.output_artifacts[1]

    # 2. Verify Output Artifacts exist in Registry
    assert reg.exists(adata_uri)
    assert reg.exists(table_uri)

    meta_adata, payload_adata = reg.get(adata_uri)
    meta_table, payload_table = reg.get(table_uri)

    res_data = payload_adata if isinstance(payload_adata, SCData) else SCData.from_dict(payload_adata)
    res_table = payload_table if isinstance(payload_table, pd.DataFrame) else pd.DataFrame(payload_table)

    # 3. Target gene expression is knocked down near zero
    trem2_idx = list(res_data.var["gene_name"].values).index("Trem2")
    orig_meta, orig_payload = reg.get("adata://AD_TEST/microglia_subset/v1")
    orig_data = orig_payload if isinstance(orig_payload, SCData) else SCData.from_dict(orig_payload)

    orig_trem2_mean = float(np.mean(orig_data.X[:, trem2_idx]))
    pert_trem2_mean = float(np.mean(res_data.X[:, trem2_idx]))

    assert pert_trem2_mean < orig_trem2_mean * 0.05
    assert pert_trem2_mean >= 0.0

    # 4. Downstream propagation shifts other connected genes
    downstream_shifts = res_table[res_table["gene"] != "Trem2"]
    assert len(downstream_shifts) > 0
    # Top downstream genes must have non-zero expression shifts
    assert np.any(np.abs(downstream_shifts["expression_shift"].values) > 1e-4)

    # 5. Network attenuation comparison: alpha=0 vs alpha=0.3
    contract_zero_alpha = TaskContract(
        task_id="task_perturb_ko_zero_alpha",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={
            "target_gene": "Trem2",
            "perturbation_type": "knockout",
            "efficiency": 0.98,
            "network_attenuation": 0.0,
        },
    )
    res_zero = cap.execute(contract_zero_alpha, reg)
    _, table_zero_payload = reg.get(res_zero.output_artifacts[1])
    table_zero_df = table_zero_payload if isinstance(table_zero_payload, pd.DataFrame) else pd.DataFrame(table_zero_payload)
    
    # With alpha=0, downstream genes should have 0 shift
    non_target_zero_shifts = table_zero_df[table_zero_df["gene"] != "Trem2"]["expression_shift"].values
    assert np.allclose(non_target_zero_shifts, 0.0, atol=1e-5)


def test_genetic_perturbation_overexpression(populated_registry):
    """Verifies in silico overexpression increases target gene expression and computes differential shifts."""
    reg = populated_registry
    cap = GeneticPerturbationCapability()

    contract = TaskContract(
        task_id="task_perturb_oe_p2ry12",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={
            "target_gene": "P2ry12",
            "perturbation_type": "overexpression",
            "overexpression_factor": 4.0,
            "network_attenuation": 0.25,
        },
    )

    result = cap.execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS

    meta_table, payload_table = reg.get(result.output_artifacts[1])
    shift_df = payload_table if isinstance(payload_table, pd.DataFrame) else pd.DataFrame(payload_table)

    # Check that P2ry12 has positive shift
    p2ry12_row = shift_df[shift_df["gene"] == "P2ry12"].iloc[0]
    assert p2ry12_row["expression_shift"] > 0.0
    assert p2ry12_row["perturbed_mean"] > p2ry12_row["baseline_mean"]
    assert p2ry12_row["is_target_gene"] == True


def test_custom_grn_adjacency_propagation(populated_registry):
    """Tests exact mathematical propagation Delta x = (I - alpha A)^(-1) v using a controlled 3-gene GRN."""
    reg = populated_registry
    
    # 3-gene synthetic dataset with baseline expression sufficient to avoid clipping at 0
    X_simple = np.array([
        [10.0, 20.0, 10.0],
        [12.0, 20.0, 10.0],
    ], dtype=np.float32)
    obs_simple = pd.DataFrame({"cell_id": ["c1", "c2"], "condition": ["AD", "control"]})
    var_simple = pd.DataFrame({"gene_name": ["GeneA", "GeneB", "GeneC"]})
    simple_data = SCData(X=X_simple, obs=obs_simple, var=var_simple)

    reg.register(
        uri_str="adata://AD_TEST/simple_grn/v1",
        payload=simple_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD_TEST",
        created_by_task="task_init_simple",
        operation="init_simple",
    )

    # Custom adjacency: GeneA -> GeneB (weight 1.0), GeneB -> GeneC (weight 1.0)
    # A = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
    custom_A = np.array([
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0],
    ], dtype=np.float32)

    alpha = 0.5
    # (I - alpha A) = [[1, -0.5, 0], [0, 1, -0.5], [0, 0, 1]]
    # (I - alpha A)^(-1) = [[1, 0.5, 0.25], [0, 1, 0.5], [0, 0, 1]]
    cap = GeneticPerturbationCapability()
    contract = TaskContract(
        task_id="task_custom_grn",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/simple_grn/v1"],
        parameters={
            "target_gene": "GeneA",
            "perturbation_type": "knockout",
            "efficiency": 1.0,
            "network_attenuation": alpha,
            "grn_adjacency": custom_A.tolist(),
        },
    )

    result = cap.execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS

    _, payload_adata = reg.get(result.output_artifacts[0])
    res_data = payload_adata if isinstance(payload_adata, SCData) else SCData.from_dict(payload_adata)

    # Base perturbation on GeneA: v = [-10, 0, 0] for cell 0, [-12, 0, 0] for cell 1
    # Propagated shift: GeneB receives -11.0 * 0.5 = -5.5
    # GeneC receives -11.0 * 0.25 = -2.75
    expected_shift_gene_b = -1.0 * float(np.mean(X_simple[:, 0])) * 0.5
    expected_shift_gene_c = -1.0 * float(np.mean(X_simple[:, 0])) * 0.25

    _, payload_table = reg.get(result.output_artifacts[1])
    shift_df = payload_table if isinstance(payload_table, pd.DataFrame) else pd.DataFrame(payload_table)

    gene_b_shift = shift_df[shift_df["gene"] == "GeneB"]["expression_shift"].iloc[0]
    gene_c_shift = shift_df[shift_df["gene"] == "GeneC"]["expression_shift"].iloc[0]

    assert pytest.approx(gene_b_shift, rel=1e-3) == expected_shift_gene_b
    assert pytest.approx(gene_c_shift, rel=1e-3) == expected_shift_gene_c


def test_compound_perturbation_therapeutic_reversal(populated_registry):
    """Verifies CMAP-style cosine discordance and counterfactual cell state transitions for therapeutic compounds."""
    reg = populated_registry
    cap = CompoundPerturbationCapability()

    contract = TaskContract(
        task_id="task_compound_bexarotene",
        capability="compound_perturbation_simulation",
        method="in_silico_compound_response_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        allowed_operations=[
            "compute_disease_signature",
            "calculate_cmap_discordance",
            "simulate_counterfactual_transitions",
            "compute_transition_matrix",
        ],
        parameters={
            "compound_name": "Bexarotene",
            "dosage": 1.0,
            "n_permutations": 300,
        },
    )

    result = cap.execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 2

    table_uri, adata_uri = result.output_artifacts[0], result.output_artifacts[1]
    assert reg.exists(table_uri)
    assert reg.exists(adata_uri)

    _, table_payload = reg.get(table_uri)
    trans_matrix = table_payload if isinstance(table_payload, pd.DataFrame) else pd.DataFrame(table_payload)

    # 1. Bexarotene reverses DAM signature -> positive discordance score
    assert result.metrics["reversal_score"] > 0.0
    assert result.metrics["cosine_similarity"] < 0.0
    assert result.metrics["p_value"] <= 1.0

    # 2. Transition matrix is row-stochastic (rows sum to ~1.0)
    row_sums = trans_matrix.sum(axis=1).values
    for r_sum in row_sums:
        assert pytest.approx(r_sum, rel=1e-2) == 1.0

    # 3. Output AnnData contains counterfactual predictions
    _, adata_payload = reg.get(adata_uri)
    drug_data = adata_payload if isinstance(adata_payload, SCData) else SCData.from_dict(adata_payload)
    assert "predicted_transition_state" in drug_data.obs.columns


def test_compound_perturbation_disease_exacerbator(populated_registry):
    """Verifies that disease-exacerbating compounds yield negative discordance scores and flag no therapeutic potential."""
    reg = populated_registry
    cap = CompoundPerturbationCapability()

    contract = TaskContract(
        task_id="task_compound_exacerbator",
        capability="compound_perturbation_simulation",
        method="in_silico_compound_response_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={
            "compound_name": "Mock_Exacerbator",
            "dosage": 1.0,
            "n_permutations": 200,
        },
    )

    result = cap.execute(contract, reg)
    assert result.status == TaskStatus.SUCCESS

    # Mock_Exacerbator aligns with disease signature -> negative discordance score
    assert result.metrics["reversal_score"] < 0.0
    assert result.metrics["therapeutic_potential"] is False


def test_calibrated_evidence_generation_and_confidence_bounds(populated_registry):
    """Verifies that EvidenceType.PERTURBATION nodes are strictly capped at causal confidence <= 0.50."""
    reg = populated_registry
    cap = GeneticPerturbationCapability()

    contract = TaskContract(
        task_id="task_perturb_ko_calib",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={"target_gene": "Trem2", "perturbation_type": "knockout", "efficiency": 1.0},
    )
    result = cap.execute(contract, reg)

    # Generate calibrated evidence node
    ev = generate_genetic_perturbation_evidence(
        contract=contract,
        result=result,
        target_gene="Trem2",
        reversion_rate=0.85,  # Even with high raw reversion rate 85%
        perturbation_type="knockout",
    )

    # 1. Type and score verification
    assert ev.type == EvidenceType.PERTURBATION
    assert ev.score <= 0.50  # Must be capped at 0.50!
    assert ev.score > 0.0
    assert ev.biological_context["causal_status"] == "in_silico_perturbed"

    # 2. Test Compound Evidence Generation
    comp_cap = CompoundPerturbationCapability()
    comp_contract = TaskContract(
        task_id="task_comp_calib",
        capability="compound_perturbation_simulation",
        method="in_silico_compound_response_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={"compound_name": "Bexarotene"},
    )
    comp_result = comp_cap.execute(comp_contract, reg)
    
    comp_ev = generate_compound_perturbation_evidence(
        contract=comp_contract,
        result=comp_result,
        compound_name="Bexarotene",
        reversal_score=0.92,
        transition_rate=0.75,
    )
    assert comp_ev.type == EvidenceType.PERTURBATION
    assert comp_ev.score <= 0.50  # Must be capped at 0.50!

    # 3. Integrate with ConfidenceCalculator
    deg_ev = EvidenceNode(
        evidence_id="E_deg_001",
        type=EvidenceType.PSEUDOBULK_DEG,
        strength=EvidenceStrength.STRONG,
        score=0.95,
        summary="Trem2 is upregulated in AD.",
        source_task_id="task_deg_001",
    )

    conf_score = ConfidenceCalculator.calculate(
        supporting_evidence=[deg_ev, ev],
        contradicting_evidence=[],
    )

    # Causal score must be strictly <= 0.50
    assert conf_score.causal > 0.0
    assert conf_score.causal <= 0.50
    assert conf_score.association > 0.60
    assert conf_score.overall > 0.0


def test_language_enforcer_guardrail_with_perturbation():
    """Verifies that in silico perturbation permits Level 4 hypotheses but prevents observational causal overclaims."""
    # 1. Observational claim using forbidden causal verb "drives" must fail Level 2
    valid_obs, error_msg = LanguageEnforcer.audit_statement(
        statement="Trem2 expression drives the transition from homeostatic to DAM microglia.",
        tier=LanguageTier.LEVEL_2_STATISTICAL_INFERENCE,
        causal_status="observational",
    )
    assert not valid_obs
    assert "EPISTEMIC VIOLATION" in error_msg

    # 2. In silico perturbation hypothesis phrased under Level 4 passes
    valid_hyp, _ = LanguageEnforcer.audit_statement(
        statement="In silico knockout simulations suggest Trem2 may participate in maintaining the DAM state.",
        tier=LanguageTier.LEVEL_4_HYPOTHESIS,
        causal_status="observational",
    )
    assert valid_hyp


def test_versioned_artifact_lineage_and_immutability(populated_registry):
    """Verifies SHA-256 content addressing, lineage tracking, and raw data immutability."""
    reg = populated_registry
    cap_reg = CapabilityRegistry()
    cap_reg.register(GeneticPerturbationCapability())
    cap_reg.register(CompoundPerturbationCapability())

    # Check original payload before execution
    orig_meta, orig_payload = reg.get("adata://AD_TEST/microglia_subset/v1")
    orig_hash = orig_meta.sha256_hash

    contract = TaskContract(
        task_id="task_lineage_test",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        allowed_operations=[
            "construct_grn_adjacency",
            "simulate_genetic_perturbation",
            "propagate_network_shift",
            "compute_state_reversion",
        ],
        forbidden_operations=["filter_cells", "recluster"],
        parameters={"target_gene": "Trem2", "perturbation_type": "knockout"},
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.SUCCESS

    # 1. Original input artifact is completely unchanged (Invariant 1: Raw data immutability)
    after_meta, _ = reg.get("adata://AD_TEST/microglia_subset/v1")
    assert after_meta.sha256_hash == orig_hash

    # 2. Output artifact has complete parent lineage (Invariant 2: Lineage tracking)
    out_uri = result.output_artifacts[0]
    out_meta, out_payload = reg.get(out_uri)
    assert "adata://AD_TEST/microglia_subset/v1" in out_meta.parent_uris
    assert out_meta.created_by_task == "task_lineage_test"
    assert out_meta.sha256_hash is not None and out_meta.sha256_hash.startswith("sha256:")

    # 3. Attempting to overwrite existing artifact URI without version increment raises error
    with pytest.raises(ArtifactAlreadyExistsError):
        reg.register(
            uri_str=out_uri,
            payload=out_payload,
            artifact_type=ArtifactType.ANNDATA,
            study_id="AD_TEST",
            created_by_task="task_overwrite_attempt",
            operation="unauthorized_overwrite",
        )


def test_edge_cases_and_error_handling(populated_registry):
    """Tests error handling for missing target genes, invalid perturbation types, and dimension mismatches."""
    reg = populated_registry
    cap = GeneticPerturbationCapability()

    # 1. Missing target gene
    bad_gene_contract = TaskContract(
        task_id="task_bad_gene",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={"target_gene": "NonExistentGene_9999", "perturbation_type": "knockout"},
    )
    with pytest.raises(KeyError) as exc_info:
        cap.execute(bad_gene_contract, reg)
    assert "NonExistentGene_9999" in str(exc_info.value)

    # 2. Invalid perturbation type
    bad_type_contract = TaskContract(
        task_id="task_bad_type",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={"target_gene": "Trem2", "perturbation_type": "unsupported_magic_shift"},
    )
    with pytest.raises(ValueError) as exc_info:
        cap.execute(bad_type_contract, reg)
    assert "Unsupported perturbation_type" in str(exc_info.value)

    # 3. Custom adjacency matrix dimension mismatch
    bad_adj_contract = TaskContract(
        task_id="task_bad_adj",
        capability="genetic_perturbation_simulation",
        method="in_silico_crispr_ko_v1",
        input_artifacts=["adata://AD_TEST/microglia_subset/v1"],
        parameters={
            "target_gene": "Trem2",
            "perturbation_type": "knockout",
            "grn_adjacency": [[0.0, 1.0], [1.0, 0.0]],  # 2x2 matrix for 50-gene dataset
        },
    )
    with pytest.raises(ValueError) as exc_info:
        cap.execute(bad_adj_contract, reg)
    assert "Custom GRN adjacency matrix shape" in str(exc_info.value)
