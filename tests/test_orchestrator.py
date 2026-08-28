"""
Unit and integration tests for Scientific Orchestrator Intent Parsing, Planning, Routing, and DAG contracts.
"""

import pytest
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.router import CapabilityRouter
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.capabilities import create_default_capability_registry
from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec, AnalysisPolicy, Hypotheses


def test_intent_parser_preserves_scientific_objectives():
    prompt = "分析 AD 小鼠脑组织，找疾病相关小胶质细胞亚群及其演化轨迹和关键调控通路。"
    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id="AD_test_99")

    assert manifest.study_id == "AD_test_99"
    assert manifest.biological_design.species == "mus_musculus"
    assert manifest.biological_design.disease == "Alzheimer"
    assert "Microglia" in manifest.biological_design.target_cell_types
    assert manifest.analysis_policy.discovery_mode is True


def test_capability_router_3_tier_resolution():
    cap_reg = create_default_capability_registry()
    router = CapabilityRouter(cap_reg)

    manifest = IntentParser.parse_prompt_to_manifest("AD mouse brain single cell study", study_id="AD_01")
    
    # 1. Routing for DEG when donor replicates >= 3 -> Pseudobulk
    deg_method = router.resolve_method("deg", manifest, current_state={"min_replicates": 6})
    assert deg_method == "deg_pseudobulk_v1"

    # 2. Routing for Trajectory without velocity -> paga_dpt
    manifest.data.has_rna_velocity = False
    traj_method = router.resolve_method("trajectory_inference", manifest, current_state={})
    assert traj_method == "paga_dpt"

    # 3. Routing for Trajectory with velocity -> cellrank
    manifest.data.has_rna_velocity = True
    traj_method_vel = router.resolve_method("trajectory_inference", manifest, current_state={})
    assert traj_method_vel == "cellrank"

    # 4. Spatial capabilities routing
    assert router.resolve_method("spatial_domain", manifest, {}) == "spatial_domain_knn_v1"
    assert router.resolve_method("spatial_deg", manifest, {}) == "spatial_deg_morans_i_v1"
    assert router.resolve_method("cell_cell_communication", manifest, {}) == "cci_ligand_receptor_v1"

    # 5. External agent adapters routing
    assert router.resolve_method("spacell_microenvironment_analysis", manifest, {}) == "spacell_agent_v1"
    assert router.resolve_method("chatcell_dialogue_prediction", manifest, {}) == "chatcell_agent_v1"
    assert router.resolve_method("gene_function_reasoning", manifest, {}) == "gene_agent_v1"

    # 6. Knowledge Engine routing (Discovery vs Prior-guided)
    manifest.analysis_policy.prior_guided_analysis = False
    manifest.hypotheses.user_provided = []
    assert router.resolve_method("knowledge_retrieval", manifest, {}) == "knowledge_engine_discovery_v1"

    manifest.analysis_policy.prior_guided_analysis = True
    manifest.hypotheses.user_provided = ["DAM TREM2-APOE axis"]
    assert router.resolve_method("knowledge_retrieval", manifest, {}) == "knowledge_engine_prior_v1"

    # 7. Perturbation routing
    assert router.resolve_method("genetic_perturbation_simulation", manifest, {"perturbation_type": "knockout"}) == "in_silico_crispr_ko_v1"
    assert router.resolve_method("genetic_perturbation_simulation", manifest, {"perturbation_type": "overexpression"}) == "in_silico_overexpression_v1"
    assert router.resolve_method("compound_perturbation_simulation", manifest, {}) == "in_silico_compound_response_v1"


def test_capability_router_unregistered_raises_keyerror():
    cap_reg = create_default_capability_registry()
    router = CapabilityRouter(cap_reg)
    manifest = IntentParser.parse_prompt_to_manifest("AD study", study_id="AD_01")

    with pytest.raises(KeyError) as exc_info:
        router.resolve_method("unregistered_magic_capability_v99", manifest, {})
    assert "Unregistered capability" in str(exc_info.value)


def test_computational_dag_planner_contracts():
    manifest = IntentParser.parse_prompt_to_manifest("AD mouse brain single cell study", study_id="AD_01")
    tasks = ComputationalDAGPlanner.build_study_plan(manifest)

    assert len(tasks) >= 9
    task_capabilities = [t.capability for t in tasks]
    assert "dataset_audit" in task_capabilities
    assert "qc" in task_capabilities
    assert "normalization" in task_capabilities
    assert "clustering" in task_capabilities
    assert "deg" in task_capabilities
    assert "trajectory_inference" in task_capabilities

    # Verify task contract boundaries
    traj_task = next(t for t in tasks if t.capability == "trajectory_inference")
    assert "filter_cells" in traj_task.forbidden_operations
    assert "recluster" in traj_task.forbidden_operations
    assert "build_neighbor_graph" in traj_task.allowed_operations


def test_dynamic_dag_planning_spatial_and_prior_guided():
    # Multi-plane manifest with spatial coordinates and prior-guided hypothesis
    manifest = StudyManifest(
        study_id="AD_spatial_prior_study",
        title="Spatial and Prior-guided AD Study",
        biological_design=BiologicalDesign(
            species="mus_musculus",
            tissue="cortex",
            disease="Alzheimer",
            target_cell_types=["Microglia"],
        ),
        data=DataSpec(
            modalities=["scRNA", "spatial", "perturbation"],
            has_spatial_coordinates=True,
        ),
        hypotheses=Hypotheses(
            user_provided=["DAM TREM2-APOE axis activation in plaque microenvironment"],
        ),
        analysis_policy=AnalysisPolicy(
            discovery_mode=False,
            prior_guided_analysis=True,
        ),
    )

    tasks = ComputationalDAGPlanner.build_study_plan(
        manifest,
        current_state={
            "include_spatial": True,
            "include_adapters": True,
            "include_knowledge": True,
            "include_perturbation": True,
            "run_compound_perturbation": True,
        }
    )

    caps = [t.capability for t in tasks]

    # Check spatial tasks
    assert "spatial_domain" in caps
    assert "spatial_deg" in caps
    assert "cell_cell_communication" in caps

    # Check adapter tasks
    assert "spacell_microenvironment_analysis" in caps
    assert "gene_function_reasoning" in caps
    assert "chatcell_dialogue_prediction" in caps

    # Check knowledge task
    assert "knowledge_retrieval" in caps
    know_task = next(t for t in tasks if t.capability == "knowledge_retrieval")
    assert know_task.parameters["prior_guided"] is True

    # Check perturbation tasks
    assert "genetic_perturbation_simulation" in caps
    assert "compound_perturbation_simulation" in caps
    perturb_task = next(t for t in tasks if t.capability == "genetic_perturbation_simulation")
    assert perturb_task.parameters["target_gene"] == "Trem2"
    assert "recluster" in perturb_task.forbidden_operations
