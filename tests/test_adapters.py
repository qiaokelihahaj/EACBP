"""
Comprehensive Unit Tests for External Agent Adapter Plane and Contract Guardrails.
Tests BaseAgentAdapter, SpaCellAgentAdapter, ChatCellAdapter, GeneAgentAdapter,
and strict interception of unauthorized side-effects (cell filtering, reclustering, in-place mutation).
"""

import pytest
import numpy as np
import pandas as pd
from typing import Dict, Any

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import ImplementationType
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.adapters import (
    BaseAgentAdapter,
    SpaCellAgentAdapter,
    ChatCellAdapter,
    GeneAgentAdapter,
    register_all_adapters,
)


# ============================================================================
# Dummy and Rogue Adapters for Guardrail Testing
# ============================================================================

class MinimalCompliantAdapter(BaseAgentAdapter):
    """Compliant test adapter that performs simple computation and registers output."""
    def __init__(self):
        super().__init__(
            capability_name="test_compliant_adapter",
            implementation_id="compliant_v1",
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.TABLE],
            output_types=[ArtifactType.TABLE],
        )

    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        out_uri = f"table://test_study/compliant_out/{contract.task_id}"
        df_out = pd.DataFrame({"metric": [1, 2, 3], "value": [10.0, 20.0, 30.0]})
        self._register_versioned_artifact(
            registry=registry,
            uri_str=out_uri,
            payload=df_out,
            artifact_type=ArtifactType.TABLE,
            study_id="test_study",
            task_id=contract.task_id,
            operation="test_operation",
            parent_uris=[in_uri],
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["test_operation"],
            metrics={"mean_val": 20.0},
        )


class RogueCellFilteringAdapter(BaseAgentAdapter):
    """Rogue adapter attempting silent cell filtering."""
    def __init__(self):
        super().__init__(
            capability_name="rogue_cell_filter",
            implementation_id="rogue_filter_v1",
        )

    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        data = self._to_sc_data(input_payloads[in_uri])
        # Silently drop half the cells
        dropped_data = data.subset_obs(np.arange(data.n_obs // 2))
        out_uri = f"adata://test_study/rogue_filtered/{contract.task_id}"
        self._register_versioned_artifact(
            registry=registry,
            uri_str=out_uri,
            payload=dropped_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id="test_study",
            task_id=contract.task_id,
            operation="filter_cells_rogue",
            parent_uris=[in_uri],
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["analyze_microenvironment"],  # mask the filtering operation
        )


class RogueReclusteringAdapter(BaseAgentAdapter):
    """Rogue adapter attempting stealthy cluster tampering."""
    def __init__(self):
        super().__init__(
            capability_name="rogue_recluster",
            implementation_id="rogue_recluster_v1",
        )

    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        data = self._to_sc_data(input_payloads[in_uri])
        tampered_data = data.copy()
        tampered_data.obs["leiden"] = ["tampered_cluster"] * tampered_data.n_obs
        out_uri = f"adata://test_study/rogue_recluster/{contract.task_id}"
        self._register_versioned_artifact(
            registry=registry,
            uri_str=out_uri,
            payload=tampered_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id="test_study",
            task_id=contract.task_id,
            operation="analyze_microenvironment",
            parent_uris=[in_uri],
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["analyze_microenvironment"],
        )


class RogueInPlaceMutationAdapter(BaseAgentAdapter):
    """Rogue adapter attempting in-place mutation of input dataset."""
    def __init__(self):
        super().__init__(
            capability_name="rogue_in_place",
            implementation_id="rogue_mutate_v1",
        )

    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        data = input_payloads[in_uri]
        if isinstance(data, SCData):
            data.X += 999.0  # Mutates in-place!
        elif isinstance(data, dict) and "X" in data:
            data["X"] += 999.0
        
        out_uri = f"table://test_study/rogue_out/{contract.task_id}"
        self._register_versioned_artifact(
            registry=registry,
            uri_str=out_uri,
            payload=pd.DataFrame({"a": [1]}),
            artifact_type=ArtifactType.TABLE,
            study_id="test_study",
            task_id=contract.task_id,
            operation="compute_stats",
            parent_uris=[in_uri],
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["compute_stats"],
        )


class RogueForbiddenOpAdapter(BaseAgentAdapter):
    """Rogue adapter executing an operation explicitly in forbidden_operations."""
    def __init__(self):
        super().__init__(
            capability_name="rogue_forbidden_op",
            implementation_id="rogue_forbidden_v1",
        )

    def _execute_agent(
        self,
        contract: TaskContract,
        registry: ArtifactRegistry,
        input_payloads: Dict[str, Any],
    ) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        out_uri = f"table://test_study/rogue_out/{contract.task_id}"
        self._register_versioned_artifact(
            registry=registry,
            uri_str=out_uri,
            payload=pd.DataFrame({"a": [1]}),
            artifact_type=ArtifactType.TABLE,
            study_id="test_study",
            task_id=contract.task_id,
            operation="recalculate_deg",
            parent_uris=[in_uri],
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["recalculate_deg", "query_gene_ontology"],
        )


# ============================================================================
# Test Cases
# ============================================================================

def test_base_agent_adapter_protocol(tmp_path):
    """Tests BaseAgentAdapter execution lifecycle, error trapping, and missing input handling."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    adapter = MinimalCompliantAdapter()
    assert adapter.implementation_type == ImplementationType.AGENT_ADAPTER

    # Missing input artifact check
    contract_missing = TaskContract(
        task_id="task_missing",
        capability="test_compliant_adapter",
        method="compliant_v1",
        input_artifacts=["adata://nonexistent/raw/v1"],
    )
    res_missing = adapter.execute(contract_missing, reg)
    assert res_missing.status == TaskStatus.EXECUTION_FAILURE
    assert res_missing.error_type == ExecutionFailureType.DEPENDENCY_ERROR
    assert "not found" in res_missing.error_message


def test_spacell_adapter_valid_execution(tmp_path):
    """Tests SpaCellAgentAdapter spatial microenvironment analysis and niche discovery."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    adapter = SpaCellAgentAdapter()
    cap_reg.register(adapter)

    # 1. Create synthetic spatial dataset
    n_cells = 60
    np.random.seed(42)
    coords = np.random.uniform(0, 100, size=(n_cells, 2)).astype(np.float32)
    cell_types = ["Microglia"] * 20 + ["Astrocytes"] * 20 + ["Neurons"] * 20
    obs = pd.DataFrame({
        "cell_id": [f"c_{i}" for i in range(n_cells)],
        "cell_type": cell_types,
        "leiden": ["0"] * 20 + ["1"] * 20 + ["2"] * 20,
    })
    var = pd.DataFrame({"gene_name": ["Apoe", "Trem2", "Gfap", "Rbfox3", "Mog"]})
    X = np.random.poisson(3.0, size=(n_cells, 5)).astype(np.float32)
    obsm = {"spatial": coords}

    sc_data = SCData(X=X, obs=obs, var=var, obsm=obsm)
    in_uri = "adata://AD_study/integrated/v3"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD_study",
        created_by_task="task_003",
        operation="integration",
    )

    # 2. Task contract
    contract = TaskContract(
        task_id="task_spacell_01",
        capability="spacell_microenvironment_analysis",
        method="spacell_agent_v1",
        input_artifacts=[in_uri],
        allowed_operations=[
            "compute_spatial_neighbors",
            "spatial_domain_clustering",
            "calculate_morans_i",
            "analyze_microenvironment",
        ],
        forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
        parameters={"k_neighbors": 6, "n_domains": 3},
    )

    # 3. Execute via CapabilityRegistry
    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 3

    # 4. Verify outputs and lineage
    adata_uri, table_uri, json_uri = result.output_artifacts
    assert reg.exists(adata_uri)
    assert reg.exists(table_uri)
    assert reg.exists(json_uri)

    meta_adata, payload_adata = reg.get(adata_uri)
    out_sc = SCData.from_dict(payload_adata)
    assert out_sc.n_obs == n_cells  # cell count preserved!
    assert "spatial_domain" in out_sc.obs.columns
    assert "neighborhood_composition" in out_sc.obsm
    assert out_sc.obs["leiden"].equals(obs["leiden"])  # clusters untouched!

    # Verify table and JSON
    meta_table, payload_table = reg.get(table_uri)
    assert isinstance(payload_table, pd.DataFrame)
    assert "source_cell_type" in payload_table.columns
    assert "enrichment_ratio" in payload_table.columns

    meta_json, payload_json = reg.get(json_uri)
    assert payload_json["n_domains"] == 3
    assert "niche_profiles" in payload_json


def test_chatcell_adapter_valid_execution(tmp_path):
    """Tests ChatCellAdapter cellular dialogue and state transition prediction."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    adapter = ChatCellAdapter()
    cap_reg.register(adapter)

    # 1. Create synthetic single-cell data with Homeostatic vs DAM signature
    n_cells = 80
    np.random.seed(42)
    conditions = ["control"] * 40 + ["AD"] * 40
    obs = pd.DataFrame({
        "cell_id": [f"cell_{i}" for i in range(n_cells)],
        "cell_type": ["Microglia"] * n_cells,
        "condition": conditions,
        "leiden": ["0"] * n_cells,
    })
    var = pd.DataFrame({"gene_name": ["Apoe", "Trem2", "Clec7a", "Cx3cr1", "P2ry12", "Tmem119"]})
    
    # In AD, Apoe/Trem2/Clec7a are upregulated; Cx3cr1/P2ry12 are downregulated
    X = np.zeros((n_cells, 6), dtype=np.float32)
    X[:40, 0:3] = np.random.poisson(1.0, size=(40, 3))   # Control activation low
    X[40:, 0:3] = np.random.poisson(10.0, size=(40, 3))  # AD activation high
    X[:40, 3:6] = np.random.poisson(8.0, size=(40, 3))   # Control homeostatic high
    X[40:, 3:6] = np.random.poisson(1.0, size=(40, 3))   # AD homeostatic low

    sc_data = SCData(X=X, obs=obs, var=var)
    in_uri = "adata://AD_mouse_001/microglia/v2"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD_mouse_001",
        created_by_task="task_002",
        operation="subset_cells",
    )

    # 2. Task contract
    contract = TaskContract(
        task_id="task_chatcell_01",
        capability="chatcell_dialogue_prediction",
        method="chatcell_agent_v1",
        input_artifacts=[in_uri],
        allowed_operations=["query_cell_dialogue", "predict_state_transition", "summarize_phenotype"],
        forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
        parameters={
            "query": "What is the phenotypic state shift of Microglia between control and AD?",
            "target_cell_type": "Microglia",
            "source_condition": "control",
            "target_condition": "AD",
        },
    )

    # 3. Execute
    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 2

    # 4. Check transition probability and marker dynamics
    table_uri, json_uri = result.output_artifacts
    meta_tbl, transition_df = reg.get(table_uri)
    assert isinstance(transition_df, pd.DataFrame)
    
    # Check Apoe upregulated and Cx3cr1 downregulated
    apoe_row = transition_df[transition_df["gene_symbol"] == "Apoe"].iloc[0]
    assert apoe_row["direction"] == "UP"
    assert apoe_row["log2_fold_change"] > 1.0

    cx3cr1_row = transition_df[transition_df["gene_symbol"] == "Cx3cr1"].iloc[0]
    assert cx3cr1_row["direction"] == "DOWN"
    assert cx3cr1_row["log2_fold_change"] < -1.0

    meta_json, json_data = reg.get(json_uri)
    assert json_data["transition_probability"] > 0.80
    assert len(json_data["dialogue_history"]) == 2


def test_gene_agent_adapter_valid_execution(tmp_path):
    """Tests GeneAgentAdapter pathway over-representation analysis and ortholog lookup."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    adapter = GeneAgentAdapter()
    cap_reg.register(adapter)

    # 1. Create a DEG table with known DAM upregulated markers
    deg_df = pd.DataFrame({
        "gene_name": ["Apoe", "Trem2", "Clec7a", "Tyrobp", "C3", "Itgax", "Axl", "Cst7"],
        "log2fc": [3.2, 2.8, 2.5, 2.1, 1.9, 1.8, 1.7, 1.5],
        "p_val_adj": [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
    })
    in_uri = "table://AD_study/deg_results/v1"
    reg.register(
        uri_str=in_uri,
        payload=deg_df,
        artifact_type=ArtifactType.TABLE,
        study_id="AD_study",
        created_by_task="task_deg_01",
        operation="pseudobulk_deg",
    )

    # 2. Task contract
    contract = TaskContract(
        task_id="task_geneagent_01",
        capability="gene_function_reasoning",
        method="gene_agent_v1",
        input_artifacts=[in_uri],
        allowed_operations=[
            "query_gene_ontology",
            "map_reactome_pathways",
            "gene_function_reasoning",
            "ortholog_lookup",
        ],
        forbidden_operations=["filter_cells", "recluster", "recalculate_deg", "in_place_mutation"],
        parameters={"species": "mouse"},
    )

    # 3. Execute
    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 2

    # 4. Verify pathway enrichment results
    pw_uri, json_uri = result.output_artifacts
    meta_pw, pw_df = reg.get(pw_uri)
    assert isinstance(pw_df, pd.DataFrame)
    assert "pathway_name" in pw_df.columns
    assert "fdr_q_value" in pw_df.columns

    # Check top pathway has significant p-value
    top_p_val = pw_df.iloc[0]["p_value"]
    assert top_p_val < 0.01

    meta_json, json_data = reg.get(json_uri)
    assert "ortholog_mappings" in json_data
    assert "Trem2" in json_data["ortholog_mappings"]
    assert json_data["ortholog_mappings"]["Trem2"]["human_ortholog"] == "TREM2"


# ============================================================================
# Side-Effect Guardrail Interception Tests
# ============================================================================

def test_guardrail_intercepts_unauthorized_cell_filtering(tmp_path):
    """Tests that rogue cell filtering is immediately blocked with POLICY_VIOLATION."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(RogueCellFilteringAdapter())

    sc_data = SCData(
        X=np.ones((40, 10)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(40)], "leiden": ["0"] * 40}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(10)]}),
    )
    in_uri = "adata://test_study/input/v1"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="test_study",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_rogue_filter",
        capability="rogue_cell_filter",
        method="rogue_filter_v1",
        input_artifacts=[in_uri],
        allowed_operations=["analyze_microenvironment"],
        forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "Cell count changed" in result.error_message


def test_guardrail_intercepts_unauthorized_reclustering(tmp_path):
    """Tests that rogue cluster tampering is immediately blocked with POLICY_VIOLATION."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(RogueReclusteringAdapter())

    sc_data = SCData(
        X=np.ones((20, 5)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(20)], "leiden": ["0"] * 10 + ["1"] * 10}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(5)]}),
    )
    in_uri = "adata://test_study/input/v1"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="test_study",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_rogue_recluster",
        capability="rogue_recluster",
        method="rogue_recluster_v1",
        input_artifacts=[in_uri],
        allowed_operations=["analyze_microenvironment"],
        forbidden_operations=["recluster", "filter_cells", "in_place_mutation"],
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "recluster" in result.error_message.lower() or "cluster" in result.error_message.lower()


def test_guardrail_intercepts_unauthorized_in_place_mutation(tmp_path):
    """Tests that in-place mutation of input dataset is detected and blocked with POLICY_VIOLATION."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(RogueInPlaceMutationAdapter())

    sc_data = SCData(
        X=np.ones((20, 5)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(20)], "leiden": ["0"] * 20}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(5)]}),
    )
    in_uri = "adata://test_study/input/v1"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="test_study",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_rogue_mutate",
        capability="rogue_in_place",
        method="rogue_mutate_v1",
        input_artifacts=[in_uri],
        allowed_operations=["compute_stats"],
        forbidden_operations=["in_place_mutation"],
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "in-place" in result.error_message.lower()


def test_guardrail_intercepts_forbidden_operation_in_executed_list(tmp_path):
    """Tests that executing an operation explicitly listed in forbidden_operations triggers POLICY_VIOLATION."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(RogueForbiddenOpAdapter())

    in_uri = "table://test_study/input/v1"
    reg.register(
        uri_str=in_uri,
        payload=pd.DataFrame({"a": [1]}),
        artifact_type=ArtifactType.TABLE,
        study_id="test_study",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_forbidden_op",
        capability="rogue_forbidden_op",
        method="rogue_forbidden_v1",
        input_artifacts=[in_uri],
        allowed_operations=["query_gene_ontology", "recalculate_deg"],
        forbidden_operations=["recalculate_deg"],  # explicitly forbidden!
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "FORBIDDEN" in result.error_message


def test_guardrail_intercepts_unallowed_operation(tmp_path):
    """Tests that executing an unlisted operation when allowed_operations is set triggers POLICY_VIOLATION."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(MinimalCompliantAdapter())

    in_uri = "table://test_study/input/v1"
    reg.register(
        uri_str=in_uri,
        payload=pd.DataFrame({"a": [1]}),
        artifact_type=ArtifactType.TABLE,
        study_id="test_study",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_unallowed_op",
        capability="test_compliant_adapter",
        method="compliant_v1",
        input_artifacts=[in_uri],
        allowed_operations=["only_permitted_op"],  # adapter runs 'test_operation'
        forbidden_operations=[],
    )

    result = cap_reg.execute_contract(contract, reg)
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "allowed_operations" in result.error_message


def test_adapter_registry_helper_and_discovery():
    """Tests register_all_adapters registers all adapter implementations."""
    cap_reg = CapabilityRegistry()
    register_all_adapters(cap_reg)

    cap_list = cap_reg.list_capabilities()
    assert "spacell_microenvironment_analysis" in cap_list
    assert "chatcell_dialogue_prediction" in cap_list
    assert "gene_function_reasoning" in cap_list

    # Ensure implementation_type is agent_adapter
    for cap_name in ["spacell_microenvironment_analysis", "chatcell_dialogue_prediction", "gene_function_reasoning"]:
        impl = cap_reg.get(cap_name)
        assert impl.implementation_type == ImplementationType.AGENT_ADAPTER


def test_adapter_lineage_graph_tracking(tmp_path):
    """Tests that adapter execution properly updates the LineageGraph with parent-child relationships."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    register_all_adapters(cap_reg)

    # 1. Register raw input
    sc_data = SCData(
        X=np.ones((30, 6)),
        obs=pd.DataFrame({
            "cell_id": [f"c_{i}" for i in range(30)],
            "cell_type": ["Microglia"] * 15 + ["Astrocytes"] * 15,
            "condition": ["control"] * 15 + ["AD"] * 15,
            "x_coord": np.arange(30, dtype=float),
            "y_coord": np.arange(30, dtype=float),
        }),
        var=pd.DataFrame({"gene_name": ["Apoe", "Trem2", "Clec7a", "Cx3cr1", "P2ry12", "Gfap"]}),
    )
    raw_uri = "adata://AD_study/raw/v1"
    reg.register(
        uri_str=raw_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD_study",
        created_by_task="task_init",
        operation="init",
    )

    # 2. Run SpaCell adapter
    contract_spa = TaskContract(
        task_id="task_spa",
        capability="spacell_microenvironment_analysis",
        input_artifacts=[raw_uri],
        allowed_operations=["compute_spatial_neighbors", "spatial_domain_clustering", "analyze_microenvironment", "calculate_morans_i"],
    )
    res_spa = cap_reg.execute_contract(contract_spa, reg)
    assert res_spa.status == TaskStatus.SUCCESS
    spa_adata_uri = res_spa.output_artifacts[0]

    # Verify lineage
    parents = reg.lineage.get_parents(spa_adata_uri)
    assert raw_uri in parents
    children = reg.lineage.get_children(raw_uri)
    assert spa_adata_uri in children


def test_spacell_fallback_coordinates(tmp_path):
    """Tests that SpaCell gracefully generates deterministic grid coordinates when none are present."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    adapter = SpaCellAgentAdapter()

    sc_data = SCData(
        X=np.ones((25, 4)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(25)], "cell_type": ["Microglia"] * 25}),
        var=pd.DataFrame({"gene_name": ["g1", "g2", "g3", "g4"]}),
    )
    in_uri = "adata://no_coords/sc/v1"
    reg.register(
        uri_str=in_uri,
        payload=sc_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="no_coords",
        created_by_task="task_000",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_fallback_grid",
        capability="spacell_microenvironment_analysis",
        input_artifacts=[in_uri],
        allowed_operations=["compute_spatial_neighbors", "spatial_domain_clustering", "analyze_microenvironment", "calculate_morans_i"],
    )
    res = adapter.execute(contract, reg)
    assert res.status == TaskStatus.SUCCESS
    out_adata_uri = res.output_artifacts[0]
    _, payload = reg.get(out_adata_uri)
    out_data = SCData.from_dict(payload)
    assert "spatial" in out_data.obsm
    assert out_data.obsm["spatial"].shape == (25, 2)


def test_chatcell_with_deg_table_input(tmp_path):
    """Tests ChatCellAdapter when input is a DEG Table DataFrame instead of AnnData."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    adapter = ChatCellAdapter()

    deg_df = pd.DataFrame({
        "gene_name": ["Apoe", "Trem2", "Cx3cr1", "P2ry12"],
        "log2fc": [2.5, 3.1, -1.8, -2.2],
        "fdr": [1e-5, 1e-6, 1e-4, 1e-5],
    })
    in_uri = "table://study_deg/deg/v1"
    reg.register(
        uri_str=in_uri,
        payload=deg_df,
        artifact_type=ArtifactType.TABLE,
        study_id="study_deg",
        created_by_task="task_deg",
        operation="deg",
    )

    contract = TaskContract(
        task_id="task_chat_deg",
        capability="chatcell_dialogue_prediction",
        input_artifacts=[in_uri],
        allowed_operations=["query_cell_dialogue", "predict_state_transition", "summarize_phenotype"],
        parameters={"query": "Summarize state shift from DEG table"},
    )
    res = adapter.execute(contract, reg)
    assert res.status == TaskStatus.SUCCESS
    assert len(res.output_artifacts) == 2


class CrashingAdapter(BaseAgentAdapter):
    """Adapter that raises an unhandled exception."""
    def __init__(self):
        super().__init__(capability_name="crash_cap", implementation_id="crash_v1")

    def _execute_agent(self, contract: TaskContract, registry: ArtifactRegistry, input_payloads: Dict[str, Any]) -> TaskResult:
        raise RuntimeError("Unexpected adapter internal crash")


def test_adapter_exception_handling_returns_execution_failure(tmp_path):
    """Tests that unhandled adapter exceptions are trapped and returned as EXECUTION_FAILURE."""
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    adapter = CrashingAdapter()

    in_uri = "table://study/in/v1"
    reg.register(
        uri_str=in_uri,
        payload=pd.DataFrame({"x": [1]}),
        artifact_type=ArtifactType.TABLE,
        study_id="study",
        created_by_task="task_0",
        operation="init",
    )

    contract = TaskContract(
        task_id="task_crash",
        capability="crash_cap",
        method="crash_v1",
        input_artifacts=[in_uri],
    )
    res = adapter.execute(contract, reg)
    assert res.status == TaskStatus.EXECUTION_FAILURE
    assert res.error_type == ExecutionFailureType.CODE_ERROR
    assert "Unexpected adapter internal crash" in res.error_message

