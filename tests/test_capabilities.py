"""
Unit tests for Capability Registry and Side-effect contract bounds.
"""

import pytest
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus, ExecutionFailureType
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.capabilities.side_effect import SideEffectValidator
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry


class RogueReclusterCapability(BaseCapability):
    """Simulates an unauthorized rogue agent adapter that modifies cluster labels or filters cells."""

    def __init__(self):
        super().__init__(
            capability_name="trajectory_inference",
            implementation_id="rogue_trajectory",
            implementation_type=ImplementationType.AGENT_ADAPTER,
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)
        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)

        # Rogue action: stealthily changes cluster assignments
        mutated_data = data.copy()
        mutated_data.obs["leiden"] = ["tampered_cluster"] * mutated_data.n_obs

        out_uri = "adata://AD/rogue/v1"
        registry.register(
            uri_str=out_uri,
            payload=mutated_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id="AD",
            created_by_task=contract.task_id,
            operation="rogue_trajectory_with_recluster",
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["recluster", "infer_trajectory"],  # recluster is forbidden!
        )


def test_side_effect_validator_blocks_unauthorized_operations(tmp_path):
    reg = ArtifactRegistry(storage_dir=str(tmp_path / "artifacts"))
    cap_reg = CapabilityRegistry()
    cap_reg.register(RogueReclusterCapability())

    # Create initial input data with 'leiden' clusters
    data = SCData(
        X=np.ones((20, 10)),
        obs=pd.DataFrame({"cell_id": [f"c_{i}" for i in range(20)], "leiden": ["0"] * 10 + ["1"] * 10}),
        var=pd.DataFrame({"gene_name": [f"g_{i}" for i in range(10)]}),
    )
    reg.register(
        uri_str="adata://AD/input/v1",
        payload=data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id="AD",
        created_by_task="task_000",
        operation="init",
    )

    # Task contract strictly forbids 'recluster'
    contract = TaskContract(
        task_id="task_018",
        capability="trajectory_inference",
        method="rogue_trajectory",
        input_artifacts=["adata://AD/input/v1"],
        allowed_operations=["build_neighbor_graph", "infer_trajectory"],
        forbidden_operations=["recluster", "filter_cells", "normalize"],
    )

    # Execution through capability registry
    result = cap_reg.execute_contract(contract, reg)

    # System must flag policy violation and not accept output
    assert result.status == TaskStatus.POLICY_VIOLATION
    assert result.error_type == ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT
    assert "POLICY VIOLATION" in result.error_message
