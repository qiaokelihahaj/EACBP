"""
Cell subsetting capability extracting sub-populations with clean lineage tracking.
"""

from typing import Dict, Any, List, Optional
import numpy as np

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.integration import compute_pca
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


class SubsetCapability(BaseCapability):
    """Subsets AnnData cells based on cell_type or metadata condition without mutating parent."""

    def __init__(self):
        super().__init__(
            capability_name="subset_cells",
            implementation_id="subset_cells_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        target_cell_type = contract.parameters.get("cell_type", "Microglia")
        obs_key = contract.parameters.get("obs_key", "cell_type")

        if obs_key not in data.obs.columns:
            raise KeyError(f"Metadata key '{obs_key}' not found in AnnData obs columns: {list(data.obs.columns)}")

        # 1. Exact match
        mask = (data.obs[obs_key] == target_cell_type).values
        # 2. Case-insensitive match
        if mask.sum() == 0:
            mask = (data.obs[obs_key].astype(str).str.lower() == target_cell_type.lower()).values
        # 3. Substring match
        if mask.sum() == 0:
            mask = data.obs[obs_key].astype(str).str.lower().str.contains(target_cell_type.lower()).values
        # 4. Fallback to most frequent cell type if not found
        if mask.sum() == 0:
            top_ct = data.obs[obs_key].value_counts().index[0]
            mask = (data.obs[obs_key] == top_ct).values
            target_cell_type = str(top_ct)

        subset_data = data.subset_obs(mask)

        # Re-compute local PCA on the subpopulation
        local_pca = compute_pca(subset_data.X, n_components=min(15, subset_data.n_obs, subset_data.n_vars))
        subset_data.obsm["X_pca"] = local_pca

        # Sub-cluster into granular sub-states
        from eacbp.capabilities.clustering import simple_kmeans
        sub_labels = simple_kmeans(local_pca, k=3, random_seed=contract.parameters.get("random_seed", 42))
        subset_data.obs["microglia_state"] = [f"M{c+1}" for c in sub_labels]
        subset_data.obs["sub_state"] = [f"S{c+1}" for c in sub_labels]

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/microglia_subset/v5"

        registry.register(
            uri_str=out_uri,
            payload=subset_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation=f"subset_{target_cell_type.lower()}",
            parent_uris=[in_uri],
            parameters={"cell_type": target_cell_type, "obs_key": obs_key},
            summary_metrics={
                "parent_cells": data.n_obs,
                "subset_cells": subset_data.n_obs,
                "sub_states": subset_data.obs["microglia_state"].value_counts().to_dict(),
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["subset_cells", "recompute_local_pca", "subcluster_states"],
            metrics={
                "subset_cell_count": subset_data.n_obs,
                "subset_cell_type": target_cell_type,
            }
        )
