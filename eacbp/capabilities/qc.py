"""
Dataset Audit and Quality Control (QC) Capabilities.
"""

from typing import Dict, Any, List
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


class DatasetAuditCapability(BaseCapability):
    """Audits biological replication, batch composition, and modality presence."""

    def __init__(self):
        super().__init__(
            capability_name="dataset_audit",
            implementation_id="sc_audit_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        obs = data.obs

        # Audit biological units and conditions
        condition_col = "condition" if "condition" in obs.columns else obs.columns[0]
        mouse_col = "mouse_id" if "mouse_id" in obs.columns else "sample_id" if "sample_id" in obs.columns else None
        batch_col = "batch" if "batch" in obs.columns else None

        n_cells = data.n_obs
        n_genes = data.n_vars
        conditions = obs[condition_col].unique().tolist()
        batches = obs[batch_col].unique().tolist() if batch_col else ["single_batch"]

        donor_counts = {}
        if mouse_col:
            for cond in conditions:
                donor_counts[cond] = int(obs[obs[condition_col] == cond][mouse_col].nunique())
        else:
            for cond in conditions:
                donor_counts[cond] = 1

        min_reps = min(donor_counts.values()) if donor_counts else 1
        sufficient_reps = min_reps >= 2
        batch_effect_possible = len(batches) > 1

        audit_table = pd.DataFrame([{
            "n_cells": n_cells,
            "n_genes": n_genes,
            "conditions": ", ".join(conditions),
            "batches": ", ".join(batches),
            "donor_replicates": str(donor_counts),
            "min_replicates_per_condition": min_reps,
            "biological_replication_sufficient": sufficient_reps,
            "batch_effect_possible": batch_effect_possible,
        }])

        # Register audit table artifact
        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"table://{uri_obj.study_id}/dataset_audit/v1"

        registry.register(
            uri_str=out_uri,
            payload=audit_table,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="dataset_audit",
            parent_uris=[in_uri],
            summary_metrics={
                "n_cells": n_cells,
                "min_replicates": min_reps,
                "batch_effect_possible": batch_effect_possible,
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["audit_metadata", "assess_replication", "assess_batches"],
            metrics={
                "n_cells": n_cells,
                "min_replicates": min_reps,
                "biological_replication_sufficient": sufficient_reps,
                "batch_effect_possible": batch_effect_possible,
            }
        )


class QCCapability(BaseCapability):
    """Filters low-quality cells, checks mitochondrial distribution, and computes QC stats."""

    def __init__(self):
        super().__init__(
            capability_name="qc",
            implementation_id="sc_qc_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.ANNDATA],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        
        min_genes = contract.parameters.get("min_genes", 100)
        max_mito_pct = contract.parameters.get("max_mito_pct", 20.0)

        obs = data.obs
        n_genes_expressed = (data.X > 0).sum(axis=1)
        mito_pct = obs["percent_mito"].values if "percent_mito" in obs.columns else np.zeros(data.n_obs)

        # Quality mask
        valid_mask = (n_genes_expressed >= min_genes) & (mito_pct <= max_mito_pct)
        n_filtered = int((~valid_mask).sum())

        filtered_data = data.subset_obs(valid_mask)
        filtered_data.obs["qc_passed"] = True

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"adata://{uri_obj.study_id}/qc/v1"

        registry.register(
            uri_str=out_uri,
            payload=filtered_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="filter_cells_qc",
            parent_uris=[in_uri],
            parameters={"min_genes": min_genes, "max_mito_pct": max_mito_pct},
            summary_metrics={
                "initial_cells": data.n_obs,
                "retained_cells": filtered_data.n_obs,
                "filtered_cells": n_filtered,
                "retention_rate": float(filtered_data.n_obs / data.n_obs),
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["filter_low_quality_cells", "mitochondrial_filtering"],
            metrics={
                "initial_cells": data.n_obs,
                "retained_cells": filtered_data.n_obs,
                "retention_rate": float(filtered_data.n_obs / data.n_obs),
            }
        )
