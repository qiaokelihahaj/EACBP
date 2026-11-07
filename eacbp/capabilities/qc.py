"""
Dataset Audit and Quality Control (QC) Capabilities.
"""

from typing import Dict, Any, List
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.biological_units import complete_identity, resolve_donor_role
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
        meta, data = registry.get_scdata(in_uri)
        obs = data.obs

        # Audit biological units and conditions
        condition_col = contract.parameters.get("condition_col", "condition")
        donor_role = resolve_donor_role(obs, contract.parameters)
        mouse_col = donor_role.column
        batch_col = contract.parameters.get("batch_col") or ("batch" if "batch" in obs.columns else None)
        for key in ("condition_col", "donor_col", "batch_col"):
            if contract.parameters.get(key) and contract.parameters[key] not in obs:
                raise ValueError(f"Explicit metadata column {key} is absent: {contract.parameters[key]}")

        n_cells = data.n_obs
        n_genes = data.n_vars
        conditions = obs[condition_col].dropna().unique().tolist() if condition_col in obs else []
        batches = obs[batch_col].dropna().unique().tolist() if batch_col else []

        donor_counts = {}
        if mouse_col:
            for cond in conditions:
                donor_counts[cond] = int(obs[obs[condition_col] == cond][mouse_col].nunique())
        else:
            for cond in conditions:
                donor_counts[cond] = 0

        min_reps = min(donor_counts.values()) if donor_counts else 0
        donor_complete = bool(mouse_col and complete_identity(obs[mouse_col]))
        sufficient_reps = min_reps >= 2 and donor_complete
        paired_donors = bool(mouse_col and condition_col in obs and (obs.groupby(mouse_col)[condition_col].nunique() > 1).any())
        batch_effect_possible = len(batches) > 1

        audit_table = pd.DataFrame([{
            "n_cells": n_cells,
            "n_genes": n_genes,
            "conditions": ", ".join(map(str, conditions)) or "unknown",
            "batches": ", ".join(map(str, batches)) or "unknown",
            "donor_replicates": str(donor_counts),
            "min_replicates_per_condition": min_reps,
            "biological_replication_sufficient": sufficient_reps,
            "batch_effect_possible": batch_effect_possible,
            **{**donor_role.as_metrics(), "donor_col": mouse_col or "unavailable"},
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
                **donor_role.as_metrics(),
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
                "donor_replicates": donor_counts,
                "observed_conditions": list(map(str, conditions)),
                "condition_metadata_complete": condition_col in obs and not obs[condition_col].isna().any(),
                "condition_col": condition_col,
                **donor_role.as_metrics(),
                "donor_metadata_complete": donor_complete,
                "paired_donors_observed": paired_donors,
                "total_samples": int(obs[mouse_col].nunique()) if mouse_col else 0,
                "batches": batches,
                "min_replicates": min_reps,
                "biological_replication_sufficient": sufficient_reps,
                "batch_effect_possible": batch_effect_possible,
            }
        )


class QCCapability(BaseCapability):
    """Applies QC thresholds and records dispositions without diagnosing damage."""

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
        meta, data = registry.get_scdata(in_uri)
        
        min_genes = float(contract.parameters.get("min_genes", 100))
        max_mito_pct = float(contract.parameters.get("max_mito_pct", 20.0))
        if not np.isfinite(min_genes) or min_genes < 0:
            raise ValueError("min_genes must be finite and nonnegative")
        if not np.isfinite(max_mito_pct) or not 0 <= max_mito_pct <= 100:
            raise ValueError("max_mito_pct must be between 0 and 100")

        obs = data.obs
        counts = data.layers.get("counts", data.X)
        counts_source = "layers.counts" if "counts" in data.layers else "X"
        values = counts.data if hasattr(counts, "tocsr") else np.asarray(counts)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("QC requires finite nonnegative counts")
        n_genes_expressed = np.asarray((counts > 0).sum(axis=1)).ravel()
        total_counts = np.asarray(counts.sum(axis=1)).ravel()
        n_mito_genes = None
        if "percent_mito" in obs:
            mito_pct = pd.to_numeric(obs["percent_mito"], errors="raise").to_numpy(dtype=float)
            mito_source = "obs.percent_mito"
        else:
            if "mt" in data.var:
                if data.var["mt"].dtype != bool or data.var["mt"].isna().any():
                    raise ValueError("var.mt must contain complete boolean mitochondrial annotations")
                mito_mask = data.var["mt"].to_numpy(dtype=bool)
                mito_source = "var.mt"
            else:
                genes = data.var["gene_name"] if "gene_name" in data.var else data.var.index
                mito_mask = np.asarray(pd.Index(genes).astype(str).str.lower().str.startswith("mt-"))
                mito_source = "gene_symbol_prefix"
            n_mito_genes = int(mito_mask.sum())
            if n_mito_genes:
                mito_counts = np.asarray(counts[:, mito_mask].sum(axis=1)).ravel()
                mito_pct = np.divide(mito_counts * 100.0, total_counts,
                                     out=np.zeros(data.n_obs, dtype=float), where=total_counts > 0)
            else:
                mito_pct = None
                mito_source = "unavailable"
        if mito_pct is not None and (not np.isfinite(mito_pct).all() or np.any((mito_pct < 0) | (mito_pct > 100))):
            raise ValueError("percent_mito must contain finite percentages between 0 and 100")

        # Quality mask
        valid_mask = n_genes_expressed >= min_genes
        if mito_pct is not None:
            valid_mask &= mito_pct <= max_mito_pct
        n_filtered = int((~valid_mask).sum())

        # Preserve every input row's threshold disposition inside the registered
        # QC artifact, including cells excluded from downstream analysis. These
        # flags describe the mask; they are not independent damage diagnoses.
        low_gene_flag = n_genes_expressed < min_genes
        high_mito_flag = mito_pct > max_mito_pct if mito_pct is not None else np.zeros(data.n_obs, dtype=bool)
        reasons = np.full(data.n_obs, "retained", dtype=object)
        reasons[low_gene_flag] = "low_gene_count"
        reasons[high_mito_flag] = "high_mito_percentage"
        reasons[low_gene_flag & high_mito_flag] = "low_gene_count_and_high_mito_percentage"
        reason_names = ("retained", "low_gene_count", "high_mito_percentage", "low_gene_count_and_high_mito_percentage")
        disposition_counts = {reason: int(np.count_nonzero(reasons == reason)) for reason in reason_names}
        mito_available = np.full(data.n_obs, mito_pct is not None, dtype=bool)
        if mito_source in ("var.mt", "gene_symbol_prefix"):
            mito_available &= total_counts > 0
        disposition_mito = np.full(data.n_obs, np.nan, dtype=float)
        if mito_pct is not None:
            disposition_mito[mito_available] = mito_pct[mito_available]
        disposition = {
            "ids": obs.index.astype(str).to_numpy(),
            "n_genes": n_genes_expressed,
            "total_counts": total_counts,
            "mito_pct": disposition_mito,
            "mito_pct_available": mito_available,
            "low_gene_flag": low_gene_flag,
            "high_mito_flag": high_mito_flag,
            "retained": valid_mask.copy(),
            "reason": reasons.astype(str),
        }
        interpretation_warning = (
            "Low RNA content alone does not independently establish cellular damage; "
            "low_gene_flag records a threshold disposition, not a damage label. "
            "Population-dependent RNA content and QC sensitivity remain untested."
        )

        filtered_data = data.subset_obs(valid_mask)
        filtered_data.obs["qc_passed"] = True
        filtered_data.obs["n_genes_by_counts"] = n_genes_expressed[valid_mask]
        filtered_data.obs["total_counts"] = total_counts[valid_mask]
        if mito_pct is not None:
            filtered_data.obs["percent_mito"] = mito_pct[valid_mask]
        metrics = {
            "initial_cells": data.n_obs,
            "retained_cells": filtered_data.n_obs,
            "filtered_cells": n_filtered,
            "retention_rate": float(filtered_data.n_obs / data.n_obs) if data.n_obs else 0.0,
            "counts_source": counts_source,
            "mitochondrial_filtering_applied": mito_pct is not None,
            "mitochondrial_metric_source": mito_source,
            "n_mitochondrial_genes": n_mito_genes,
            "cell_disposition_schema_version": 1,
            "cell_disposition_counts": disposition_counts,
            "interpretation_warning": interpretation_warning,
        }
        filtered_data.uns["qc"] = {
            **metrics, "min_genes": min_genes, "max_mito_pct": max_mito_pct,
            "cell_disposition": disposition,
            "cell_disposition_id_source": "obs.index",
            "input_artifact_uri": in_uri,
            "input_sha256": meta.sha256_hash,
        }

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
            summary_metrics=metrics,
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["filter_low_quality_cells"] + (["mitochondrial_filtering"] if mito_pct is not None else []),
            metrics=metrics,
        )
