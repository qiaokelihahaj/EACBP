"""
Computational Validator auditing matrix integrity, non-finite values (NaN/Inf), and memory/dimensions.
"""

from typing import Dict, Any
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.schemas.artifact import ArtifactType
from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData


class ComputationalValidator(BaseAuditor):
    """Audits computational execution correctness, non-finite values, and matrix dimensions."""

    def __init__(self):
        super().__init__(auditor_name="computational_validator")

    def audit(
        self,
        contract: TaskContract,
        result: TaskResult,
        registry: ArtifactRegistry,
    ) -> ValidationReport:
        target_uri = result.output_artifacts[0] if result.output_artifacts else None
        report = ValidationReport(
            auditor_name=self.auditor_name,
            target_task_id=contract.task_id,
            target_artifact_uri=target_uri,
        )

        if not target_uri or not registry.exists(target_uri):
            report.add_check(
                name="artifact_existence",
                passed=False,
                severity=ValidationSeverity.ERROR,
                message=f"Output artifact '{target_uri}' does not exist in registry.",
            )
            return report

        meta, payload = registry.get(target_uri)

        # 1. Check for AnnData / SCData / SpatialData payloads
        if meta.type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
            data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
            
            # Check shape
            n_cells, n_genes = data.shape
            has_cells = n_cells > 0
            has_genes = n_genes > 0
            report.add_check(
                name="matrix_non_empty",
                passed=has_cells and has_genes,
                severity=ValidationSeverity.ERROR,
                message=f"Matrix dimensions are valid: {n_cells} cells x {n_genes} genes.",
                metrics={"n_cells": n_cells, "n_genes": n_genes},
            )

            # Check NaNs / Infs in X
            if hasattr(data.X, "tocsr"):
                nan_count = int(np.isnan(data.X.data).sum())
                inf_count = int(np.isinf(data.X.data).sum())
            else:
                x_arr = np.asarray(data.X, dtype=np.float32)
                nan_count = int(np.isnan(x_arr).sum())
                inf_count = int(np.isinf(x_arr).sum())
            report.add_check(
                name="expression_finite_values",
                passed=(nan_count == 0 and inf_count == 0),
                severity=ValidationSeverity.ERROR,
                message=f"Expression matrix contains {nan_count} NaNs and {inf_count} Infs.",
                metrics={"nan_count": nan_count, "inf_count": inf_count},
            )

            # Check spatial coordinates dimension match if present
            if "spatial" in data.obsm:
                spatial_coords = np.asarray(data.obsm["spatial"], dtype=np.float32)
                spatial_nan = int((~np.isfinite(spatial_coords)).sum())
                spatial_match = spatial_coords.shape[0] == n_cells
                report.add_check(
                    name="embedding_spatial_finite",
                    passed=(spatial_nan == 0 and spatial_match),
                    severity=ValidationSeverity.ERROR,
                    message=f"Spatial coordinates contain {spatial_nan} NaNs (shape: {spatial_coords.shape}).",
                    metrics={"spatial_nan": spatial_nan, "spatial_shape": list(spatial_coords.shape)},
                )

            # Check embeddings finite values
            for emb_name, emb_val in data.obsm.items():
                if emb_name == "spatial":
                    continue
                try:
                    emb_arr = np.asarray(emb_val, dtype=np.float32)
                    emb_nan = int((~np.isfinite(emb_arr)).sum())
                except Exception:
                    emb_nan = 1

                report.add_check(
                    name=f"embedding_{emb_name}_finite",
                    passed=(emb_nan == 0),
                    severity=ValidationSeverity.ERROR,
                    message=f"Embedding '{emb_name}' contains {emb_nan} NaNs.",
                    metrics={"emb_name": emb_name, "nan_count": emb_nan},
                )

        # 2. Check Table payloads
        elif meta.type == ArtifactType.TABLE:
            df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            is_empty = df.empty
            report.add_check(
                name="table_non_empty",
                passed=not is_empty or (len(df.columns) > 0 and contract.capability in ("trajectory_inference", "deg", "spatial_deg", "knowledge_retrieval", "cell_cell_communication", "liana_communication")),
                severity=ValidationSeverity.ERROR,
                message=f"Table artifact has {len(df)} rows and {len(df.columns)} columns.",
                metrics={"rows": len(df), "cols": len(df.columns)},
            )

            # Check for non-finite values in numeric columns
            if not is_empty:
                numeric_cols = df.select_dtypes(include=[np.number]).columns
                check_df = df
                if contract.capability == "knowledge_retrieval" and "category" in df:
                    # Literature records have no statistical test; their null p/FDR values are intentional.
                    check_df = df[df["category"] != "Local_Curated_Literature"]
                total_nans = int((~np.isfinite(check_df[numeric_cols].to_numpy(dtype=float))).sum())
                from eacbp.auditor.advanced_statistics import METHODS, NULLABLE_STATISTICS
                if result.method_used in METHODS:
                    # These methods retain unestimated statistical quantities.
                    # Their independent validator still rejects infinities,
                    # invalid probabilities, designs and fabricated intervals.
                    numeric = check_df[numeric_cols]
                    total_nans = int(np.isinf(numeric.to_numpy(dtype=float)).sum())
                    required_numeric = [c for c in numeric_cols if c not in NULLABLE_STATISTICS]
                    if result.method_used == "pydeseq2_leave_one_donor_out_v1" and "status" in check_df and check_df.status.eq("skipped").all():
                        required_numeric = [c for c in required_numeric if c not in {"gene", "left_out_donor"}]
                    total_nans += int(numeric[required_numeric].isna().sum().sum())
                if result.method_used == "liana_rank_aggregate_v1" and "comparison_status" in check_df:
                    nullable = [c for c in numeric_cols if c.startswith(("mean_", "comparison_p_value_", "comparison_fdr_"))]
                    total_nans = int(np.isinf(check_df[numeric_cols].to_numpy(dtype=float)).sum())
                    total_nans += int(check_df[[c for c in numeric_cols if c not in nullable]].isna().sum().sum())
                report.add_check(
                    name="table_finite_values",
                    passed=(total_nans == 0),
                    severity=ValidationSeverity.ERROR if total_nans > 0 else ValidationSeverity.INFO,
                    message=f"Table numeric columns contain {total_nans} non-finite values.",
                    metrics={"total_nans": total_nans},
                )

        # 3. Check JSON payloads
        elif meta.type == ArtifactType.JSON:
            has_content = bool(payload)
            report.add_check(
                name="json_valid_payload",
                passed=has_content,
                severity=ValidationSeverity.ERROR,
                message="JSON artifact has valid structured payload." if has_content else "JSON payload is empty.",
                metrics={"has_content": has_content},
            )

        return report
