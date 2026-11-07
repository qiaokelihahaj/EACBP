"""Preserve excluded raw-count cells for descriptive, audited QC review.

These outputs are a review branch. They never restore cells to the primary
analysis or establish cell damage, cell-type identity or differential expression.
"""

from __future__ import annotations

from copy import deepcopy
import json

import numpy as np
import pandas as pd
from pydantic import Field
from scipy import sparse

from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.capabilities.parameters import MethodParameters, NonNegativeInt
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


REASON_NAMES = (
    "retained", "low_gene_count", "high_mito_percentage",
    "low_gene_count_and_high_mito_percentage",
)
INTERPRETATION_WARNING = (
    "Low RNA content alone does not independently establish cellular damage; "
    "QC flags describe threshold dispositions, not damage labels. "
    "Raw-count detection summaries do not establish differential expression "
    "or cell-type identity. Candidate cells remain excluded from the primary "
    "analysis; no restoration was performed and no reference labels were used."
)
CELL_COLUMNS = (
    "cell_id", "n_genes", "total_counts", "mito_pct", "mito_pct_available",
    "low_gene_flag", "high_mito_flag", "retained", "reason",
)
FEATURE_COLUMNS = (
    "reason", "gene_index", "gene_id", "gene", "n_cells", "n_detected",
    "total_counts", "detection_fraction", "mean_raw_count", "status",
)


class QCReviewParameters(MethodParameters):
    min_genes: NonNegativeInt = 100
    max_mito_pct: float = Field(20.0, ge=0, le=100)


def _output_uris(study_id):
    return {
        "summary": f"json://{study_id}/qc_review/v1",
        "candidates": f"adata://{study_id}/qc_candidates/v1",
        "cells": f"table://{study_id}/qc_review_cells/v1",
        "features": f"table://{study_id}/qc_review_features/v1",
    }


def _plan_qc_review(*, manifest, parameters, tasks):
    if any(task.capability == "qc_review" for task in tasks):
        return []
    qc = next((task for task in tasks if task.capability == "qc"), None)
    if qc is None or len(qc.input_artifacts) != 1:
        raise ValueError("QC review requires an existing QC task with one raw input")
    if qc.method != "sc_qc_v1":
        raise ValueError("QC review requires the native sc_qc_v1 threshold method")
    supplied = QCReviewParameters.model_validate(parameters)
    inherited = QCReviewParameters.model_validate({
        "min_genes": qc.parameters.get("min_genes", 100),
        "max_mito_pct": qc.parameters.get("max_mito_pct", 20.0),
    })
    for key in ("min_genes", "max_mito_pct"):
        if key in parameters and getattr(supplied, key) != getattr(inherited, key):
            raise ValueError("QC review thresholds must match the primary QC task")
    values = {**parameters, "min_genes": inherited.min_genes,
              "max_mito_pct": inherited.max_mito_pct}
    return TaskContract(
        task_id="task_001b_qc_review", capability="qc_review",
        method="raw_count_qc_review_v1", depends_on=list(qc.depends_on),
        input_artifacts=list(qc.input_artifacts), parameters=values,
        allowed_operations=["summarize_qc_candidates"],
        forbidden_operations=["restore_cells", "normalize", "recluster", "infer_damage",
                              "replace_cell_type_labels"],
        expected_outputs=list(_output_uris(manifest.study_id).values()),
        validation_requirements=["qc_review_integrity"],
    )


def _extract_qc_review_evidence(contract, result, report, registry):
    """The admitted review is reported directly and creates no scientific claim."""
    return []


class QCReviewCapability(BaseCapability):
    def __init__(self):
        from eacbp.auditor.qc_review import QCReviewValidator

        super().__init__(
            "qc_review", "raw_count_qc_review_v1",
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            output_types=[ArtifactType.JSON, ArtifactType.ANNDATA, ArtifactType.TABLE],
            descriptor=CapabilityDescriptor(
                "qc_review", "raw_count_qc_review_v1",
                parameter_model=QCReviewParameters,
                input_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
                output_types=[ArtifactType.JSON, ArtifactType.ANNDATA, ArtifactType.TABLE],
                required_audit_ids=["qc_review_integrity"], validator=QCReviewValidator(),
                evidence_extractor=_extract_qc_review_evidence,
                plan_factory=_plan_qc_review, scope="shared",
                description="Preserve excluded raw counts and descriptive QC dispositions without restoring cells.",
            ),
        )

    def execute(self, contract, registry):
        parameters = QCReviewParameters.model_validate(contract.parameters)
        if len(contract.input_artifacts) != 1:
            raise ValueError("QC review requires exactly one raw single-cell input")
        in_uri = contract.input_artifacts[0]
        source_meta, data = registry.get_scdata(in_uri)
        if data.n_obs == 0 or data.n_vars == 0:
            raise ValueError("QC review input must contain cells and genes")
        ids = data.obs.index.astype(str).to_numpy()
        if not pd.Index(ids).is_unique:
            raise ValueError("QC review cell identifiers collide after conversion to strings")
        counts = data.layers["counts"] if "counts" in data.layers else data.X
        counts_source = "layers.counts" if "counts" in data.layers else "X"
        values = counts.data if sparse.issparse(counts) else np.asarray(counts)
        if (not np.isfinite(values).all() or np.any(values < 0)
                or not np.equal(values, np.rint(values)).all()):
            raise ValueError("QC review requires finite nonnegative integer raw counts")
        if np.issubdtype(values.dtype, np.integer) and values.size:
            # NumPy/SciPy integer sums can silently wrap before a finite-value
            # check sees them. Bound both row and feature reductions first.
            limit = (np.iinfo(np.uint64).max if values.dtype == np.uint64
                     else np.iinfo(np.int64).max)
            if int(values.max()) * max(data.n_obs, data.n_vars) > int(limit):
                raise ValueError("QC review raw-count reductions could overflow")
        n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if not np.isfinite(totals).all() or np.any(totals < 0):
            raise ValueError("QC review raw-count totals are invalid")
        mito_pct, mito_source, n_mito = None, "unavailable", None
        if "percent_mito" in data.obs:
            mito_pct = pd.to_numeric(data.obs["percent_mito"], errors="raise").to_numpy(dtype=float)
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
            n_mito = int(mito_mask.sum())
            if n_mito:
                mito_counts = np.asarray(counts[:, mito_mask].sum(axis=1)).ravel()
                mito_pct = np.divide(mito_counts * 100.0, totals,
                                     out=np.zeros(data.n_obs, dtype=float), where=totals > 0)
            else:
                mito_source = "unavailable"
        if mito_pct is not None and (not np.isfinite(mito_pct).all()
                or np.any((mito_pct < 0) | (mito_pct > 100))):
            raise ValueError("percent_mito must contain finite percentages between 0 and 100")
        low = n_genes < parameters.min_genes
        high = mito_pct > parameters.max_mito_pct if mito_pct is not None else np.zeros(data.n_obs, dtype=bool)
        retained = ~(low | high)
        reasons = np.full(data.n_obs, "retained", dtype=object)
        reasons[low] = "low_gene_count"
        reasons[high] = "high_mito_percentage"
        reasons[low & high] = "low_gene_count_and_high_mito_percentage"
        available = np.full(data.n_obs, mito_pct is not None, dtype=bool)
        if mito_source in ("var.mt", "gene_symbol_prefix"):
            available &= totals > 0
        recorded_mito = np.full(data.n_obs, np.nan, dtype=float)
        if mito_pct is not None:
            recorded_mito[available] = mito_pct[available]
        cells = pd.DataFrame({
            "cell_id": ids, "n_genes": n_genes, "total_counts": totals,
            "mito_pct": recorded_mito, "mito_pct_available": available,
            "low_gene_flag": low, "high_mito_flag": high,
            "retained": retained, "reason": reasons,
        }, columns=CELL_COLUMNS)
        outputs = _output_uris(ArtifactURI.parse(in_uri).study_id)
        summary = {
            "schema_version": 1, "source_artifact_uri": in_uri,
            "source_sha256": source_meta.sha256_hash, "counts_source": counts_source,
            "min_genes": parameters.min_genes, "max_mito_pct": parameters.max_mito_pct,
            "initial_cells": data.n_obs, "candidate_cells": int((~retained).sum()),
            "retained_cells": int(retained.sum()), "n_genes": data.n_vars,
            "mitochondrial_metric_source": mito_source,
            "mitochondrial_filtering_applied": mito_pct is not None,
            "n_mitochondrial_genes": n_mito,
            "cell_disposition_counts": {reason: int((reasons == reason).sum()) for reason in REASON_NAMES},
            "zero_count_cells": int((totals == 0).sum()), "outputs": outputs,
            "restoration_performed": False, "damage_diagnosis": False,
            "reference_labels_used": False, "analysis_scope": "within_dataset_descriptive",
            "interpretation_warning": INTERPRETATION_WARNING,
        }
        # Fail before registering outputs if a future summary adds NaN/Infinity.
        json.dumps(summary, allow_nan=False)
        selected_counts = counts[~retained].copy()
        candidate_obs = data.obs.iloc[np.flatnonzero(~retained)].copy(deep=True)
        candidate_obs["qc_review_reason"] = reasons[~retained]
        candidate_obs["qc_review_candidate"] = True
        candidate_obs["qc_passed"] = False
        candidates = SCData(
            sparse.csr_matrix(selected_counts, copy=True), candidate_obs,
            data.var.copy(deep=True), layers={"counts": selected_counts},
            uns={"qc_review": deepcopy(summary)},
        )
        gene_ids = data.var.index.astype(str).to_numpy()
        gene_names = (data.var["gene_name"] if "gene_name" in data.var else data.var.index).astype(str).to_numpy()
        feature_rows = []
        for reason in REASON_NAMES:
            selected = counts[reasons == reason]
            n_cells = selected.shape[0]
            detected = np.asarray((selected > 0).sum(axis=0)).ravel()
            raw_totals = np.asarray(selected.sum(axis=0)).ravel()
            for gene_index in range(data.n_vars):
                feature_rows.append({
                    "reason": reason, "gene_index": gene_index,
                    "gene_id": gene_ids[gene_index], "gene": gene_names[gene_index],
                    "n_cells": n_cells, "n_detected": int(detected[gene_index]),
                    "total_counts": raw_totals[gene_index].item(),
                    "detection_fraction": float(detected[gene_index] / n_cells) if n_cells else np.nan,
                    "mean_raw_count": float(raw_totals[gene_index] / n_cells) if n_cells else np.nan,
                    "status": "described" if n_cells else "no_cells",
                })
        features = pd.DataFrame(feature_rows, columns=FEATURE_COLUMNS)
        finite_columns = ["total_counts", "n_detected", "n_cells"]
        if not np.isfinite(features[finite_columns].to_numpy(dtype=float)).all():
            raise ValueError("QC review feature raw-count totals are invalid")
        artifact_parameters = {"min_genes": parameters.min_genes,
                               "max_mito_pct": parameters.max_mito_pct}
        for key, payload, artifact_type, operation in (
            ("summary", summary, ArtifactType.JSON, "raw_count_qc_review_summary"),
            ("candidates", candidates, ArtifactType.ANNDATA, "preserve_qc_candidates"),
            ("cells", cells, ArtifactType.TABLE, "qc_review_cell_dispositions"),
            ("features", features, ArtifactType.TABLE, "qc_review_feature_counts"),
        ):
            registry.register(
                outputs[key], payload, artifact_type, ArtifactURI.parse(in_uri).study_id,
                contract.task_id, operation, parent_uris=[in_uri],
                parameters=artifact_parameters, summary_metrics=summary,
            )
        return TaskResult(
            task_id=contract.task_id, capability=self.capability_name,
            method_used=self.implementation_id, status=TaskStatus.SUCCESS,
            input_artifacts=[in_uri], output_artifacts=list(outputs.values()),
            executed_operations=["summarize_qc_candidates"], metrics=summary,
        )
