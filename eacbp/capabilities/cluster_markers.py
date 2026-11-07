"""Reference-free, descriptive markers for existing single-cell clusters.

The effect is the difference of mean log1p library-size-normalized expression,
not a fold change or an independent statistical test. Ranking the same cells
used for clustering is exploratory; neither donor DE nor cell-type identity is
established by this table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field
from scipy import sparse

from eacbp.artifact.uri import ArtifactURI
from eacbp.capabilities.base import BaseCapability, CapabilityDescriptor
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceStrength, EvidenceType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


class ClusterMarkerParameters(BaseModel):
    model_config = ConfigDict(extra="allow")

    cluster_col: str = Field("cluster", min_length=1)
    # Keep at least the standard top-20 audit scope, or emit the full ranking.
    top_n: int | None = Field(20, ge=20)
    target_sum: float = Field(10000.0, gt=0, allow_inf_nan=False)


def _expression_for_summary(data, target_sum):
    """Choose a proven expression source and preserve sparse representation."""
    normalization = data.uns.get("normalization", {})
    if "counts" in data.layers:
        matrix, source, counts = data.layers["counts"], "layers.counts", True
    elif normalization.get("log1p") is True:
        matrix = data.layers.get("normalized", data.X)
        source = "layers.normalized" if "normalized" in data.layers else "X_normalized"
        target_sum = float(normalization.get("target_sum", target_sum))
        counts = False
    else:
        matrix, source, counts = data.X, "X_integer_counts", True
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
        if not np.equal(values, np.rint(values)).all():
            raise ValueError("Cluster markers require counts or explicitly documented log1p-normalized expression")
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Cluster markers require finite nonnegative expression")
    if not np.isfinite(target_sum) or target_sum <= 0:
        raise ValueError("Expression target_sum must be finite and positive")
    matrix = matrix.astype(np.float64, copy=True)
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
        matrix.sum_duplicates()
    if counts:
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        factors = float(target_sum) / np.where(totals > 0, totals, 1.0)
        if sparse.issparse(matrix):
            matrix = matrix.multiply(factors[:, None]).tocsr()
            matrix.data = np.log1p(matrix.data)
        else:
            matrix = np.log1p(matrix * factors[:, None])
    return matrix, source, float(target_sum)


def _plan_cluster_markers(*, manifest, parameters, tasks):
    if any(task.capability == "cluster_markers" for task in tasks):
        return []
    clustering = next((task for task in reversed(tasks) if task.capability == "clustering"), None)
    if clustering is None:
        raise ValueError("Cluster-marker discovery requires an existing clustering task")
    return TaskContract(
        task_id="task_005b_cluster_markers", capability="cluster_markers",
        method="cluster_vs_rest_descriptive_v1", depends_on=[clustering.task_id],
        input_artifacts=[clustering.expected_outputs[0]], parameters=parameters,
        allowed_operations=["summarize_cluster_expression", "rank_cluster_markers"],
        forbidden_operations=["filter_cells", "recluster", "replace_cell_type_labels"],
        expected_outputs=[f"table://{manifest.study_id}/cluster_markers/v1"],
        validation_requirements=["cluster_marker_integrity"],
    )


def _extract_marker_evidence(contract, result, report, registry):
    table = registry.get(result.output_artifacts[0])[1]
    input_data = registry.get_scdata(contract.input_artifacts[0])[1]
    top_markers = {}
    for cluster, rows in table.groupby("cluster", sort=True, observed=True):
        top_markers[str(cluster)] = rows.sort_values("rank").head(20)[[
            "gene", "gene_index", "rank", "mean_log1p_difference", "detection_difference"
        ]].to_dict(orient="records")
    return [EvidenceNode(
        evidence_id=f"E_{contract.task_id}_cluster_markers", type=EvidenceType.CLUSTER_MARKERS,
        strength=EvidenceStrength.WEAK, score=0.5,
        summary=(f"Descriptive cluster-versus-rest expression rankings cover {len(top_markers)} clusters. "
                 "Effects are differences of mean log1p normalized expression; rankings do not establish "
                 "cell-type identity, donor-level differential expression, or independent replication."),
        source_task_id=contract.task_id, source_artifact_uris=list(result.output_artifacts),
        audit_passed=True, is_simulated=bool(input_data.uns.get("is_simulated", False)),
        validation_scope="within_dataset_descriptive", data_origin_uris=list(contract.input_artifacts),
        metrics={**result.metrics, "top_markers_by_cluster": top_markers,
                 "independent_replication": False, "confirmatory_allowed": False},
    )]


class ClusterMarkersCapability(BaseCapability):
    def __init__(self):
        from eacbp.auditor.cluster_markers import ClusterMarkersValidator

        super().__init__(
            "cluster_markers", "cluster_vs_rest_descriptive_v1",
            accepts_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
            output_types=[ArtifactType.TABLE],
            descriptor=CapabilityDescriptor(
                "cluster_markers", "cluster_vs_rest_descriptive_v1",
                parameter_model=ClusterMarkerParameters,
                input_types=[ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA],
                output_types=[ArtifactType.TABLE], required_audit_ids=["cluster_marker_integrity"],
                validator=ClusterMarkersValidator(), evidence_extractor=_extract_marker_evidence,
                plan_factory=_plan_cluster_markers, scope="shared",
                description="Reference-free descriptive gene-expression and detection rankings for existing clusters.",
            ),
        )

    def execute(self, contract, registry):
        parameters = ClusterMarkerParameters.model_validate(contract.parameters)
        if len(contract.input_artifacts) != 1:
            raise ValueError("Cluster markers require exactly one annotated single-cell input")
        in_uri = contract.input_artifacts[0]
        data = registry.get_scdata(in_uri)[1]
        if data.n_obs == 0 or data.n_vars == 0:
            raise ValueError("Cluster-marker input must contain cells and genes")
        if parameters.cluster_col not in data.obs or data.obs[parameters.cluster_col].isna().any():
            raise ValueError("Cluster-marker input requires a complete cluster column")
        labels = data.obs[parameters.cluster_col].astype(str).to_numpy()
        clusters = sorted(set(labels))
        if len(clusters) < 2:
            raise ValueError("Cluster-versus-rest markers require at least two clusters")
        if data.obs[parameters.cluster_col].nunique() != len(clusters):
            raise ValueError("Cluster labels collide after conversion to strings")
        genes = (data.var["gene_name"] if "gene_name" in data.var else data.var.index).astype(str).to_numpy()
        matrix, source, target_sum = _expression_for_summary(data, parameters.target_sum)
        rows = []
        for cluster in clusters:
            inside = labels == cluster
            n_in, n_rest = int(inside.sum()), int((~inside).sum())
            mean_in = np.asarray(matrix[inside].mean(axis=0)).ravel()
            mean_rest = np.asarray(matrix[~inside].mean(axis=0)).ravel()
            detection_in = np.asarray((matrix[inside] > 0).sum(axis=0)).ravel() / n_in
            detection_rest = np.asarray((matrix[~inside] > 0).sum(axis=0)).ravel() / n_rest
            effects, detection_effects = mean_in - mean_rest, detection_in - detection_rest
            # Quantize only ranking keys to make mathematical ties stable
            # across dense/sparse reduction order; retain unrounded effects.
            order = np.lexsort((np.arange(data.n_vars), genes, -np.round(detection_effects, 12), -np.round(effects, 12)))
            if parameters.top_n is not None:
                order = order[:parameters.top_n]
            for rank, gene_index in enumerate(order, start=1):
                rows.append({
                    "cluster": cluster, "gene": genes[gene_index], "gene_index": int(gene_index), "rank": rank,
                    "mean_in": float(mean_in[gene_index]), "mean_rest": float(mean_rest[gene_index]),
                    "detection_in": float(detection_in[gene_index]), "detection_rest": float(detection_rest[gene_index]),
                    "mean_log1p_difference": float(effects[gene_index]),
                    "detection_difference": float(detection_effects[gene_index]), "n_in": n_in, "n_rest": n_rest,
                    "expression_source": source, "target_sum": target_sum,
                    "effect_name": "mean_log1p_expression_difference", "effect_units": "natural_log1p_normalized_expression",
                    "analysis_scope": "within_dataset_descriptive",
                })
        table = pd.DataFrame(rows)
        study_id = ArtifactURI.parse(in_uri).study_id
        out_uri = f"table://{study_id}/cluster_markers/v1"
        metrics = {"n_cells": data.n_obs, "n_genes_evaluated": data.n_vars,
                   "n_clusters": len(clusters), "n_marker_rows": len(table), "expression_source": source,
                   "effect_name": "mean_log1p_expression_difference", "effect_units": "natural_log1p_normalized_expression",
                   "analysis_scope": "within_dataset_descriptive", "confirmatory_allowed": False}
        registry.register(
            out_uri, table, ArtifactType.TABLE, study_id, contract.task_id,
            "summarize_and_rank_existing_cluster_expression", parent_uris=[in_uri],
            parameters={**parameters.model_dump(), "expression_source": source, "effective_target_sum": target_sum},
            summary_metrics=metrics,
        )
        return TaskResult(
            task_id=contract.task_id, capability=self.capability_name, method_used=self.implementation_id,
            status=TaskStatus.SUCCESS, input_artifacts=[in_uri], output_artifacts=[out_uri],
            executed_operations=["summarize_cluster_expression", "rank_cluster_markers"], metrics=metrics,
        )
