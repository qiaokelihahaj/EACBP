"""Independent reconstruction of descriptive cluster-marker tables.

No executor helper or marker reference is used here. The audit reconstructs
the full gene universe from the input and checks every emitted rank, including
whether a top-20 table actually contains the top 20 input-derived features.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.schemas.artifact import ArtifactType


class ClusterMarkersValidator(BaseAuditor):
    def __init__(self):
        super().__init__("cluster_markers_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(
            auditor_name=self.auditor_name, target_task_id=contract.task_id,
            target_artifact_uri=result.output_artifacts[0] if result.output_artifacts else None,
        )
        if contract.capability != "cluster_markers":
            return report
        try:
            self._audit_table(contract, result, registry, report)
        except (ValueError, TypeError, KeyError, IndexError, AttributeError) as exc:
            report.add_check("cluster_marker_integrity", False, ValidationSeverity.ERROR,
                             f"Cannot independently reconstruct cluster markers: {exc}")
        return report

    def _audit_table(self, contract, result, registry, report):
        if result.method_used != "cluster_vs_rest_descriptive_v1":
            raise ValueError("Unrecognized cluster-marker method")
        if len(contract.input_artifacts) != 1 or len(result.output_artifacts) != 1:
            raise ValueError("Exactly one single-cell input and one marker table are required")
        data = registry.get_scdata(contract.input_artifacts[0])[1]
        meta, table = registry.get(result.output_artifacts[0])
        if meta.type != ArtifactType.TABLE or not isinstance(table, pd.DataFrame) or table.empty:
            raise ValueError("Cluster-marker output must be a nonempty table")
        cluster_col = contract.parameters.get("cluster_col", "cluster")
        top_n = contract.parameters.get("top_n", 20)
        if top_n is not None and (isinstance(top_n, bool) or int(top_n) != top_n or top_n < 20):
            raise ValueError("top_n must be at least 20, or None for all genes")
        if data.n_vars == 0 or cluster_col not in data.obs or data.obs[cluster_col].isna().any():
            raise ValueError("Input has no gene universe or complete cluster column")
        labels = data.obs[cluster_col].map(str).to_numpy()
        clusters = sorted(np.unique(labels))
        if len(clusters) < 2 or data.obs[cluster_col].nunique() != len(clusters):
            raise ValueError("Input must have at least two unambiguous clusters")
        genes = (data.var["gene_name"] if "gene_name" in data.var else data.var.index).astype(str).to_numpy()
        target_sum = float(contract.parameters.get("target_sum", 10000.0))
        normalization = data.uns.get("normalization", {})
        if "counts" in data.layers:
            expression = data.layers["counts"]
            expression_source, transform = "layers.counts", True
        elif normalization.get("log1p") is True:
            expression = data.layers["normalized"] if "normalized" in data.layers else data.X
            expression_source = "layers.normalized" if "normalized" in data.layers else "X_normalized"
            target_sum = float(normalization.get("target_sum", target_sum))
            transform = False
        else:
            expression, expression_source, transform = data.X, "X_integer_counts", True
            values = expression.data if sparse.issparse(expression) else np.asarray(expression)
            if not (values == np.floor(values)).all():
                raise ValueError("Undocumented non-count expression source")
        values = expression.data if sparse.issparse(expression) else np.asarray(expression)
        if not np.isfinite(values).all() or (values < 0).any() or not np.isfinite(target_sum) or target_sum <= 0:
            raise ValueError("Non-finite or negative expression, or invalid library-size target")
        if sparse.issparse(expression):
            expression = sparse.csr_matrix(expression, dtype=np.float64, copy=True)
            expression.sum_duplicates()
            if transform:
                totals = np.asarray(expression.sum(axis=1)).ravel()
                scale = np.divide(target_sum, totals, out=np.full(len(totals), target_sum), where=totals > 0)
                expression = (sparse.diags(scale) @ expression).tocsr()
                expression.data = np.log1p(expression.data)
        else:
            expression = np.asarray(expression, dtype=np.float64).copy()
            if transform:
                totals = expression.sum(axis=1)
                expression = np.log1p(expression * (target_sum / np.where(totals > 0, totals, 1.0))[:, None])

        numeric = ["gene_index", "rank", "mean_in", "mean_rest", "detection_in", "detection_rest",
                   "mean_log1p_difference", "detection_difference", "n_in", "n_rest", "target_sum"]
        text = ["cluster", "gene", "expression_source", "effect_name", "effect_units", "analysis_scope"]
        if not set(numeric + text).issubset(table.columns):
            raise ValueError("Marker table is missing its descriptive effect/provenance columns")
        converted = table[numeric].apply(pd.to_numeric, errors="coerce")
        scope_ok = bool(
            np.isfinite(converted.to_numpy(dtype=float)).all()
            and table["effect_name"].eq("mean_log1p_expression_difference").all()
            and table["effect_units"].eq("natural_log1p_normalized_expression").all()
            and table["analysis_scope"].eq("within_dataset_descriptive").all()
            and contract.inference_contract is None and result.inference_contract is None and result.scientific_result is None
            and not any(column.lower() in {"p", "p_value", "pvalue", "p_val", "p_val_adj", "pvals", "pvals_adj", "fdr", "fdr_q_value",
                                          "q", "qvalue", "q_value", "padj", "logfc", "log2fc", "logfoldchanges", "log2foldchange"}
                        for column in table.columns)
        )
        report.add_check("cluster_marker_descriptive_scope", scope_ok, ValidationSeverity.ERROR,
                         "Effects must be descriptive mean-log1p differences with no fabricated significance or fold changes.")
        provenance_ok = bool(
            meta.parent_uris == contract.input_artifacts and meta.created_by_task == contract.task_id
            and meta.operation == "summarize_and_rank_existing_cluster_expression"
            and table["expression_source"].eq(expression_source).all()
            and np.allclose(converted["target_sum"], target_sum, rtol=0, atol=1e-12)
            and meta.parameters.get("cluster_col") == cluster_col
            and meta.parameters.get("top_n") == top_n
            and meta.parameters.get("expression_source") == expression_source
            and meta.parameters.get("effective_target_sum") == target_sum
        )
        report.add_check("cluster_marker_input_provenance", provenance_ok, ValidationSeverity.ERROR,
                         "Expression source, target scale, cluster settings and parent lineage must match the audited input.")
        expected_rows = data.n_vars if top_n is None else min(int(top_n), data.n_vars)
        coverage_ok = bool(set(table["cluster"].astype(str)) == set(clusters)
                           and len(table) == len(clusters) * expected_rows
                           and not table.duplicated(["cluster", "gene_index"]).any())
        numbers_ok, ranking_ok = True, True
        for cluster in clusters:
            inside = labels == cluster
            n_in, n_rest = int(inside.sum()), int((~inside).sum())
            # Reconstruct group sums directly rather than reading reported means.
            mean_in = np.asarray(expression[inside].sum(axis=0)).ravel() / n_in
            mean_rest = np.asarray(expression[~inside].sum(axis=0)).ravel() / n_rest
            detection_in = np.asarray((expression[inside] > 0).sum(axis=0)).ravel() / n_in
            detection_rest = np.asarray((expression[~inside] > 0).sum(axis=0)).ravel() / n_rest
            effect, detection = mean_in - mean_rest, detection_in - detection_rest
            expected_order = np.lexsort((np.arange(data.n_vars), genes, -np.round(detection, 12), -np.round(effect, 12)))[:expected_rows]
            rows = table.loc[table["cluster"].astype(str) == cluster].sort_values("rank")
            actual = converted.loc[rows.index]
            coverage_ok = coverage_ok and len(rows) == expected_rows
            if len(rows) != expected_rows:
                numbers_ok = ranking_ok = False
                continue
            ranking_ok = ranking_ok and bool(
                np.array_equal(actual["rank"].to_numpy(), np.arange(1, expected_rows + 1))
                and np.array_equal(actual["gene_index"].to_numpy(), expected_order)
                and np.array_equal(rows["gene"].astype(str).to_numpy(), genes[expected_order])
            )
            references = {
                "mean_in": mean_in[expected_order], "mean_rest": mean_rest[expected_order],
                "detection_in": detection_in[expected_order], "detection_rest": detection_rest[expected_order],
                "mean_log1p_difference": effect[expected_order], "detection_difference": detection[expected_order],
                "n_in": np.repeat(n_in, expected_rows), "n_rest": np.repeat(n_rest, expected_rows),
            }
            numbers_ok = numbers_ok and all(np.allclose(actual[column], values, rtol=1e-7, atol=1e-9)
                                           for column, values in references.items())
        report.add_check("cluster_marker_cluster_coverage", bool(coverage_ok), ValidationSeverity.ERROR,
                         "Every observed cluster must have the complete requested ranking coverage.",
                         {"n_clusters": len(clusters), "rows_per_cluster": expected_rows, "n_genes_evaluated": data.n_vars})
        report.add_check("cluster_marker_values_and_ranks", bool(numbers_ok and ranking_ok), ValidationSeverity.ERROR,
                         "Recomputed mean/detection effects and full-universe top-gene ordering from the input matrix.")
        summary_ok = all(result.metrics.get(key) == value and meta.summary_metrics.get(key) == value for key, value in {
            "n_cells": data.n_obs, "n_clusters": len(clusters), "n_genes_evaluated": data.n_vars,
            "n_marker_rows": len(table), "expression_source": expression_source,
            "effect_name": "mean_log1p_expression_difference", "effect_units": "natural_log1p_normalized_expression",
            "analysis_scope": "within_dataset_descriptive", "confirmatory_allowed": False,
        }.items())
        report.add_check("cluster_marker_integrity", bool(scope_ok and provenance_ok and coverage_ok
                                                        and numbers_ok and ranking_ok and summary_ok), ValidationSeverity.ERROR,
                         "Independent cluster-marker reconstruction passed; interpretation remains within-dataset and descriptive.")
