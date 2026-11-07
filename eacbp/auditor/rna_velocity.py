"""Audit native RNA velocity without treating model output as observed fate.

The auditor compares input counts, identities, neighborhoods and per-cell
summary with persisted model output. It validates kinetic matrices and graph
dimensions but does not independently refit scVelo or verify the biological
kinetic assumptions; those limits remain explicit in evidence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType


def _same_matrix(left, right):
    if left.shape != right.shape:
        return False
    if sparse.issparse(left) or sparse.issparse(right):
        delta = sparse.csr_matrix(left) - sparse.csr_matrix(right)
        delta.eliminate_zeros()
        return delta.nnz == 0
    return np.array_equal(left, right)


def _finite(matrix):
    return np.isfinite(matrix.data if sparse.issparse(matrix) else np.asarray(matrix)).all()


def _counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    return bool(np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all())


def _close_in_gene_blocks(left, right):
    """Bound dense audit memory to 128 genes, even for a complete gene space."""
    if left.shape != right.shape:
        return False
    for start in range(0, left.shape[1], 128):
        a, b = left[:, start:start + 128], right[:, start:start + 128]
        a = a.toarray() if sparse.issparse(a) else np.asarray(a)
        b = b.toarray() if sparse.issparse(b) else np.asarray(b)
        if not np.allclose(a, b, rtol=2e-5, atol=1e-4):
            return False
    return True


class RNAVelocityValidator(BaseAuditor):
    def __init__(self):
        super().__init__("rna_velocity_validator")

    def audit(self, contract, result, registry):
        report = ValidationReport(
            auditor_name=self.auditor_name, target_task_id=contract.task_id,
            target_artifact_uri=result.output_artifacts[0] if result.output_artifacts else None,
        )
        if contract.capability != "rna_velocity":
            return report
        try:
            self._audit(contract, result, registry, report)
        except (ValueError, TypeError, KeyError, IndexError, AttributeError) as exc:
            report.add_check("rna_velocity_integrity", False, ValidationSeverity.ERROR,
                             f"Cannot audit persisted native RNA velocity: {exc}")
        return report

    def _audit(self, contract, result, registry, report):
        if result.method_used != "scvelo_stochastic_v1" or len(contract.input_artifacts) != 1 or len(result.output_artifacts) != 2:
            raise ValueError("Native stochastic RNA velocity requires one input and paired AnnData/table outputs")
        original = registry.get_scdata(contract.input_artifacts[0])[1]
        out_meta, out_payload = registry.get(result.output_artifacts[0])
        table_meta, table = registry.get(result.output_artifacts[1])
        if out_meta.type != ArtifactType.ANNDATA or table_meta.type != ArtifactType.TABLE or not isinstance(table, pd.DataFrame):
            raise ValueError("RNA velocity outputs must be AnnData plus a per-cell summary table")
        data = SCData.from_payload(out_payload)
        ids = original.obs["cell_id"].astype(str).to_numpy() if "cell_id" in original.obs else original.obs.index.astype(str).to_numpy()
        identity_ok = bool(
            original.shape == data.shape and np.array_equal(original.obs.index.astype(str), data.obs.index.astype(str))
            and np.array_equal(original.var.index.astype(str), data.var.index.astype(str)) and len(set(ids)) == original.n_obs
            and original.obs["cluster"].astype(str).equals(data.obs["cluster"].astype(str))
            and _same_matrix(original.X, data.X)
            and _same_matrix(original.obsm["X_pca"], data.obsm["X_pca"])
            and all(_same_matrix(original.obsp[name], data.obsp[name]) for name in ("connectivities", "distances"))
        )
        counts_ok = all(
            name in original.layers and saved in data.layers and _counts(original.layers[name])
            and _counts(data.layers[saved]) and _same_matrix(original.layers[name], data.layers[saved])
            for name, saved in (("counts", "counts"), ("spliced", "spliced_counts"), ("unspliced", "unspliced_counts"))
        )
        if "gene_name" in original.var:
            identity_ok = identity_ok and "gene_name" in data.var and original.var["gene_name"].equals(data.var["gene_name"])
        report.add_check("rna_velocity_count_and_identity_preservation", bool(identity_ok and counts_ok), ValidationSeverity.ERROR,
                         "Compared every raw count, cell/gene identity, project cluster, PCA and neighbor matrix against the input.")
        layers_ok = bool(
            all(name in data.layers and data.layers[name].shape == data.shape and _finite(data.layers[name])
                for name in ("spliced", "unspliced", "Ms", "Mu", "velocity"))
            and "velocity_genes" in data.var and data.var["velocity_genes"].isin([True, False]).all()
            and int(data.var["velocity_genes"].sum()) > 0
        )
        if layers_ok:
            layers_ok = all(((data.layers[name].data if sparse.issparse(data.layers[name]) else np.asarray(data.layers[name])) >= 0).all()
                            for name in ("spliced", "unspliced", "Ms", "Mu"))
        report.add_check("rna_velocity_kinetic_layer_integrity", bool(layers_ok), ValidationSeverity.ERROR,
                         "Stochastic RNA moments, signed velocity, and usable kinetic genes must retain the full input dimensions and finite values.")
        moments_ok = bool(layers_ok and counts_ok)
        if moments_ok:
            target_sum = float(contract.parameters.get("target_sum", 10000.0))
            neighbors = (sparse.csr_matrix(original.obsp["connectivities"]) > 0).astype(np.float32)
            neighbors.setdiag(1)
            neighbors = neighbors.multiply(1.0 / np.asarray(neighbors.sum(axis=1))).tocsr().astype(np.float32)
            for source, moment in (("spliced", "Ms"), ("unspliced", "Mu")):
                raw = sparse.csr_matrix(original.layers[source], dtype=np.float64)
                totals = np.asarray(raw.sum(axis=1)).ravel()
                normalized = sparse.diags(target_sum / np.where(totals > 0, totals, 1.0)) @ raw
                moments_ok = moments_ok and _close_in_gene_blocks(normalized, data.layers[source])
                # scVelo 0.3.x uses positive graph support including self,
                # row-normalized in float32, for first RNA moments.
                for start in range(0, original.n_vars, 128):
                    expected = (neighbors @ sparse.csr_matrix(data.layers[source][:, start:start + 128])).astype(np.float32)
                    moments_ok = moments_ok and _close_in_gene_blocks(expected, data.layers[moment][:, start:start + 128])
        report.add_check("rna_velocity_normalization_and_moments", bool(moments_ok), ValidationSeverity.ERROR,
                         "Independently reconstructed library-size RNA normalization and project-neighbor first moments in bounded gene blocks.")
        fit_parameters = data.uns.get("velocity_params", {})
        residual_ok = bool(layers_ok and fit_parameters.get("mode") == "stochastic"
                           and fit_parameters.get("fit_offset") is False)
        residual_evaluated = False
        if residual_ok:
            selected = data.var["velocity_genes"].to_numpy(dtype=bool)
            for name in ("velocity_r2", "velocity_gamma", "velocity_beta", "velocity_offset", "velocity_offset2"):
                if name in data.var:
                    residual_ok = residual_ok and bool(np.isfinite(pd.to_numeric(data.var.loc[selected, name], errors="coerce")).all())
            if "velocity_gamma" in data.var:
                # Preserve the stored coefficient dtype so cancellation in
                # Mu - gamma*Ms follows the actual scVelo arithmetic.
                gamma = pd.to_numeric(data.var["velocity_gamma"], errors="coerce").to_numpy()
                residual_evaluated = True
                for start in range(0, data.n_vars, 128):
                    ms, mu = data.layers["Ms"][:, start:start + 128], data.layers["Mu"][:, start:start + 128]
                    ms = ms.toarray() if sparse.issparse(ms) else np.asarray(ms)
                    mu = mu.toarray() if sparse.issparse(mu) else np.asarray(mu)
                    expected = mu - ms * gamma[start:start + 128]
                    residual_ok = residual_ok and _close_in_gene_blocks(expected, data.layers["velocity"][:, start:start + 128])
        report.add_check("rna_velocity_parameter_and_residual_consistency", bool(residual_ok), ValidationSeverity.ERROR,
                         ("Checked finite retained kinetic parameters and reconstructed signed velocity as Mu - gamma * Ms without refitting."
                          if residual_evaluated else
                          "Checked stochastic mode and available finite kinetic parameters; scVelo omits constant gamma columns, so the residual equation was not evaluated."),
                         {"residual_equation_evaluated": residual_evaluated, "kinetic_fit_refitted": False})
        graph_ok = True
        graph_edges = 0
        for name, lower, upper in (("velocity_graph", 0.0, 1.0), ("velocity_graph_neg", -1.0, 0.0)):
            graph = data.uns.get(name)
            if graph is None or not hasattr(graph, "shape") or graph.shape != (data.n_obs, data.n_obs):
                graph_ok = False
                continue
            values = graph.data if sparse.issparse(graph) else np.asarray(graph)
            graph_ok = graph_ok and bool(np.isfinite(values).all() and (values >= lower - 1e-6).all() and (values <= upper + 1e-6).all())
            if name == "velocity_graph":
                graph_edges = int(np.count_nonzero(values))
        graph_ok = graph_ok and graph_edges > 0
        report.add_check("rna_velocity_graph_integrity", bool(graph_ok), ValidationSeverity.ERROR,
                         "Validated positive/negative signed kinetic neighbor cosine graphs and nonempty positive edges.",
                         {"positive_graph_edges": graph_edges})
        columns = ["cell_id", "cluster", "velocity_pseudotime", "root_score", "end_score", "velocity_confidence",
                   "analysis_scope", "direction_source", "time_units"]
        if not set(columns).issubset(table.columns):
            raise ValueError("Per-cell summary is missing its model-derived quantities or scope")
        summary_ok = bool(
            len(table) == original.n_obs and table["cell_id"].is_unique
            and np.array_equal(table["cell_id"].astype(str).to_numpy(), ids)
            and np.array_equal(table["cluster"].astype(str).to_numpy(), original.obs["cluster"].astype(str).to_numpy())
            and table["analysis_scope"].eq("within_dataset_model_inferred").all()
            and table["direction_source"].eq("stochastic_kinetic_model").all()
            and table["time_units"].eq("unitless_model_pseudotime").all()
        )
        for column, obs_column, lower, upper in (
            ("velocity_pseudotime", "velocity_pseudotime", 0.0, 1.0),
            ("root_score", "root_cells", 0.0, 1.0), ("end_score", "end_points", 0.0, 1.0),
            ("velocity_confidence", "velocity_confidence", 0.0, 1.0),
        ):
            values = pd.to_numeric(table[column], errors="coerce").to_numpy()
            if obs_column not in data.obs:
                summary_ok = False
                continue
            reference = pd.to_numeric(data.obs[obs_column], errors="coerce").to_numpy()
            summary_ok = summary_ok and bool(
                np.isfinite(values).all() and np.isfinite(reference).all()
                and (values >= lower - 1e-6).all() and (values <= upper + 1e-6).all()
                and np.allclose(values, reference, rtol=1e-7, atol=1e-8)
            )
        summary_ok = summary_ok and bool(table["velocity_pseudotime"].max() > table["velocity_pseudotime"].min())
        report.add_check("rna_velocity_per_cell_summary", bool(summary_ok), ValidationSeverity.ERROR,
                         "Reconstructed every summary row from persisted scVelo observations and checked scores/pseudotime ranges.")
        settings = data.uns.get("rna_velocity", {})
        assumptions = settings.get("model_assumptions", [])
        assumptions_ok = isinstance(assumptions, (list, tuple, np.ndarray)) and len(assumptions) >= 5
        if assumptions_ok:
            assumption_text = " ".join(str(item) for item in assumptions).lower()
            assumptions_ok = all(token in assumption_text for token in ("stochastic", "splicing", "neighborhood", "model mismatch", "not measured time", "causal"))
        forbidden_parameters = {"root", "root_key", "root_cell_id", "iroot", "end_key", "end_cell_id", "latent_time", "mode",
                                "filter_genes", "use_highly_variable", "spliced_layer", "unspliced_layer",
                                "fit_offset", "fit_offset2", "min_r2", "min_likelihood", "perc", "use_raw",
                                "groups", "tkey", "gene_subset", "vkey", "n_neighbors", "n_pcs", "basis", "approx"}
        scope_ok = bool(
            settings.get("method") == "scvelo_stochastic_v1" and assumptions_ok
            and settings.get("caller_root_used") is False and settings.get("root_source") == "kinetic_terminal_states"
            and settings.get("neighbors_source") == "native_project_clustering"
            and settings.get("kinetic_groupby") == "none_global_graph" and settings.get("summary_groupby") == "cluster"
            and settings.get("counts_layers") == {"spliced": "spliced_counts", "unspliced": "unspliced_counts", "expression": "counts"}
            and settings.get("mode") == "stochastic"
            and settings.get("feature_filtering") is False and settings.get("n_jobs") == 1
            and settings.get("confirmatory_allowed") is False and settings.get("analysis_scope") == "within_dataset_model_inferred"
            and not forbidden_parameters.intersection(contract.parameters)
            and contract.parameters.get("groupby", "cluster") == "cluster" and contract.parameters.get("n_jobs", 1) == 1
            and "latent_time" not in data.obs and "iroot" not in data.uns
            and result.inference_contract is None and result.scientific_result is None and contract.inference_contract is None
        )
        provenance_ok = bool(
            out_meta.created_by_task == table_meta.created_by_task == contract.task_id
            and out_meta.parent_uris == contract.input_artifacts and table_meta.parent_uris == [result.output_artifacts[0]]
            and out_meta.operation == "scvelo_stochastic_v1" and table_meta.operation == "summarize_stochastic_rna_velocity"
            and out_meta.software_versions.get("scvelo") not in (None, "unknown", "unavailable")
            and out_meta.software_versions == table_meta.software_versions
            and original.uns.get("clustering", {}).get("method") == "scanpy_leiden_umap_v1"
        )
        effective_parameters = {"groupby": "cluster", "n_jobs": 1, "n_dcs": 10,
                                "target_sum": 10000.0, "random_seed": 42, **contract.parameters}
        provenance_ok = provenance_ok and bool(
            out_meta.parameters == table_meta.parameters == effective_parameters
            and out_meta.random_seed == table_meta.random_seed == effective_parameters["random_seed"]
            and settings.get("target_sum") == effective_parameters["target_sum"]
            and settings.get("n_dcs") == min(effective_parameters["n_dcs"], original.n_obs - 2)
        )
        metrics = {
            "n_cells": original.n_obs, "n_genes": original.n_vars,
            "n_velocity_genes": int(data.var["velocity_genes"].sum()) if "velocity_genes" in data.var else 0,
            "mode": "stochastic", "direction_source": "stochastic_kinetic_model", "root_source": "kinetic_terminal_states",
            "analysis_scope": "within_dataset_model_inferred", "counts_preserved": True,
            "feature_filtering": False, "n_jobs": 1, "summary_groupby": "cluster",
            "kinetic_groupby": "none_global_graph", "confirmatory_allowed": False,
            "mean_velocity_confidence": float(pd.to_numeric(table["velocity_confidence"], errors="coerce").mean()),
            "pseudotime_min": float(pd.to_numeric(table["velocity_pseudotime"], errors="coerce").min()),
            "pseudotime_max": float(pd.to_numeric(table["velocity_pseudotime"], errors="coerce").max()),
        }
        for key, expected in metrics.items():
            for actual in (result.metrics.get(key), out_meta.summary_metrics.get(key), table_meta.summary_metrics.get(key), settings.get(key)):
                provenance_ok = provenance_ok and (bool(np.isclose(actual, expected, rtol=1e-7, atol=1e-8))
                                                    if isinstance(expected, float) and isinstance(actual, (int, float)) else actual == expected)
        report.add_check("rna_velocity_model_scope_and_provenance", bool(scope_ok and provenance_ok), ValidationSeverity.ERROR,
                         "Verified native stochastic model scope, kinetic root provenance, assumptions and library metadata; no independent biological validation is claimed.")
        report.add_check("rna_velocity_integrity", bool(identity_ok and counts_ok and layers_ok and moments_ok and residual_ok and graph_ok and summary_ok and scope_ok and provenance_ok),
                         ValidationSeverity.ERROR, "Persisted native RNA-velocity structural and summary checks passed; kinetic model fit and biological assumptions are not independently refitted here.")
