"""Trajectory-like exploratory analysis with honest method provenance.

The local implementation estimates a root-relative Euclidean distance in an
existing embedding and tests gene/rank associations.  It is intentionally
not labelled PAGA, DPT, CellRank, or any other graph/velocity algorithm.
Those historical IDs are accepted as request aliases by the capability
registry; returned results always carry ``root_distance_pseudotime_v1``.
"""

from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.numerics import benjamini_hochberg


def _to_dense(value) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    return np.asarray(value, dtype=np.float32)


def _bh(p_values: np.ndarray) -> np.ndarray:
    """Compatibility wrapper for the historical private trajectory helper."""

    return benjamini_hochberg(p_values)


EMPTY_DYNAMIC_COLUMNS = ["gene", "spearman_rho", "p_value", "fdr_q_value", "trend"]


class TrajectoryCapability(BaseCapability):
    """Estimate root-distance pseudotime and exploratory gene dynamics."""

    def __init__(self, implementation_id: str = "root_distance_pseudotime_v1"):
        requested_id = implementation_id
        super().__init__(
            capability_name="trajectory_inference",
            implementation_id="root_distance_pseudotime_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )
        self.requested_implementation_id = requested_id
        self.legacy_aliases = {
            "paga_dpt": self.implementation_id,
            "cellrank": self.implementation_id,
        }

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        _, payload = registry.get(in_uri)
        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)

        emb = data.obsm.get("X_pca")
        if emb is None:
            emb = data.X[:, : min(10, data.n_vars)]
        emb = _to_dense(emb)
        if emb.ndim != 2 or emb.shape[0] != data.n_obs or emb.shape[1] == 0:
            raise ValueError("Trajectory requires a non-empty cell-by-coordinate embedding")

        gene_names = (
            data.var["gene_name"].astype(str).values
            if "gene_name" in data.var.columns
            else np.array([f"G_{i}" for i in range(data.n_vars)], dtype=str)
        )

        params = contract.parameters
        seed = int(params.get("random_seed", 42))
        rng = np.random.default_rng(seed)
        root_idx = int(params.get("root_cell_index", 0))
        root_gene = params.get("root_gene")
        if root_gene is None and "Apoe" in set(gene_names):
            root_gene = "Apoe"
        if root_gene in set(gene_names):
            gene_idx = int(np.flatnonzero(gene_names == root_gene)[0])
            root_expr = _to_dense(data.X[:, gene_idx]).reshape(-1)
            root_idx = int(np.nanargmin(root_expr))
        root_idx = min(max(root_idx, 0), max(0, data.n_obs - 1))

        distances = cdist(emb[[root_idx]], emb, metric="euclidean")[0]
        span = float(np.max(distances) - np.min(distances)) if len(distances) else 0.0
        pseudotime = np.zeros(data.n_obs, dtype=np.float32) if span <= 1e-12 else ((distances - np.min(distances)) / span).astype(np.float32)

        # Stability is an estimate, not a prior.  Keep it missing when no
        # valid split has a defined rank correlation.
        correlations = []
        n_cells = len(pseudotime)
        n_repeats = max(0, int(params.get("stability_repeats", 5)))
        fraction = float(np.clip(params.get("stability_fraction", 0.8), 0.1, 1.0))
        if n_cells >= 3:
            sample_size = min(n_cells, max(2, int(round(fraction * n_cells))))
            for _ in range(n_repeats):
                sub_idx = np.sort(rng.choice(n_cells, size=sample_size, replace=False))
                # The root may be absent from a subsample.  In that case use
                # the nearest sampled cell to the original root coordinate.
                sub_root = int(np.argmin(cdist(emb[[root_idx]], emb[sub_idx])[0]))
                sub_dist = cdist(emb[[sub_idx[sub_root]]], emb[sub_idx])[0]
                sub_span = float(np.max(sub_dist) - np.min(sub_dist)) if len(sub_dist) else 0.0
                if sub_span <= 1e-12:
                    continue
                sub_pt = (sub_dist - np.min(sub_dist)) / sub_span
                corr, _ = stats.spearmanr(pseudotime[sub_idx], sub_pt)
                if np.isfinite(corr):
                    correlations.append(float(corr))
        stability_score: Optional[float] = float(np.mean(correlations)) if correlations else None

        # Test every gene, then apply BH to all tested genes.  Returning an
        # empty table with the same schema is a valid no-discovery result.
        all_rows = []
        X = data.X
        for j, gene in enumerate(gene_names):
            g_expr = _to_dense(X[:, j]).reshape(-1)
            if np.ptp(pseudotime) <= 1e-12 or np.ptp(g_expr) <= 1e-12:
                rho, p_value = 0.0, 1.0
            else:
                rho, p_value = stats.spearmanr(pseudotime, g_expr)
            rho = float(rho) if np.isfinite(rho) else 0.0
            p_value = float(p_value) if np.isfinite(p_value) else 1.0
            all_rows.append({"gene": str(gene), "spearman_rho": rho, "p_value": p_value})

        all_df = pd.DataFrame(all_rows, columns=["gene", "spearman_rho", "p_value"])
        if all_df.empty:
            dyn_df = pd.DataFrame(columns=EMPTY_DYNAMIC_COLUMNS)
        else:
            all_df["fdr_q_value"] = _bh(all_df["p_value"].values)
            min_abs_rho = float(params.get("min_abs_rho", 0.3))
            fdr_threshold = float(params.get("fdr_threshold", 0.05))
            discoveries = all_df[(all_df["fdr_q_value"] <= fdr_threshold) & (all_df["spearman_rho"].abs() >= min_abs_rho)].copy()
            discoveries["trend"] = np.where(discoveries["spearman_rho"] >= 0, "upregulated_along_pseudotime", "downregulated_along_pseudotime")
            dyn_df = discoveries[EMPTY_DYNAMIC_COLUMNS].sort_values("spearman_rho", ascending=False).reset_index(drop=True)

        res_data = data.copy()
        res_data.obs["pseudotime"] = pseudotime
        res_data.uns["trajectory"] = {
            "method": self.implementation_id,
            "requested_method": self.requested_implementation_id,
            "algorithm": "root_relative_euclidean_distance",
            "root_cell": root_idx,
            "stability_score": stability_score,
            "stability_repeats": len(correlations),
            "random_seed": seed,
        }

        uri_obj = ArtifactURI.parse(in_uri)
        # Target branches provide an explicit output namespace.  Keep the
        # historical URI as the fallback used by legacy single-target plans.
        out_uri = contract.expected_outputs[0] if contract.expected_outputs else f"table://{uri_obj.study_id}/trajectory_results/v1"
        registry.register(
            uri_str=out_uri,
            payload=dyn_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="estimate_root_distance_pseudotime",
            parent_uris=[in_uri],
            parameters={
                "method": self.implementation_id,
                "root_cell": root_idx,
                "random_seed": seed,
                "fdr_threshold": float(params.get("fdr_threshold", 0.05)),
                "min_abs_rho": float(params.get("min_abs_rho", 0.3)),
                "fdr_family": params.get("fdr_family", "study:trajectory"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
            summary_metrics={
                "stability_score": stability_score,
                "dynamic_genes_count": len(dyn_df),
                "top_dynamic_gene": dyn_df.iloc[0]["gene"] if not dyn_df.empty else None,
                "discoveries_supported_by_fdr": True,
            },
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["estimate_root_distance_pseudotime", "test_subsample_stability", "spearman_gene_association", "benjamini_hochberg_correction"],
            metrics={
                "stability_score": stability_score,
                "stability_repeats": len(correlations),
                "top_dynamic_genes": dyn_df["gene"].head(5).tolist() if not dyn_df.empty else [],
                "dynamic_genes_count": len(dyn_df),
                "fdr_family": params.get("fdr_family", "study:trajectory"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
        )
