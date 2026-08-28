"""
Differential Expression (DEG) and Differential Abundance Capabilities.
Implements pseudobulk analysis and pseudoreplication detection.
"""

from typing import Dict, Any, List, Optional
import numpy as np
import pandas as pd
from scipy import stats

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    """Calculates Benjamini-Hochberg FDR adjusted p-values."""
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.array([])
    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    fdr = np.zeros(n)
    curr_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        adj = (sorted_p[i] * n) / rank
        curr_min = min(curr_min, adj)
        fdr[i] = min(curr_min, 1.0)
    rev_indices = np.argsort(sorted_indices)
    return fdr[rev_indices]


class DifferentialAbundanceCapability(BaseCapability):
    """Tests for state enrichment (e.g. M3 enriched in AD) across biological replicates."""

    def __init__(self):
        super().__init__(
            capability_name="differential_abundance",
            implementation_id="state_abundance_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        obs = data.obs

        cond_col = "condition" if "condition" in obs.columns else obs.columns[0]
        mouse_col = "mouse_id" if "mouse_id" in obs.columns else "sample_id"
        state_col = contract.parameters.get("state_col", "microglia_state" if "microglia_state" in obs.columns else "leiden")

        # Compute state proportions per donor mouse
        ct_table = pd.crosstab([obs[mouse_col], obs[cond_col]], obs[state_col], normalize="index")
        ct_df = ct_table.reset_index()

        results = []
        conditions = obs[cond_col].unique()
        if len(conditions) >= 2:
            c1, c2 = conditions[0], conditions[1]
            for state in obs[state_col].unique():
                vals1 = ct_df[ct_df[cond_col] == c1][state].values
                vals2 = ct_df[ct_df[cond_col] == c2][state].values
                mean1 = float(np.mean(vals1)) if len(vals1) > 0 else 0.0
                mean2 = float(np.mean(vals2)) if len(vals2) > 0 else 0.0
                
                # t-test across mice
                if len(vals1) >= 2 and len(vals2) >= 2:
                    tt = stats.ttest_ind(vals1, vals2, equal_var=False)
                    p_val = float(tt.pvalue) if not np.isnan(tt.pvalue) else 1.0
                else:
                    p_val = 1.0

                log2_ratio = float(np.log2((mean1 + 1e-4) / (mean2 + 1e-4)))
                results.append({
                    "state": state,
                    f"mean_prop_{c1}": mean1,
                    f"mean_prop_{c2}": mean2,
                    "log2_ratio": log2_ratio,
                    "p_value": p_val,
                    "enriched_in": c1 if log2_ratio > 0 else c2,
                })

        res_df = pd.DataFrame(results)
        if not res_df.empty:
            res_df["fdr"] = benjamini_hochberg(res_df["p_value"].values)

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"table://{uri_obj.study_id}/differential_abundance/v1"

        registry.register(
            uri_str=out_uri,
            payload=res_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="test_differential_abundance",
            parent_uris=[in_uri],
            summary_metrics={
                "tested_states": len(res_df),
                "top_enriched_state": res_df.sort_values("p_value").iloc[0]["state"] if not res_df.empty else "none",
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["aggregate_donor_proportions", "welch_t_test", "benjamini_hochberg"],
            metrics={
                "abundance_results": res_df.to_dict(orient="records"),
            }
        )


class DifferentialExpressionCapability(BaseCapability):
    """Performs Differential Expression (preferring pseudobulk when donor replicates >= 3)."""

    def __init__(self):
        super().__init__(
            capability_name="deg",
            implementation_id="deg_pseudobulk_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        meta, payload = registry.get(in_uri)

        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        obs = data.obs
        gene_names = data.var["gene_name"].values if "gene_name" in data.var.columns else np.array([f"G_{i}" for i in range(data.n_vars)])

        cond_col = "condition" if "condition" in obs.columns else obs.columns[0]
        mouse_col = "mouse_id" if "mouse_id" in obs.columns else None
        
        conditions = obs[cond_col].unique()
        if len(conditions) < 2:
            raise ValueError(f"Need at least 2 conditions for DEG, found: {conditions}")

        cond_ad = "AD" if "AD" in conditions else conditions[0]
        cond_ctrl = "control" if "control" in conditions else conditions[1]

        # Audit biological replicates
        if mouse_col:
            n_reps_ad = obs[obs[cond_col] == cond_ad][mouse_col].nunique()
            n_reps_ctrl = obs[obs[cond_col] == cond_ctrl][mouse_col].nunique()
        else:
            n_reps_ad, n_reps_ctrl = 1, 1

        is_pseudobulk = (n_reps_ad >= 3 and n_reps_ctrl >= 3)
        
        results = []
        if is_pseudobulk and mouse_col:
            # 1. Pseudobulk Aggregation per donor
            donors_ad = obs[obs[cond_col] == cond_ad][mouse_col].unique()
            donors_ctrl = obs[obs[cond_col] == cond_ctrl][mouse_col].unique()

            bulk_ad = np.array([data.X[obs[mouse_col] == d].mean(axis=0) for d in donors_ad])
            bulk_ctrl = np.array([data.X[obs[mouse_col] == d].mean(axis=0) for d in donors_ctrl])

            mean_ad = np.mean(bulk_ad, axis=0)
            mean_ctrl = np.mean(bulk_ctrl, axis=0)
            log2_fc = np.log2((mean_ad + 1e-4) / (mean_ctrl + 1e-4))

            # Two-sample t-test across donor means
            t_stat, p_vals = stats.ttest_ind(bulk_ad, bulk_ctrl, axis=0, equal_var=False)
            p_vals = np.nan_to_num(p_vals, nan=1.0)
            fdr = benjamini_hochberg(p_vals)

            for i in range(len(gene_names)):
                results.append({
                    "gene": str(gene_names[i]),
                    "mean_AD": float(mean_ad[i]),
                    "mean_Ctrl": float(mean_ctrl[i]),
                    "log2_fold_change": float(log2_fc[i]),
                    "p_value": float(p_vals[i]),
                    "fdr_q_value": float(fdr[i]),
                    "statistical_unit": "donor_pseudobulk",
                    "n_donors_AD": int(n_reps_ad),
                    "n_donors_Ctrl": int(n_reps_ctrl),
                })
        else:
            # 2. Cell-level Wilcoxon / t-test (with pseudoreplication warning)
            mask_ad = (obs[cond_col] == cond_ad).values
            mask_ctrl = (obs[cond_col] == cond_ctrl).values

            cells_ad = data.X[mask_ad]
            cells_ctrl = data.X[mask_ctrl]

            mean_ad = np.mean(cells_ad, axis=0)
            mean_ctrl = np.mean(cells_ctrl, axis=0)
            log2_fc = np.log2((mean_ad + 1e-4) / (mean_ctrl + 1e-4))

            # Wilcoxon / Mann-Whitney U test per gene
            p_vals = []
            for j in range(data.n_vars):
                try:
                    u_stat, p = stats.mannwhitneyu(cells_ad[:, j], cells_ctrl[:, j], alternative="two-sided")
                    p_vals.append(p)
                except Exception:
                    p_vals.append(1.0)
            p_vals = np.array(p_vals)
            fdr = benjamini_hochberg(p_vals)

            for i in range(len(gene_names)):
                results.append({
                    "gene": str(gene_names[i]),
                    "mean_AD": float(mean_ad[i]),
                    "mean_Ctrl": float(mean_ctrl[i]),
                    "log2_fold_change": float(log2_fc[i]),
                    "p_value": float(p_vals[i]),
                    "fdr_q_value": float(fdr[i]),
                    "statistical_unit": "single_cell_exploratory",
                    "pseudoreplication_warning": True,
                })

        deg_df = pd.DataFrame(results).sort_values("fdr_q_value").reset_index(drop=True)
        sig_degs = deg_df[(deg_df["fdr_q_value"] < 0.05) & (deg_df["log2_fold_change"].abs() > 0.5)]

        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = f"table://{uri_obj.study_id}/deg_results/v1"

        registry.register(
            uri_str=out_uri,
            payload=deg_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="differential_expression",
            parent_uris=[in_uri],
            parameters={"is_pseudobulk": is_pseudobulk, "cond_ad": cond_ad, "cond_ctrl": cond_ctrl},
            summary_metrics={
                "total_genes": len(deg_df),
                "significant_degs_fdr05": len(sig_degs),
                "statistical_unit": "donor_pseudobulk" if is_pseudobulk else "single_cell",
                "top_upregulated": sig_degs[sig_degs["log2_fold_change"] > 0]["gene"].head(5).tolist(),
            }
        )

        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=self.implementation_id,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=["pseudobulk_aggregation" if is_pseudobulk else "cell_level_mannwhitney", "fdr_correction"],
            metrics={
                "significant_degs": len(sig_degs),
                "is_pseudobulk": is_pseudobulk,
                "top_genes": sig_degs["gene"].head(5).tolist(),
            }
        )
