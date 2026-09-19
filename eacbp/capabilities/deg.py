"""Differential abundance/expression with explicit statistical units.

The primary DEG path aggregates raw counts by donor and compares donor-level
normalised pseudobulk values with Welch's t-test.  If donor metadata or
replicate support is insufficient, the capability emits a clearly labelled
``single_cell_exploratory`` result.  It never invents condition or donor
labels and never calls a Welch test Mann-Whitney.
"""

from typing import Any, Optional

import numpy as np
import pandas as pd
from scipy import stats

from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus
from eacbp.schemas.artifact import ArtifactType
from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.sc_data import SCData
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.uri import ArtifactURI
from eacbp.numerics import benjamini_hochberg


def _dense(value) -> np.ndarray:
    if hasattr(value, "toarray"):
        value = value.toarray()
    return np.asarray(value, dtype=np.float64)


def _counts_layer(data: SCData):
    layers = getattr(data, "layers", {}) or {}
    return layers.get("counts")


def _resolve_conditions(obs: pd.DataFrame, params: dict) -> tuple[str, str]:
    cond_col = params.get("condition_col", "condition")
    if cond_col not in obs.columns:
        raise ValueError(
            f"No condition metadata column '{cond_col}' is available; provide condition_col explicitly."
        )
    conditions = list(pd.unique(obs[cond_col].dropna()))
    requested_a = params.get("condition_a")
    requested_b = params.get("condition_b")
    if requested_a is not None or requested_b is not None:
        if requested_a is None or requested_b is None:
            raise ValueError("condition_a and condition_b must be supplied together")
        if requested_a not in conditions or requested_b not in conditions:
            raise ValueError(f"Requested conditions {requested_a!r}, {requested_b!r} are not present: {conditions}")
        if str(requested_a) == str(requested_b):
            raise ValueError("condition_a and condition_b must be different observed conditions")
        return str(requested_a), str(requested_b)
    if len(conditions) != 2:
        raise ValueError(
            f"DEG requires exactly two observed conditions when condition_a/condition_b are omitted; found {conditions}"
        )
    return str(conditions[0]), str(conditions[1])


def _resolve_donor_column(obs: pd.DataFrame, params: dict) -> Optional[str]:
    explicit = params.get("donor_col")
    if explicit is not None:
        if explicit not in obs.columns:
            raise ValueError(f"Donor metadata column '{explicit}' is not present")
        return str(explicit)
    for name in ("donor", "donor_id", "mouse_id", "sample_id", "sample"):
        if name in obs.columns:
            return name
    return None


def _donor_pseudobulk(counts, obs, condition_col: str, donor_col: str, condition: str) -> tuple[np.ndarray, list[str]]:
    donor_ids = [str(v) for v in pd.unique(obs.loc[obs[condition_col].astype(str) == str(condition), donor_col].dropna())]
    rows = []
    for donor in donor_ids:
        mask = (obs[condition_col].astype(str).values == str(condition)) & (obs[donor_col].astype(str).values == donor)
        selected = counts[mask]
        if hasattr(selected, "tocsr"):
            row = np.asarray(selected.sum(axis=0)).ravel()
        else:
            row = np.asarray(selected, dtype=float).sum(axis=0)
        rows.append(row)
    if not rows:
        return np.empty((0, counts.shape[1]), dtype=float), donor_ids
    return np.vstack(rows), donor_ids


def _library_size_normalize(bulk: np.ndarray, target_sum: float = 1_000_000.0) -> np.ndarray:
    if bulk.size == 0:
        return bulk.astype(float)
    totals = bulk.sum(axis=1)
    safe = np.where(totals > 0, totals, 1.0)
    return bulk / safe[:, None] * float(target_sum)


class DifferentialAbundanceCapability(BaseCapability):
    """Compare cell-state proportions across donor-level biological units."""

    def __init__(self, implementation_id: str = "state_abundance_donor_welch_v1"):
        super().__init__(
            capability_name="differential_abundance",
            implementation_id="state_abundance_donor_welch_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )
        self.requested_implementation_id = implementation_id

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        _, payload = registry.get(in_uri)
        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        obs = data.obs
        params = contract.parameters
        condition_col = str(params.get("condition_col", "condition"))
        if condition_col not in obs.columns:
            raise ValueError(f"Condition metadata column '{condition_col}' is required")
        donor_col = _resolve_donor_column(obs, params)
        state_col = str(params.get("state_col", "cell_type" if "cell_type" in obs.columns else "microglia_state" if "microglia_state" in obs.columns else "cluster"))
        if donor_col is None:
            raise ValueError("Differential abundance requires donor/sample metadata; no donor column was found")
        if state_col not in obs.columns:
            raise ValueError(f"State metadata column '{state_col}' is required")
        cond_a, cond_b = _resolve_conditions(obs, params)

        donors_a = set(obs.loc[obs[condition_col].astype(str) == cond_a, donor_col].dropna().astype(str))
        donors_b = set(obs.loc[obs[condition_col].astype(str) == cond_b, donor_col].dropna().astype(str))
        if donors_a & donors_b:
            raise ValueError("Independent abundance Welch test cannot analyze paired donors")
        ct_table = pd.crosstab([obs[donor_col], obs[condition_col]], obs[state_col], normalize="index").reset_index()
        states = list(pd.unique(obs[state_col]))
        rows = []
        for state in states:
            vals_a = ct_table.loc[ct_table[condition_col].astype(str) == cond_a, state].to_numpy(dtype=float)
            vals_b = ct_table.loc[ct_table[condition_col].astype(str) == cond_b, state].to_numpy(dtype=float)
            p_val = 1.0
            if len(vals_a) >= 2 and len(vals_b) >= 2:
                p_val = float(stats.ttest_ind(vals_a, vals_b, equal_var=False).pvalue)
                if not np.isfinite(p_val):
                    p_val = 1.0
            mean_a = float(np.mean(vals_a)) if len(vals_a) else np.nan
            mean_b = float(np.mean(vals_b)) if len(vals_b) else np.nan
            ratio = float(np.log2((mean_a + 1e-8) / (mean_b + 1e-8))) if np.isfinite(mean_a) and np.isfinite(mean_b) else np.nan
            rows.append({
                "state": state,
                "condition_a": cond_a,
                "condition_b": cond_b,
                "mean_prop_condition_a": mean_a,
                "mean_prop_condition_b": mean_b,
                "log2_ratio_condition_a_vs_b": ratio,
                "p_value": p_val,
                "donor_col": donor_col,
                "statistical_unit": "donor_proportion",
                "enriched_in": cond_a if np.isfinite(ratio) and ratio > 0 else cond_b if np.isfinite(ratio) else None,
            })
        res_df = pd.DataFrame(rows)
        if not res_df.empty:
            res_df["fdr_q_value"] = benjamini_hochberg(res_df["p_value"].values)
        uri_obj = ArtifactURI.parse(in_uri)
        out_uri = (contract.expected_outputs[0] if contract.expected_outputs
                   else f"table://{uri_obj.study_id}/differential_abundance/v1")
        registry.register(
            uri_str=out_uri,
            payload=res_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="test_donor_state_proportions_welch",
            parent_uris=[in_uri],
            summary_metrics={
                "tested_states": len(res_df),
                "statistical_unit": "donor_proportion",
                "condition_a": cond_a,
                "condition_b": cond_b,
                "fdr_family": params.get("fdr_family", "study:differential_abundance"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
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
                "condition_a": cond_a,
                "condition_b": cond_b,
                "donor_col": donor_col,
                "fdr_family": params.get("fdr_family", "study:differential_abundance"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
        )


class DifferentialExpressionCapability(BaseCapability):
    """Donor pseudobulk Welch DEG with an explicitly exploratory fallback."""

    def __init__(self, implementation_id: str = "donor_pseudobulk_welch_v1"):
        super().__init__(
            capability_name="deg",
            implementation_id="donor_pseudobulk_welch_v1",
            implementation_type=ImplementationType.PYTHON_TOOL,
            accepts_types=[ArtifactType.ANNDATA],
            output_types=[ArtifactType.TABLE],
        )
        self.requested_implementation_id = implementation_id
        self.legacy_aliases = {"deg_pseudobulk_v1": self.implementation_id}

    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        in_uri = contract.input_artifacts[0]
        _, payload = registry.get(in_uri)
        data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
        obs = data.obs
        params = contract.parameters
        condition_col = str(params.get("condition_col", "condition"))
        cond_a, cond_b = _resolve_conditions(obs, params)
        donor_col = _resolve_donor_column(obs, params)
        gene_names = (
            data.var["gene_name"].astype(str).values
            if "gene_name" in data.var.columns
            else np.array([f"G_{i}" for i in range(data.n_vars)], dtype=str)
        )
        counts = _counts_layer(data)
        if counts is None:
            counts = data.X if params.get("allow_x_as_counts", False) else None

        n_donor_a = int(obs.loc[obs[condition_col].astype(str) == cond_a, donor_col].nunique()) if donor_col else 0
        n_donor_b = int(obs.loc[obs[condition_col].astype(str) == cond_b, donor_col].nunique()) if donor_col else 0
        if donor_col:
            donors_a_set = set(obs.loc[obs[condition_col].astype(str) == cond_a, donor_col].dropna().astype(str))
            donors_b_set = set(obs.loc[obs[condition_col].astype(str) == cond_b, donor_col].dropna().astype(str))
            overlap = donors_a_set.intersection(donors_b_set)
            if overlap:
                raise ValueError(
                    "The same donor IDs occur in both conditions; independent Welch pseudobulk is invalid. "
                    "Use a paired model; the allow_paired_donors flag cannot bypass this restriction."
                )
        min_replicates = max(2, int(params.get("min_replicates", 2)))
        pseudobulk_supported = donor_col is not None and counts is not None and n_donor_a >= min_replicates and n_donor_b >= min_replicates
        require_pseudobulk = bool(params.get("require_pseudobulk", False))
        if require_pseudobulk and not pseudobulk_supported:
            missing = "counts layer" if counts is None else "at least two donors per condition"
            raise ValueError(f"Donor pseudobulk requested but {missing} is unavailable")

        if pseudobulk_supported:
            bulk_a, donors_a = _donor_pseudobulk(counts, obs, condition_col, donor_col, cond_a)
            bulk_b, donors_b = _donor_pseudobulk(counts, obs, condition_col, donor_col, cond_b)
            expr_a = _library_size_normalize(bulk_a)
            expr_b = _library_size_normalize(bulk_b)
            unit = "donor_pseudobulk"
            method_used = "donor_pseudobulk_welch_v1"
            rows = []
            pseudocount = float(params.get("pseudocount", 0.5))
            mean_a = np.mean(expr_a, axis=0)
            mean_b = np.mean(expr_b, axis=0)
            log2_fc = np.log2((mean_a + pseudocount) / (mean_b + pseudocount))
            for j, gene in enumerate(gene_names):
                p_val = float(stats.ttest_ind(expr_a[:, j], expr_b[:, j], equal_var=False).pvalue)
                rows.append({
                    "gene": str(gene),
                    "condition_a": cond_a,
                    "condition_b": cond_b,
                    "mean_condition_a": float(mean_a[j]),
                    "mean_condition_b": float(mean_b[j]),
                    "log2_fold_change": float(log2_fc[j]),
                    "effect_definition": "log2_ratio_of_mean_donor_CPM_plus_pseudocount",
                    "p_value": p_val if np.isfinite(p_val) else 1.0,
                    "statistical_unit": unit,
                    "donor_col": donor_col,
                    "n_donors_condition_a": len(donors_a),
                    "n_donors_condition_b": len(donors_b),
                    "counts_source": "layers.counts",
                })
            operations = ["aggregate_counts_by_donor", "library_size_normalize_donor_pseudobulk", "welch_t_test_donor_pseudobulk", "benjamini_hochberg"]
            is_pseudobulk = True
        else:
            # This path is intentionally exploratory: cells are independent
            # observations for the test only, so no donor-level claim is made.
            if counts is None:
                counts = data.X
                counts_source = "X_unverified_scale"
            else:
                counts_source = "layers.counts"
            mask_a = obs[condition_col].astype(str).values == cond_a
            mask_b = obs[condition_col].astype(str).values == cond_b
            values_a = _dense(counts[mask_a])
            values_b = _dense(counts[mask_b])
            mean_a = values_a.mean(axis=0) if len(values_a) else np.full(data.n_vars, np.nan)
            mean_b = values_b.mean(axis=0) if len(values_b) else np.full(data.n_vars, np.nan)
            log2_fc = np.log2((mean_a + 0.5) / (mean_b + 0.5))
            rows = []
            for j, gene in enumerate(gene_names):
                p_val = float(stats.ttest_ind(values_a[:, j], values_b[:, j], equal_var=False).pvalue) if len(values_a) and len(values_b) else 1.0
                rows.append({
                    "gene": str(gene),
                    "condition_a": cond_a,
                    "condition_b": cond_b,
                    "mean_condition_a": float(mean_a[j]),
                    "mean_condition_b": float(mean_b[j]),
                    "log2_fold_change": float(log2_fc[j]),
                    "effect_definition": "log2_ratio_of_mean_cell_values_plus_pseudocount_exploratory",
                    "p_value": p_val if np.isfinite(p_val) else 1.0,
                    "statistical_unit": "single_cell_exploratory",
                    "donor_col": donor_col,
                    "n_donors_condition_a": n_donor_a,
                    "n_donors_condition_b": n_donor_b,
                    "counts_source": counts_source,
                    "pseudoreplication_warning": True,
                })
            operations = ["compute_cell_level_effect", "welch_t_test_single_cell_exploratory", "benjamini_hochberg"]
            method_used = "single_cell_exploratory_welch_v1"
            unit = "single_cell_exploratory"
            is_pseudobulk = False

        deg_df = pd.DataFrame(rows)
        if not deg_df.empty:
            deg_df["fdr_q_value"] = benjamini_hochberg(deg_df["p_value"].values)
            deg_df = deg_df.sort_values("fdr_q_value").reset_index(drop=True)
        sig_degs = deg_df[(deg_df["fdr_q_value"] < 0.05) & (deg_df["log2_fold_change"].abs() > 0.5)] if not deg_df.empty else deg_df

        uri_obj = ArtifactURI.parse(in_uri)
        # The planner supplies a branch-local URI for multi-target studies.
        # Falling back to the historical URI preserves single-target and
        # broad-study behavior.
        out_uri = contract.expected_outputs[0] if contract.expected_outputs else f"table://{uri_obj.study_id}/deg_results/v1"
        registry.register(
            uri_str=out_uri,
            payload=deg_df,
            artifact_type=ArtifactType.TABLE,
            study_id=uri_obj.study_id,
            created_by_task=contract.task_id,
            operation="differential_expression_welch",
            parent_uris=[in_uri],
            parameters={
                "condition_col": condition_col,
                "condition_a": cond_a,
                "condition_b": cond_b,
                "donor_col": donor_col,
                "statistical_unit": unit,
                "counts_required_for_pseudobulk": True,
                "fdr_family": params.get("fdr_family", "study:deg"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
            summary_metrics={
                "total_genes": len(deg_df),
                "significant_degs_fdr05": len(sig_degs),
                "statistical_unit": unit,
                "top_upregulated": sig_degs[sig_degs["log2_fold_change"] > 0]["gene"].head(5).tolist() if not sig_degs.empty else [],
            },
        )
        return TaskResult(
            task_id=contract.task_id,
            status=TaskStatus.SUCCESS,
            capability=self.capability_name,
            method_used=method_used,
            input_artifacts=[in_uri],
            output_artifacts=[out_uri],
            executed_operations=operations,
            metrics={
                "significant_degs": len(sig_degs),
                "is_pseudobulk": is_pseudobulk,
                "statistical_unit": unit,
                "condition_a": cond_a,
                "condition_b": cond_b,
                "top_genes": sig_degs["gene"].head(5).tolist() if not sig_degs.empty else [],
                "fdr_family": params.get("fdr_family", "study:deg"),
                "target_cell_type": params.get("target_cell_type"),
                "target_branch": params.get("target_branch"),
            },
        )
