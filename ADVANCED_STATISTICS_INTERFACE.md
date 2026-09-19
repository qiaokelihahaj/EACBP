# Advanced statistics interface

The optional `advanced-statistics` extra provides the three capabilities in
`eacbp.capabilities.advanced_statistics`.  The capabilities are intentionally
not installed into the default registry by this module; a profile or factory
can register them explicitly.

## `deg` / `pydeseq2_pseudobulk_v1`

`PyDESeq2PseudobulkCapability` consumes one AnnData/SCData artifact.  It
requires an integer, non-negative raw-count layer (`counts_layer`, default
`counts`), `condition_col` (default `condition`) and `donor_col` (or one of
`donor`, `donor_id`, `mouse_id`, `sample_id`, `sample`).  `condition_a` and
`condition_b` are required when more than two levels are observed.  If
`paired=True`, every donor must have both conditions and the formula includes a
donor blocking term.  Independent designs reject donor IDs that occur in both
conditions.

Useful parameters are:

* `design_formula` / `design`: a PyDESeq2 formula such as
  `~condition + donor + batch`.  When omitted, the formula is built from
  `paired` and `covariates`.
* `covariates`: donor-condition-constant metadata columns to append to the
  default formula.
* `contrast`: `['condition', 'tested_level', 'reference_level']` or a numeric
  contrast vector.  The default reports `condition_a` versus `condition_b`.
* `min_donors` (default 2), `min_replicates` (default 2), `alpha` (default
  0.05), `n_cpus` (default 1), `fit_type`, `size_factors_fit_type`,
  `cooks_filter`, and `independent_filter`.

The first `expected_outputs` URI is used verbatim; without one the default is
`table://<study_id>/pydeseq2_deg/v1`.  The output is a complete gene table,
including non-significant and PyDESeq2-not-estimated rows.  It contains
`gene`, `log2_fold_change`, `lfc_se`, `ci_low`, `ci_high`, `statistic`,
`p_value`, `fdr_q_value`/`padj`, design metadata and sample counts.  A missing
PyDESeq2 p-value remains `NaN` with no fabricated replacement.

`contract_operations` exposes the normal and underpowered branches, including
`skip_insufficient_replicates`.  Executed operations are
`validate_raw_integer_counts`, `validate_metadata`,
`validate_replicate_qualification`, `aggregate_counts_by_donor_condition`,
`validate_design_rank`,
`fit_pydeseq2_negative_binomial`, `wald_contrast`,
`benjamini_hochberg_fdr`, and `retain_all_genes`.  Metrics report the design,
rank, conditions, donor counts, counts source and numbers of estimated,
non-estimated/FDR-missing and significant genes.  Required validation requirements are raw integer counts,
metadata completeness, donor replication, contrast validity, and full-rank
design; callers may add `multiple_testing_correction` and
`pseudoreplication_audit` to the task contract.

## `functional_activity` / `decoupler_ulm_v2`

`DecouplerFunctionalAnalysisCapability` consumes the raw-count artifact and,
optionally, a network or DEG table artifact.  A network must be supplied as a
local `network_path`, a `network` DataFrame/mapping, or an additional table
artifact with `source`, `target`, and optional `weight` columns.  The contract
must explicitly provide `species`, `network_source`, and `network_version`.
Remote URLs and decoupler resource downloads are rejected.  `network_sha256`
can pin a local resource hash.  `tmin` defaults to 1 for small local networks.

The capability applies decoupler's current `dc.mt.ulm` API to donor-condition
log1p-CPM expression, then compares donor activity scores with an OLS model
materialized from the same formula as the DEG contract.  Paired donor fixed
effects and requested covariates therefore remain in the activity contrast.
OLS tests and confidence intervals use the same residual degrees of freedom
and Student t distribution; the independent auditor verifies those intervals
and rejects missing or invalid degrees of freedom. The OLS p-values are
corrected with BH FDR. The first expected URI (default
`table://<study_id>/functional_activity/v1`) contains one row per retained
network source and preserves `p_value`, `fdr_q_value`, donor means, test name,
species, network source/version/hash, and feature type (`analysis_kind` or
`network_kind`).  A second expected URI receives the per-donor activity and
decoupler ULM p-value table.  No causal interpretation is emitted.

`contract_operations` includes both normal execution and
`skip_insufficient_donors`.  Executed operations are
`validate_local_network_provenance`, `aggregate_counts_by_donor_condition`,
`log1p_cpm_transform_for_decoupler`, `run_decoupler_ulm`,
`build_formula_activity_design`, `compare_donor_activities`,
`benjamini_hochberg_fdr`, and `retain_all_network_sources`.

## `donor_sensitivity` / `pydeseq2_leave_one_donor_out_v1`

`PyDESeq2LeaveOneDonorOutCapability` consumes the same raw-count artifact and
parameters as the DEG capability.  It requires at least three donors per
condition (`min_loo_donors` can raise this threshold).  It fits the full
model, then re-fits after each donor is removed.  In paired designs removing
one donor removes both condition pseudobulk rows.  If the threshold is not
met, it registers a one-row `status=skipped` table with `skip_reason` and
returns `scientific_robustness_claim_supported=false`.

The first expected URI (default
`table://<study_id>/donor_sensitivity/v1`) contains the full-model row and
all successful leave-one-donor-out gene tables, retaining the complete gene
set and NaN values for unestimated tests.  Per-gene summaries include
`estimated_coverage`, `direction_consistency`, min/max/median LFC and
`significance_retention`; a second expected URI receives the one-row-per-gene
summary.  Metrics include requested, successful and failed fits plus
`complete_fit_coverage`.  These are sensitivity descriptors and never a
scientific robustness claim.  Executed operations are
`validate_leave_one_donor_qualification`, `fit_full_pydeseq2_model`,
`remove_complete_donor_for_each_fit`, `refit_pydeseq2_per_donor`,
`retain_all_genes_per_fit`, and
`record_fit_failures_without_claiming_robustness`.
The exposed `contract_operations` list also includes the explicit
`skip_insufficient_donors` branch.
