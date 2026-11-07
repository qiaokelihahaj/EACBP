"""Real-library smoke tests for the advanced donor-level statistics.

These tests are intentionally skipped when the opt-in dependencies are not
installed.  When enabled, they execute PyDESeq2 and decoupler on a small real
count matrix; no algorithm is mocked.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

pytest.importorskip("pydeseq2")
pytest.importorskip("decoupler")

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.advanced_statistics import (
    AdvancedStatisticsInputError,
    DecouplerFunctionalAnalysisCapability,
    PyDESeq2LeaveOneDonorOutCapability,
    PyDESeq2PseudobulkCapability,
)
from eacbp.auditor.advanced_statistics import _finite_ci_errors
from eacbp.auditor.advanced_statistics import AdvancedStatisticsValidator
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import TaskContract, TaskStatus


def _data(n_pairs=3):
    rng = np.random.default_rng(19)
    rows = []
    metadata = []
    for donor_idx in range(n_pairs):
        donor = f"d{donor_idx + 1}"
        for condition in ("A", "B"):
            for cell_idx in range(4):
                # Gene g0 is higher in A and g1 is higher in B.  The other
                # genes are deliberately unchanged, so the output must retain
                # nonsignificant rows as well.
                center = np.array([80 if condition == "A" else 8, 8 if condition == "A" else 80, 25, 25, 25, 25])
                rows.append(rng.poisson(center).astype(np.int64))
                batch = f"b{((donor_idx + (condition == 'B')) % 2) + 1}"
                metadata.append({"cell_id": f"{donor}_{condition}_{cell_idx}", "condition": condition, "donor_id": donor, "batch": batch})
    x = np.asarray(rows, dtype=np.float32) / 10.0
    counts = np.asarray(rows, dtype=np.int64)
    return SCData(
        X=x,
        obs=pd.DataFrame(metadata),
        var=pd.DataFrame({"gene_name": [f"g{i}" for i in range(counts.shape[1])]}),
        layers={"counts": counts},
    )


def _registry(tmp_path, data):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://advanced/input/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "advanced", "input", "input")
    return registry, uri


def test_pydeseq2_pseudobulk_real_fit_preserves_all_genes_and_expected_output(tmp_path):
    registry, uri = _registry(tmp_path, _data())
    out_uri = "table://advanced/pydeseq2_expected/v1"
    contract = TaskContract(
        task_id="advanced_deg",
        capability="deg",
        method="pydeseq2_pseudobulk_v1",
        input_artifacts=[uri],
        expected_outputs=[out_uri],
        parameters={
            "condition_a": "A",
            "condition_b": "B",
            "donor_col": "donor_id",
            "paired": True,
            "covariates": ["batch"],
            "n_cpus": 1,
        },
    )
    result = PyDESeq2PseudobulkCapability().execute(contract, registry)
    assert result.status == TaskStatus.SUCCESS
    assert result.output_artifacts == [out_uri]
    table = registry.load_payload(out_uri)
    assert len(table) == 6
    assert {"log2_fold_change", "lfc_se", "ci_low", "ci_high", "statistic", "p_value", "fdr_q_value"}.issubset(table.columns)
    assert table["gene"].tolist() == [f"g{i}" for i in range(6)]
    assert table["statistical_unit"].eq("donor_pseudobulk").all()
    assert table["pseudobulk_unit"].eq("donor_condition").all()
    assert result.metrics["all_genes_retained"] is True
    assert "donor" in result.metrics["design"]
    from eacbp.auditor import ScientificAuditor
    contract.validation_requirements = ["advanced_statistics_integrity", "multiple_testing_correction", "pseudoreplication_audit"]
    audited = ScientificAuditor().audit_task(contract, result, registry)
    assert audited.overall_passed, [(c.check_name, c.message) for c in audited.checks if not c.passed]


def test_custom_source_condition_column_uses_canonical_model_contrast(tmp_path):
    from eacbp.auditor import ScientificAuditor

    data = _data()
    data.obs["label"] = data.obs["condition"]
    registry, uri = _registry(tmp_path, data)
    out_uri = "table://advanced/custom_condition_alias/v1"
    contract = TaskContract(
        task_id="custom_condition_alias_deg",
        capability="deg",
        method="pydeseq2_pseudobulk_v1",
        input_artifacts=[uri],
        expected_outputs=[out_uri],
        parameters={
            "condition_col": "label",
            "condition_a": "A",
            "condition_b": "B",
            "contrast": ["label", "B", "A"],
            "donor_col": "donor_id",
            "paired": True,
            "counts_layer": "counts",
            "n_cpus": 1,
        },
    )

    result = PyDESeq2PseudobulkCapability().execute(contract, registry)

    assert result.status == TaskStatus.SUCCESS
    assert result.metrics["contrast_spec"]["factor"] == "condition"
    assert result.metrics["contrast_spec"]["vector"] == ["condition", "B", "A"]
    table = registry.load_payload(out_uri)
    assert table["contrast_factor"].eq("condition").all()
    assert table.loc[table["gene"] == "g0", "log2_fold_change"].iloc[0] < 0
    assert table.loc[table["gene"] == "g1", "log2_fold_change"].iloc[0] > 0

    contract.validation_requirements = ["advanced_statistics_integrity"]
    auditor = ScientificAuditor()
    audited = auditor.audit_task(contract, result, registry)
    assert audited.overall_passed, [(check.check_name, check.message) for check in audited.checks if not check.passed]

    tampered = result.model_copy(deep=True)
    tampered.metrics["contrast_spec"]["factor"] = "label"
    tampered.metrics["contrast_spec"]["vector"][0] = "label"
    rejected = auditor.audit_task(contract, tampered, registry)
    assert not rejected.overall_passed
    assert any(
        not check.passed and "reported contrast_spec differs from contract" in check.message
        for check in rejected.checks
    )


def test_pydeseq2_rejects_non_integer_counts(tmp_path):
    data = _data()
    data.layers["counts"] = data.layers["counts"].astype(float)
    data.layers["counts"][0, 0] += 0.25
    registry, uri = _registry(tmp_path, data)
    contract = TaskContract(
        task_id="advanced_bad_counts",
        capability="deg",
        input_artifacts=[uri],
        parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True},
    )
    with pytest.raises(AdvancedStatisticsInputError, match="integer-valued"):
        PyDESeq2PseudobulkCapability().execute(contract, registry)


def test_decoupler_ulm_uses_local_network_and_compares_donors(tmp_path):
    registry, uri = _registry(tmp_path, _data())
    network = pd.DataFrame(
        {
            "source": ["pathway_A", "pathway_A", "pathway_A", "tf_B", "tf_B"],
            "target": ["g0", "g2", "g3", "g1", "g4"],
            "weight": [1.0, 1.0, 1.0, 1.0, 1.0],
        }
    )
    out_uri = "table://advanced/activity_expected/v1"
    contract = TaskContract(
        task_id="advanced_activity",
        capability="functional_activity",
        method="decoupler_ulm_v2",
        input_artifacts=[uri],
        expected_outputs=[out_uri],
        parameters={
            "condition_a": "A",
            "condition_b": "B",
            "donor_col": "donor_id",
            "paired": True,
            "covariates": ["batch"],
            "network": network,
            "species": "Mus musculus",
            "network_source": "local_smoke_network",
            "network_version": "1.0",
            "tmin": 1,
        },
    )
    result = DecouplerFunctionalAnalysisCapability().execute(contract, registry)
    table = registry.load_payload(out_uri)
    assert result.status == TaskStatus.SUCCESS
    assert set(table["source"]) == {"pathway_A", "tf_B"}
    assert table["p_value"].notna().all()
    assert table["fdr_q_value"].between(0, 1).all()
    assert table["network_source"].eq("local_smoke_network").all()
    assert result.metrics["network_sha256"]
    residual_df = int(table["residual_degrees_of_freedom"].iloc[0])
    critical_t = stats.t.ppf(1.0 - 0.05 / 2.0, residual_df)
    critical_z = stats.norm.ppf(1.0 - 0.05 / 2.0)
    np.testing.assert_allclose(
        table["activity_ci_low"],
        table["activity_effect_condition_a_vs_b"] - critical_t * table["activity_se"],
    )
    np.testing.assert_allclose(
        table["activity_ci_high"],
        table["activity_effect_condition_a_vs_b"] + critical_t * table["activity_se"],
    )
    assert critical_t > critical_z

    from eacbp.auditor import ScientificAuditor
    contract.validation_requirements = ["advanced_statistics_integrity"]
    audited = ScientificAuditor().audit_task(contract, result, registry)
    assert audited.overall_passed, [(c.check_name, c.message) for c in audited.checks if not c.passed]

    original_get = registry.get

    def corrupted_activity(uri):
        metadata, payload = original_get(uri)
        if uri == out_uri:
            payload = payload.copy()
            payload["activity_ci_low"] = payload["activity_effect_condition_a_vs_b"] - critical_z * payload["activity_se"]
            payload["activity_ci_high"] = payload["activity_effect_condition_a_vs_b"] + critical_z * payload["activity_se"]
        return metadata, payload

    from unittest.mock import patch
    with patch.object(registry, "get", side_effect=corrupted_activity):
        tampered = ScientificAuditor().audit_task(contract, result, registry)
    assert not tampered.overall_passed
    assert any(
        c.check_name == "advanced_statistics_result_statistics"
        and not c.passed
        and "activity" in c.message
        for c in tampered.checks
    )


def test_activity_ci_audit_rejects_missing_or_invalid_residual_degrees_of_freedom():
    frame = pd.DataFrame(
        {
            "activity_ci_low": [-1.0],
            "activity_ci_high": [1.0],
            "activity_se": [0.5],
            "activity_effect_condition_a_vs_b": [0.0],
        }
    )
    assert any("residual degrees of freedom" in error for error in _finite_ci_errors(frame, 0.05))

    for invalid_df in (0.0, -1.0, np.nan, np.inf):
        invalid = frame.assign(residual_degrees_of_freedom=invalid_df)
        errors = _finite_ci_errors(invalid, 0.05)
        assert any("finite positive residual degrees of freedom" in error for error in errors)


def test_decoupler_rejects_changed_local_network_hash(tmp_path):
    registry, uri = _registry(tmp_path, _data())
    network_path = tmp_path / "network.csv"
    network = pd.DataFrame({"source": ["s"], "target": ["g0"], "weight": [1.0]})
    network.to_csv(network_path, index=False)
    contract = TaskContract(
        task_id="advanced_activity_bad_hash",
        capability="functional_activity",
        input_artifacts=[uri],
        parameters={
            "condition_a": "A",
            "condition_b": "B",
            "donor_col": "donor_id",
            "paired": True,
            "network_path": str(network_path),
            "network_sha256": "0" * 64,
            "species": "Mus musculus",
            "network_source": "local_smoke_network",
            "network_version": "1.0",
            "tmin": 1,
        },
    )
    with pytest.raises(AdvancedStatisticsInputError, match="hash mismatch"):
        DecouplerFunctionalAnalysisCapability().execute(contract, registry)


def test_leave_one_donor_out_real_refits_and_skips_when_underpowered(tmp_path):
    registry, uri = _registry(tmp_path, _data(n_pairs=8))
    out_uri = "table://advanced/sensitivity_expected/v1"
    summary_uri = "table://advanced/sensitivity_summary_expected/v1"
    contract = TaskContract(
        task_id="advanced_sensitivity",
        capability="donor_sensitivity",
        method="pydeseq2_leave_one_donor_out_v1",
        input_artifacts=[uri],
        expected_outputs=[out_uri, summary_uri],
        parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True, "n_cpus": 1},
    )
    result = PyDESeq2LeaveOneDonorOutCapability().execute(contract, registry)
    table = registry.load_payload(out_uri)
    assert result.metrics["skipped"] is False
    assert result.metrics["n_successful_fits"] == 8
    assert {"__full_model__", *(f"d{i}" for i in range(1, 9))}.issubset(set(table["left_out_donor"]))
    assert result.metrics["paired_leaveout_removes_complete_donor"] is True
    full_rows = table.loc[table["left_out_donor"] == "__full_model__"]
    assert set(full_rows["n_pseudobulk_samples"]) == {16}
    assert set(full_rows["n_donors_condition_a"]) == {8}
    assert set(full_rows["n_donors_condition_b"]) == {8}
    for donor in (f"d{i}" for i in range(1, 9)):
        donor_rows = table.loc[table["left_out_donor"] == donor]
        assert set(donor_rows["n_pseudobulk_samples"]) == {14}
        assert set(donor_rows["n_donors_condition_a"]) == {7}
        assert set(donor_rows["n_donors_condition_b"]) == {7}

    validator = AdvancedStatisticsValidator()
    report = validator.audit(contract, result, registry)
    loo_check = next(check for check in report.checks if check.check_name == "advanced_statistics_loo_coverage")
    assert loo_check.passed, loo_check.message
    rank_check = next(check for check in report.checks if check.check_name == "advanced_statistics_design_rank")
    assert rank_check.passed, rank_check.message
    assert result.metrics["design_rank"] == 9
    assert result.metrics["residual_degrees_of_freedom"] == 7
    assert set(full_rows["design_rank"]) == {9}
    assert set(full_rows["residual_degrees_of_freedom"]) == {7}
    for donor in (f"d{i}" for i in range(1, 9)):
        donor_rows = table.loc[table["left_out_donor"] == donor]
        assert set(donor_rows["design_rank"]) == {8}
        assert set(donor_rows["residual_degrees_of_freedom"]) == {6}
    for check_name in ("multiple_testing_correction", "pseudoreplication_audit"):
        check = next(item for item in report.checks if item.check_name == check_name)
        assert check.passed, check.message
    pseudoreplication_check = next(item for item in report.checks if item.check_name == "pseudoreplication_audit")
    assert pseudoreplication_check.metrics["donor_independence"] == "unknown"
    assert pseudoreplication_check.metrics["confirmatory_allowed"] is False

    class TamperedPrimaryRegistry:
        def __init__(self, field, value):
            self.field = field
            self.value = value

        def __getattr__(self, name):
            return getattr(registry, name)

        def get(self, uri):
            metadata, payload = registry.get(uri)
            if uri == out_uri:
                payload = payload.copy()
                mask = payload["left_out_donor"].astype(str) == "d1"
                payload.loc[mask, self.field] = self.value
            return metadata, payload

    for field, value in (("n_pseudobulk_samples", 16), ("n_donors_condition_a", 8)):
        altered = validator.audit(contract, result, TamperedPrimaryRegistry(field, value))
        altered_loo_check = next(
            check for check in altered.checks if check.check_name == "advanced_statistics_loo_coverage"
        )
        assert not altered_loo_check.passed
        assert field in altered_loo_check.message

    for field, value in (("design_rank", 9), ("residual_degrees_of_freedom", 7)):
        altered = validator.audit(contract, result, TamperedPrimaryRegistry(field, value))
        altered_rank_check = next(
            check for check in altered.checks if check.check_name == "advanced_statistics_design_rank"
        )
        assert not altered_rank_check.passed
        assert "LOO fit d1" in altered_rank_check.message

    summary = registry.load_payload(summary_uri)
    assert len(summary) == 6
    assert {"estimated_coverage", "direction_consistency", "significance_retention"}.issubset(summary.columns)
    assert result.metrics["scientific_robustness_claim_supported"] is False

    small_registry, small_uri = _registry(tmp_path / "small", _data(n_pairs=2))
    small_out = "table://advanced/sensitivity_skipped/v1"
    small_summary = "table://advanced/sensitivity_skipped_summary/v1"
    small_contract = TaskContract(
        task_id="advanced_sensitivity_skip",
        capability="donor_sensitivity",
        input_artifacts=[small_uri],
        expected_outputs=[small_out, small_summary],
        parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True, "n_cpus": 1},
    )
    small_result = PyDESeq2LeaveOneDonorOutCapability().execute(small_contract, small_registry)
    assert small_result.metrics["skipped"] is True
    assert "at least 3 donors" in small_result.metrics["skip_reason"]
    assert small_registry.load_payload(small_out).iloc[0]["status"] == "skipped"
    assert small_result.output_artifacts == [small_out, small_summary]
    assert small_registry.load_payload(small_summary).iloc[0]["status"] == "skipped"


def test_leave_one_donor_out_summary_and_audit_keep_nullable_partial_fit_semantics():
    from types import SimpleNamespace

    from eacbp.auditor.advanced_statistics import _check_summary_against_fits
    from eacbp.capabilities.advanced_statistics import _leave_one_out_summary

    full_table = pd.DataFrame({
        "gene": ["g_a", "g_b", "g_c"],
        "condition_a": ["stim"] * 3,
        "condition_b": ["ctrl"] * 3,
        "contrast_kind": ["categorical"] * 3,
        "contrast_factor": ["condition"] * 3,
        "contrast_tested_level": ["stim"] * 3,
        "contrast_reference_level": ["ctrl"] * 3,
        "contrast_label": ["stim versus ctrl"] * 3,
        "effect_definition": ["PyDESeq2 log2 fold change"] * 3,
        "log2_fold_change": [2.0, 0.0, np.nan],
        "fdr_q_value": [0.08, 0.4, np.nan],
    })
    fit_tables = [
        pd.DataFrame({
            "gene": ["g_a", "g_b", "g_c"],
            "left_out_donor": ["d1"] * 3,
            "log2_fold_change": [1.1, -0.2, np.nan],
            "fdr_q_value": [0.09, 0.6, np.nan],
        }),
        pd.DataFrame({
            "gene": ["g_a", "g_b", "g_c"],
            "left_out_donor": ["d2"] * 3,
            "log2_fold_change": [np.nan, 0.3, 1.2],
            "fdr_q_value": [np.nan, 0.03, np.nan],
        }),
    ]

    # Only two of three requested donor fits are present. Missing feature
    # estimates remain missing and coverage keeps the requested-fit denominator.
    summary = _leave_one_out_summary(full_table, fit_tables, n_requested_fits=3, alpha=0.1)
    assert summary.columns.tolist() == [
        "gene", "condition_a", "condition_b", "contrast_kind", "contrast_factor",
        "contrast_tested_level", "contrast_reference_level", "contrast_label",
        "effect_definition", "alpha", "confidence_level", "baseline_log2_fold_change",
        "baseline_fdr_q_value", "baseline_significant_fdr05", "baseline_significant_at_alpha",
        "n_estimated_loo", "estimated_coverage", "direction_consistency",
        "min_log2_fold_change", "max_log2_fold_change", "median_log2_fold_change",
        "n_significant_loo", "n_significant_loo_fdr05", "significance_retention",
        "summary_status",
    ]
    by_gene = summary.set_index("gene")
    assert by_gene.loc["g_a", "n_estimated_loo"] == 1
    assert np.isclose(by_gene.loc["g_a", "estimated_coverage"], 1 / 3)
    assert not bool(by_gene.loc["g_a", "baseline_significant_fdr05"])
    assert bool(by_gene.loc["g_a", "baseline_significant_at_alpha"])
    assert by_gene.loc["g_a", "n_significant_loo"] == 1
    assert by_gene.loc["g_a", "n_significant_loo_fdr05"] == 0
    assert np.isnan(by_gene.loc["g_b", "direction_consistency"])
    assert np.isclose(by_gene.loc["g_b", "estimated_coverage"], 2 / 3)
    assert by_gene.loc["g_b", "n_significant_loo"] == 1
    assert by_gene.loc["g_b", "n_significant_loo_fdr05"] == 1
    assert np.isnan(by_gene.loc["g_c", "direction_consistency"])
    assert np.isnan(by_gene.loc["g_c", "significance_retention"])
    assert pd.isna(by_gene.loc["g_c", "baseline_significant_at_alpha"])
    assert pd.isna(by_gene.loc["g_c", "baseline_significant_fdr05"])

    primary = full_table[["gene", "log2_fold_change", "fdr_q_value"]].copy()
    primary["left_out_donor"] = "__full_model__"
    primary["n_requested_fits"] = 3
    primary = pd.concat(
        [primary, *(fit.assign(n_requested_fits=3) for fit in fit_tables)],
        ignore_index=True,
    )
    details = SimpleNamespace(genes=["g_a", "g_b", "g_c"], paired=True)
    assert _check_summary_against_fits(summary, primary, details, alpha=0.1) == []

    tampered = summary.copy()
    tampered.loc[tampered["gene"] == "g_a", "estimated_coverage"] = 0.5
    errors = _check_summary_against_fits(tampered, primary, details, alpha=0.1)
    assert "gene g_a has incorrect estimated coverage" in errors


def test_loo_significance_count_survives_na_baseline_and_artifact_round_trip(tmp_path):
    from types import SimpleNamespace

    from eacbp.auditor.advanced_statistics import _check_summary_against_fits
    from eacbp.capabilities.advanced_statistics import _leave_one_out_summary

    full_table = pd.DataFrame({
        "gene": ["g_na"],
        "condition_a": ["stim"],
        "condition_b": ["ctrl"],
        "log2_fold_change": [1.5],
        "fdr_q_value": [np.nan],
    })
    fit_tables = [
        pd.DataFrame({
            "gene": ["g_na"],
            "left_out_donor": ["d1"],
            "log2_fold_change": [1.7],
            "fdr_q_value": [0.01],
        }),
        pd.DataFrame({
            "gene": ["g_na"],
            "left_out_donor": ["d2"],
            "log2_fold_change": [1.3],
            "fdr_q_value": [np.nan],
        }),
    ]
    summary = _leave_one_out_summary(full_table, fit_tables, n_requested_fits=2, alpha=0.05)
    assert summary.loc[0, "n_significant_loo"] == 1
    assert summary.loc[0, "n_significant_loo_fdr05"] == 1
    assert pd.isna(summary.loc[0, "significance_retention"])
    assert pd.isna(summary.loc[0, "baseline_significant_at_alpha"])

    primary = pd.DataFrame({
        "gene": ["g_na"] * 3,
        "left_out_donor": ["__full_model__", "d1", "d2"],
        "log2_fold_change": [1.5, 1.7, 1.3],
        "fdr_q_value": [np.nan, 0.01, np.nan],
        "n_requested_fits": [2] * 3,
    })
    registry = ArtifactRegistry(str(tmp_path / "persisted_loo_na"))
    registry.register(
        "table://loo_na/primary/v1", primary, ArtifactType.TABLE, "loo_na", "test", "persist_primary"
    )
    registry.register(
        "table://loo_na/summary/v1", summary, ArtifactType.TABLE, "loo_na", "test", "persist_summary"
    )
    stored_primary = registry.load_payload("table://loo_na/primary/v1")
    stored_summary = registry.load_payload("table://loo_na/summary/v1")
    assert pd.isna(stored_primary.loc[0, "fdr_q_value"])
    assert pd.isna(stored_primary.loc[2, "fdr_q_value"])
    assert pd.isna(stored_summary.loc[0, "baseline_fdr_q_value"])
    assert stored_summary.loc[0, "n_significant_loo"] == 1
    assert stored_summary.loc[0, "n_significant_loo_fdr05"] == 1
    assert pd.isna(stored_summary.loc[0, "significance_retention"])
    errors = _check_summary_against_fits(
        stored_summary, stored_primary, SimpleNamespace(genes=["g_na"]), alpha=0.05
    )
    assert errors == []

    corrupted = stored_summary.copy()
    corrupted.loc[0, "n_significant_loo"] = 0
    errors = _check_summary_against_fits(
        corrupted, stored_primary, SimpleNamespace(genes=["g_na"]), alpha=0.05
    )
    assert "gene g_na has incorrect n_significant_loo" in errors

    corrupted = stored_summary.copy()
    corrupted.loc[0, "significance_retention"] = 1.0
    errors = _check_summary_against_fits(
        corrupted, stored_primary, SimpleNamespace(genes=["g_na"]), alpha=0.05
    )
    assert "gene g_na has unsupported significance retention" in errors


def test_advanced_plan_declares_donor_sensitivity_summary_output():
    from eacbp.orchestrator.advanced_plan import extend_plan

    manifest = StudyManifest(
        study_id="planner_contract",
        biological_design=BiologicalDesign(species="human", tissue="test"),
        data=DataSpec(raw_artifact_uri="adata://planner_contract/raw/v1"),
    )
    tasks = [
        TaskContract(
            task_id="deg",
            capability="deg",
            method="welch_ttest_v1",
            input_artifacts=["adata://planner_contract/raw/v1"],
            parameters={"condition_a": "A", "condition_b": "B", "donor_col": "donor"},
        ),
        TaskContract(task_id="audit", capability="dataset_audit", input_artifacts=["adata://planner_contract/raw/v1"]),
    ]
    planned = extend_plan(tasks, manifest, {"advanced_analysis": True, "analysis_extensions": {}})
    loo = next(task for task in planned if task.capability == "donor_sensitivity")
    assert loo.expected_outputs == [
        "table://planner_contract/donor_sensitivity/v1",
        "table://planner_contract/donor_sensitivity_summary/v1",
    ]
