from copy import deepcopy
import sys
from types import ModuleType
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from eacbp.auditor.advanced_statistics import _check_summary_against_fits, _inference_errors


def _fixture():
    details = SimpleNamespace(condition_a="A", condition_b="B", condition_col="condition", formula="~condition")
    spec = dict(kind="categorical", vector=["condition", "B", "A"], factor="condition",
                tested_level="B", reference_level="A", label="B versus A",
                effect_definition="PyDESeq2 log2 fold change for B versus A")
    metrics = dict(contrast_spec=spec, contrast_label=spec["label"], effect_definition=spec["effect_definition"],
                   alpha=0.1, inference_settings={"alpha": 0.1, "confidence_level": 0.9})
    table = pd.DataFrame([dict(condition_a="B", condition_b="A", contrast_kind="categorical",
                              contrast_factor="condition", contrast_tested_level="B", contrast_reference_level="A",
                              contrast_label=spec["label"], effect_definition=spec["effect_definition"],
                              alpha=0.1, confidence_level=0.9, fdr_q_value=0.08,
                              log2_fold_change=1.0, lfc_se=0.5, ci_low=0.1775731865, ci_high=1.8224268135,
                              significant_at_alpha=True, significant_fdr05=False)])
    params = {"alpha": 0.1, "contrast": ["condition", "B", "A"]}
    return table, metrics, params, details


@pytest.mark.parametrize("field,value", [("condition_a", "A"), ("alpha", 0.05),
                                            ("confidence_level", 0.95), ("significant_at_alpha", False),
                                            ("significant_fdr05", True), ("effect_definition", "A versus B")])
def test_auditor_rejects_result_semantics_that_disagree_with_contract(field, value):
    table, metrics, params, details = _fixture()
    assert not _inference_errors([table], metrics, params, details, "pydeseq2_pseudobulk_v1")
    changed = table.copy()
    changed[field] = value
    assert _inference_errors([changed], metrics, params, details, "pydeseq2_pseudobulk_v1")


def test_auditor_does_not_trust_reported_contrast_spec():
    table, metrics, params, details = _fixture()
    metrics = deepcopy(metrics)
    metrics["contrast_spec"]["tested_level"] = "A"
    assert _inference_errors([table], metrics, params, details, "pydeseq2_pseudobulk_v1")


def test_auditor_rejects_deleted_uncertainty_columns():
    table, metrics, params, details = _fixture()
    table = table.drop(columns=["lfc_se", "ci_low", "ci_high"])
    assert _inference_errors([table], metrics, params, details, "pydeseq2_pseudobulk_v1")


@pytest.mark.parametrize("alpha", ["invalid", float("nan"), 0.0, 1.0])
def test_auditor_returns_failure_report_for_invalid_alpha(alpha):
    from eacbp.auditor.advanced_statistics import AdvancedStatisticsValidator
    from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus

    contract = TaskContract(task_id="invalid_alpha", capability="deg", parameters={"alpha": alpha})
    result = TaskResult(task_id=contract.task_id, capability="deg", method_used="pydeseq2_pseudobulk_v1", status=TaskStatus.SUCCESS)
    report = AdvancedStatisticsValidator().audit(contract, result, None)
    assert not report.overall_passed
    assert any("alpha" in check.message for check in report.checks if not check.passed)


def test_loo_auditor_uses_contract_alpha_for_decision_retention():
    primary = pd.DataFrame({"gene": ["g"] * 3, "left_out_donor": ["__full_model__", "d1", "d2"],
                            "log2_fold_change": [1.0, 1.0, 2.0], "fdr_q_value": [0.08, 0.03, 0.07],
                            "n_requested_fits": [2] * 3})
    summary = pd.DataFrame([dict(gene="g", baseline_log2_fold_change=1.0, baseline_fdr_q_value=0.08,
                                 n_estimated_loo=2, estimated_coverage=1.0, direction_consistency=1.0,
                                 min_log2_fold_change=1.0, max_log2_fold_change=2.0, median_log2_fold_change=1.5,
                                 significance_retention=1.0, n_significant_loo=2, n_significant_loo_fdr05=1,
                                 baseline_significant_at_alpha=True, baseline_significant_fdr05=False)])
    details = SimpleNamespace(genes=["g"])
    assert not _check_summary_against_fits(summary, primary, details, 0.1)
    assert _check_summary_against_fits(summary, primary, details, 0.05)


def _loo_audit_case(tmp_path, monkeypatch, tamper=None):
    """Build a small persisted LOO result without running the optional fitter."""
    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.auditor import ScientificAuditor
    from eacbp.capabilities.advanced_statistics import (
        InferenceSettings,
        _build_inference_contract,
        _leave_one_out_summary,
        _record_scientific_result,
        _register_table,
        _resolve_contrast,
        _result_table,
        _scientific_result,
        _prepare_pseudobulk,
    )
    from eacbp.capabilities.sc_data import SCData
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus

    def model_matrix(formula, frame):
        terms = [term.strip() for term in formula.removeprefix("~").split("+")]
        matrix = pd.DataFrame({"Intercept": np.ones(len(frame))}, index=frame.index)
        for factor in terms:
            if factor == "1":
                continue
            if factor not in frame:
                raise ValueError(f"unsupported fixture design factor: {factor}")
            levels = sorted(frame[factor].astype(str).unique())
            for level in levels[1:]:
                matrix[f"{factor}[T.{level}]"] = frame[factor].astype(str).eq(level).astype(float)
        return matrix

    formulaic = ModuleType("formulaic")
    formulaic.model_matrix = model_matrix
    monkeypatch.setitem(sys.modules, "formulaic", formulaic)

    obs_rows = []
    count_rows = []
    for donor_index in range(1, 9):
        for condition in ("ctrl", "stim"):
            for cell_index in range(2):
                obs_rows.append({"donor_id": f"d{donor_index}", "condition": condition})
                count_rows.append([12 + donor_index, 18 + cell_index, 7 + donor_index + cell_index])
    counts = np.asarray(count_rows, dtype=np.int64)
    data = SCData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(obs_rows),
        var=pd.DataFrame({"gene_name": ["g1", "g2", "g_filtered"]}),
        layers={"counts": counts},
    )
    input_uri = "adata://loo_contract/input/v1"
    primary_uri = "table://loo_contract/primary/v1"
    summary_uri = "table://loo_contract/summary/v1"
    registry = ArtifactRegistry(str(tmp_path / "registry"))
    registry.register(
        input_uri, data, ArtifactType.ANNDATA, "loo_contract", "fixture", "input"
    )
    params = {
        "condition_a": "stim",
        "condition_b": "ctrl",
        "donor_col": "donor_id",
        "paired": True,
        "counts_layer": "counts",
        "alpha": 0.05,
        "contrast": ["condition", "stim", "ctrl"],
        "min_loo_donors": 3,
    }
    contract = TaskContract(
        task_id="loo_contract",
        capability="donor_sensitivity",
        method="pydeseq2_leave_one_donor_out_v1",
        input_artifacts=[input_uri],
        expected_outputs=[primary_uri, summary_uri],
        parameters=params,
        validation_requirements=[
            "advanced_statistics_integrity",
            "multiple_testing_correction",
            "pseudoreplication_audit",
        ],
    )
    contract.inference_contract = _build_inference_contract(contract, registry, data)
    prepared = _prepare_pseudobulk(data, params)
    inference = InferenceSettings.from_params(params)
    contrast = _resolve_contrast(prepared, params)
    genes = ["g1", "g2", "g_filtered"]
    raw_results = pd.DataFrame(
        {
            "baseMean": [30.0, 24.0, 12.0],
            "log2FoldChange": [1.0, 0.2, 0.4],
            "lfcSE": [0.25, 0.25, 0.25],
            "stat": [4.0, 0.8, 1.6],
            "pvalue": [0.01, 0.2, 0.8],
            # Independent filtering leaves the last adjusted p-value unknown.
            "padj": [0.02, 0.2, np.nan],
        },
        index=genes,
    )

    def fit_table(fit_id, fit_metadata, rank, residual_df):
        table = _result_table(raw_results, prepared, params, fit_metadata=fit_metadata)
        table.insert(1, "left_out_donor", fit_id)
        table["fit_status"] = "full_model" if fit_id == "__full_model__" else "estimated"
        table["removed_sample_count"] = 0 if fit_id == "__full_model__" else 2
        table["design_rank"] = rank
        table["residual_degrees_of_freedom"] = residual_df
        return table

    full = fit_table("__full_model__", prepared.metadata, rank=9, residual_df=7)
    fits = [
        fit_table(
            f"d{donor_index}",
            prepared.metadata.loc[prepared.metadata["donor"] != f"d{donor_index}"],
            rank=8,
            residual_df=6,
        )
        for donor_index in range(1, 9)
    ]
    if tamper == "q_value":
        altered = fits[0]["gene"].eq("g1")
        fits[0].loc[altered, ["fdr_q_value", "padj"]] = 0.99
        fits[0].loc[altered, ["significant_at_alpha", "significant_fdr05"]] = False
    elif tamper == "statistical_unit":
        fits[0]["statistical_unit"] = "cell"
    elif tamper == "sample_count":
        fits[0]["n_pseudobulk_samples"] = 16

    summary = _leave_one_out_summary(full, fits, n_requested_fits=8, alpha=inference.alpha)
    summary_for_merge = summary.drop(
        columns=[
            "condition_a", "condition_b", "contrast_kind", "contrast_tested_level",
            "contrast_factor", "contrast_reference_level", "contrast_label",
            "effect_definition", "alpha", "confidence_level",
        ],
        errors="ignore",
    )
    primary = pd.concat([full, *fits], ignore_index=True).merge(
        summary_for_merge, on="gene", how="left"
    )
    primary["n_requested_fits"] = 8
    primary["n_successful_fits"] = 8
    primary["paired_leaveout_removes_complete_donor"] = True
    primary["complete_fit_coverage"] = 1.0
    primary["all_requested_fits_completed"] = True
    primary["scientific_robustness_claim_supported"] = False

    settings = inference.as_dict()
    metrics = {
        "skipped": False,
        "skip_reason": None,
        "alpha": inference.alpha,
        "inference_settings": settings,
        "contrast_spec": contrast.as_dict(),
        "contrast_label": contrast.label,
        "effect_definition": contrast.effect_definition,
        "condition_a": contrast.tested_level,
        "condition_b": contrast.reference_level,
        "design": prepared.design,
        "design_rank": 9,
        "design_columns": ["Intercept", "condition[T.stim]", *[f"donor[T.d{i}]" for i in range(2, 9)]],
        "residual_degrees_of_freedom": 7,
        "n_requested_fits": 8,
        "n_successful_fits": 8,
        "n_failed_fits": 0,
        "failed_fits": {},
        "complete_fit_coverage": 1.0,
        "all_requested_fits_completed": True,
        "n_donors_condition_a": 8,
        "n_donors_condition_b": 8,
        "n_pseudobulk_samples": 16,
        "paired": True,
        "paired_leaveout_removes_complete_donor": True,
        "scientific_robustness_claim_supported": False,
        "scientific_robustness_limit": "leave-one-donor-out is a sensitivity analysis, not independent validation",
        "n_genes_with_direction_consistency": 3,
        "n_genes_with_significance_retention": 2,
        "n_genes": len(genes),
        "all_genes_retained_per_successful_fit": True,
    }
    scientific_result = _scientific_result(contract, primary, metrics)
    _record_scientific_result(metrics, contract, scientific_result)
    artifact_parameters = {
        "paired": True,
        "condition_a": contrast.tested_level,
        "condition_b": contrast.reference_level,
        "design": prepared.design,
        "contrast_spec": contrast.as_dict(),
        "inference_settings": settings,
    }
    for uri, table in ((primary_uri, primary), (summary_uri, summary)):
        _register_table(
            registry,
            uri,
            table,
            [input_uri],
            contract,
            "test_persisted_loo_result",
            artifact_parameters,
            metrics,
        )
    result = TaskResult(
        task_id=contract.task_id,
        status=TaskStatus.SUCCESS,
        capability=contract.capability,
        method_used=contract.method,
        input_artifacts=[input_uri],
        output_artifacts=[primary_uri, summary_uri],
        metrics=metrics,
        inference_contract=contract.inference_contract,
        scientific_result=scientific_result,
    )

    # Reopen the on-disk registry so this exercises persisted, hash-verified
    # payloads and inference provenance rather than only in-memory objects.
    persisted_registry = ArtifactRegistry(str(tmp_path / "registry"))
    report = ScientificAuditor().audit_task(contract, result, persisted_registry)
    return report, primary, summary


def _named_check(report, name, predicate):
    return next(
        check for check in report.checks
        if check.check_name == name and predicate(check.metrics)
    )


def test_scientific_auditor_accepts_persisted_eight_pair_loo_contract(tmp_path, monkeypatch):
    report, primary, summary = _loo_audit_case(tmp_path, monkeypatch)

    assert report.overall_passed, [(check.check_name, check.message) for check in report.checks if not check.passed]
    integrity = _named_check(report, "advanced_statistics_integrity", lambda _metrics: True)
    assert integrity.passed
    correction = _named_check(report, "multiple_testing_correction", lambda metrics: "fits" in metrics)
    assert correction.passed
    assert correction.metrics["n_fit_groups"] == 9
    pseudoreplication = _named_check(
        report,
        "pseudoreplication_audit",
        lambda metrics: "donor_independence" in metrics,
    )
    assert pseudoreplication.passed
    assert pseudoreplication.metrics["donor_independence"] == "unknown"
    assert pseudoreplication.metrics["confirmatory_allowed"] is False
    filtered = summary.loc[summary["gene"] == "g_filtered"].iloc[0]
    assert filtered["n_estimated_loo"] == 8
    assert filtered["estimated_coverage"] == 1.0
    assert filtered["n_significant_loo"] == 0
    filtered_primary = primary.loc[primary["gene"] == "g_filtered"]
    assert filtered_primary["status"].eq("estimated_fdr_missing").all()
    assert filtered_primary["status_reason"].str.contains("independent_filter").all()


@pytest.mark.parametrize(
    ("tamper", "check_name", "metric_predicate"),
    [
        ("q_value", "multiple_testing_correction", lambda metrics: "fits" in metrics),
        ("statistical_unit", "pseudoreplication_audit", lambda metrics: "donor_independence" in metrics),
        ("sample_count", "pseudoreplication_audit", lambda metrics: "donor_independence" in metrics),
    ],
)
def test_scientific_auditor_rejects_persisted_loo_mutations(
    tmp_path, monkeypatch, tamper, check_name, metric_predicate
):
    report, _primary, _summary = _loo_audit_case(tmp_path, monkeypatch, tamper=tamper)

    check = _named_check(report, check_name, metric_predicate)
    assert not check.passed
    assert not report.overall_passed
