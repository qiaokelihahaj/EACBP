from copy import deepcopy
from types import SimpleNamespace

import pandas as pd
import pytest

from eacbp.auditor.advanced_statistics import _check_summary_against_fits, _inference_errors


def _fixture():
    details = SimpleNamespace(condition_a="A", condition_b="B", condition_col="condition")
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
