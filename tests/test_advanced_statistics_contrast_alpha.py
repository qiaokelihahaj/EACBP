import numpy as np
import pandas as pd
import pytest
from scipy import stats

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport
from eacbp.auditor import ScientificAuditor
from eacbp.capabilities.advanced_statistics import (
    AdvancedStatisticsInputError,
    PyDESeq2PseudobulkCapability,
    _prepare_pseudobulk,
    _resolve_contrast,
    _validate_design_and_make_dds,
)
from eacbp.capabilities.sc_data import SCData
from eacbp.evidence.claim import ClaimEngine
from eacbp.evidence.extraction import extract_evidence
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.evidence import EvidenceNode, EvidenceType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _data(n_donors=4):
    rng = np.random.default_rng(1907)
    rows, metadata = [], []
    for donor_index in range(n_donors):
        donor = f"d{donor_index + 1}"
        for condition in ("A", "B"):
            center = np.array(
                [90 if condition == "A" else 10,
                 10 if condition == "A" else 90,
                 35, 35, 35, 35],
                dtype=int,
            )
            for cell_index in range(4):
                rows.append(rng.poisson(center).astype(np.int64))
                metadata.append(
                    {
                        "cell_id": f"{donor}_{condition}_{cell_index}",
                        "condition": condition,
                        "donor_id": donor,
                    }
                )
    counts = np.asarray(rows, dtype=np.int64)
    return SCData(
        X=counts.astype(np.float32),
        obs=pd.DataFrame(metadata),
        var=pd.DataFrame({"gene_name": [f"g{i}" for i in range(counts.shape[1])]}),
        layers={"counts": counts},
    )


def _registered_data(tmp_path):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://contrast/input/v1"
    registry.register(input_uri, _data(), ArtifactType.ANNDATA, "contrast", "input", "input")
    return registry, input_uri


def _run(registry, input_uri, task_id, output_uri, params):
    contract = TaskContract(
        task_id=task_id,
        capability="deg",
        method="pydeseq2_pseudobulk_v1",
        input_artifacts=[input_uri],
        expected_outputs=[output_uri],
        parameters=params,
    )
    result = PyDESeq2PseudobulkCapability().execute(contract, registry)
    return contract, result, registry.load_payload(output_uri)


def _base_params(**overrides):
    return {
        "condition_a": "A",
        "condition_b": "B",
        "donor_col": "donor_id",
        "paired": True,
        "n_cpus": 1,
        "independent_filter": False,
        "cooks_filter": False,
        **overrides,
    }


def test_reverse_contrast_reverses_effect_ci_and_result_labels(tmp_path):
    pytest.importorskip("pydeseq2")
    registry, input_uri = _registered_data(tmp_path)
    _, forward_result, forward = _run(
        registry, input_uri, "contrast_a_vs_b", "table://contrast/a_vs_b/v1", _base_params()
    )
    _, reverse_result, reverse = _run(
        registry,
        input_uri,
        "contrast_b_vs_a",
        "table://contrast/b_vs_a/v1",
        _base_params(contrast=["condition", "B", "A"]),
    )

    assert forward_result.status == reverse_result.status == TaskStatus.SUCCESS
    assert reverse["condition_a"].eq("B").all()
    assert reverse["condition_b"].eq("A").all()
    assert reverse["contrast_label"].eq("B versus A").all()
    assert reverse["effect_definition"].eq("PyDESeq2 log2 fold change for B versus A").all()
    assert reverse_result.metrics["condition_a"] == "B"
    assert reverse_result.metrics["condition_b"] == "A"
    assert reverse_result.metrics["contrast_label"] == "B versus A"

    np.testing.assert_allclose(
        reverse["log2_fold_change"], -forward["log2_fold_change"], rtol=1e-6, atol=1e-8, equal_nan=True
    )
    np.testing.assert_allclose(reverse["lfc_se"], forward["lfc_se"], rtol=1e-6, atol=1e-8, equal_nan=True)
    np.testing.assert_allclose(reverse["ci_low"], -forward["ci_high"], rtol=1e-6, atol=1e-8, equal_nan=True)
    np.testing.assert_allclose(reverse["ci_high"], -forward["ci_low"], rtol=1e-6, atol=1e-8, equal_nan=True)
    np.testing.assert_allclose(reverse["p_value"], forward["p_value"], rtol=1e-6, atol=1e-8, equal_nan=True)


def test_nondefault_alpha_controls_ci_and_keeps_legacy_fdr05_meaning(tmp_path):
    pytest.importorskip("pydeseq2")
    alpha = 0.2
    registry, input_uri = _registered_data(tmp_path)
    _, result, table = _run(
        registry,
        input_uri,
        "contrast_alpha_020",
        "table://contrast/alpha_020/v1",
        _base_params(alpha=alpha),
    )

    assert result.metrics["alpha"] == alpha
    assert result.metrics["inference_settings"]["confidence_level"] == pytest.approx(1 - alpha)
    assert table["alpha"].eq(alpha).all()
    assert table["confidence_level"].eq(1 - alpha).all()
    critical = stats.norm.ppf(1 - alpha / 2)
    np.testing.assert_allclose(
        table["ci_low"], table["log2_fold_change"] - critical * table["lfc_se"], rtol=1e-6, atol=1e-8
    )
    np.testing.assert_allclose(
        table["ci_high"], table["log2_fold_change"] + critical * table["lfc_se"], rtol=1e-6, atol=1e-8
    )
    valid = table["fdr_q_value"].notna()
    assert table.loc[valid, "significant_at_alpha"].astype(bool).tolist() == (
        table.loc[valid, "fdr_q_value"].lt(alpha).tolist()
    )
    assert table.loc[valid, "significant_fdr05"].astype(bool).tolist() == (
        table.loc[valid, "fdr_q_value"].lt(0.05).tolist()
    )
    assert result.metrics["n_significant_at_alpha"] == int(table["fdr_q_value"].lt(alpha).sum())
    assert result.metrics["n_significant_fdr05"] == int(table["fdr_q_value"].lt(0.05).sum())


def test_numeric_contrast_requires_an_effect_definition():
    data = _data()
    params = _base_params()
    prepared = _prepare_pseudobulk(data, params)
    with pytest.raises(AdvancedStatisticsInputError, match="contrast_effect_definition"):
        _resolve_contrast(prepared, {**params, "contrast": [1.0, 0.0, 0.0, 0.0]})


def test_numeric_contrast_with_definition_runs_and_labels_its_estimand(tmp_path):
    pytest.importorskip("pydeseq2")
    registry, input_uri = _registered_data(tmp_path)
    params = _base_params()
    prepared = _prepare_pseudobulk(_data(), params)
    _validate_design_and_make_dds(prepared, params)
    coefficient_names = prepared.design_columns
    condition_coefficient = next(
        index
        for index, name in enumerate(coefficient_names)
        if "condition" in name and "B" in name
    )
    vector = [0.0] * len(coefficient_names)
    vector[condition_coefficient] = 1.0
    definition = "Condition B log2 fold-change coefficient relative to condition A"
    contract, result, table = _run(
        registry,
        input_uri,
        "contrast_numeric_defined",
        "table://contrast/numeric_defined/v1",
        {
            **params,
            "contrast": vector,
            "contrast_effect_definition": definition,
            "contrast_label": "Condition B coefficient",
        },
    )

    assert result.status == TaskStatus.SUCCESS
    assert table["contrast_kind"].eq("numeric").all()
    assert table["contrast_label"].eq("Condition B coefficient").all()
    assert table["effect_definition"].eq(definition).all()
    assert table["contrast_tested_level"].isna().all()
    assert table["contrast_reference_level"].isna().all()
    assert result.metrics["contrast_spec"]["kind"] == "numeric"
    assert result.metrics["contrast_spec"]["effect_definition"] == definition
    contract.validation_requirements = ["advanced_statistics_integrity"]
    audit = ScientificAuditor().audit_task(contract, result, registry)
    assert audit.overall_passed, [(check.check_name, check.message) for check in audit.checks if not check.passed]


def test_custom_condition_factor_and_formula_preserve_reverse_direction(tmp_path):
    pytest.importorskip("pydeseq2")
    data = _data()
    data.obs["treatment"] = data.obs["condition"]
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://contrast/custom_condition/v1"
    registry.register(input_uri, data, ArtifactType.ANNDATA, "contrast", "input", "input")
    params = _base_params(
        condition_col="treatment",
        design_formula="~treatment + donor",
        contrast=["treatment", "B", "A"],
    )
    contract, result, table = _run(
        registry,
        input_uri,
        "contrast_custom_factor",
        "table://contrast/custom_factor/v1",
        params,
    )

    assert result.status == TaskStatus.SUCCESS
    assert result.metrics["contrast_spec"]["factor"] == "treatment"
    assert table["contrast_factor"].eq("treatment").all()
    assert table["condition_a"].eq("B").all()
    assert table["condition_b"].eq("A").all()
    assert table["log2_fold_change"].loc[table["gene"] == "g0"].iloc[0] < 0
    contract.validation_requirements = ["advanced_statistics_integrity"]
    audit = ScientificAuditor().audit_task(contract, result, registry)
    assert audit.overall_passed, [(check.check_name, check.message) for check in audit.checks if not check.passed]


@pytest.mark.parametrize("alpha", [0.0, 1.0, float("nan"), "invalid"])
def test_claim_engine_fails_closed_for_invalid_explicit_alpha(alpha):
    node = EvidenceNode(
        evidence_id="E",
        type=EvidenceType.PSEUDOBULK_DEG,
        summary="Test evidence",
        source_task_id="T",
        source_artifact_uris=["table://study/deg/v1"],
        audit_passed=True,
        metrics={"fdr_q_value": 0.01, "alpha": alpha},
    )
    assert not ClaimEngine.has_valid_statistics(node)


def test_evidence_extraction_uses_recorded_alpha_and_contrast_label(tmp_path):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    output_uri = "table://evidence/deg/v1"
    registry.register(
        output_uri,
        pd.DataFrame(
            [
                {
                    "gene": "GeneA",
                    "condition_a": "B",
                    "condition_b": "A",
                    "contrast_kind": "categorical",
                    "contrast_label": "B versus A",
                    "effect_definition": "PyDESeq2 log2 fold change for B versus A",
                    "alpha": 0.1,
                    "log2_fold_change": 1.2,
                    "fdr_q_value": 0.08,
                }
            ]
        ),
        ArtifactType.TABLE,
        "evidence",
        "deg",
        "differential_expression_pydeseq2_donor_pseudobulk",
        summary_metrics={"alpha": 0.1},
    )
    contract = TaskContract(task_id="evidence_alpha", capability="deg", parameters={"alpha": 0.1})
    result = TaskResult(
        task_id=contract.task_id,
        capability="deg",
        method_used="pydeseq2_pseudobulk_v1",
        status=TaskStatus.SUCCESS,
        output_artifacts=[output_uri],
        metrics={"alpha": 0.1, "is_pseudobulk": True},
    )
    report = ValidationReport(
        auditor_name="test",
        target_task_id=contract.task_id,
        overall_passed=True,
        stop_rule_triggered=False,
    )

    evidence = extract_evidence(contract, result, report, registry)
    assert len(evidence) == 1
    assert "B versus A" in evidence[0].summary
    assert "alpha=0.1" in evidence[0].summary
    assert ClaimEngine.has_valid_statistics(evidence[0])
