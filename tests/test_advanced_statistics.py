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
    registry, uri = _registry(tmp_path, _data(n_pairs=3))
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
    assert result.metrics["n_successful_fits"] == 3
    assert {"__full_model__", "d1", "d2", "d3"}.issubset(set(table["left_out_donor"]))
    assert result.metrics["paired_leaveout_removes_complete_donor"] is True
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
