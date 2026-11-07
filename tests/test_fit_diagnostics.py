"""Real fits plus adversarial checks for per-model diagnostic admission."""
from copy import deepcopy
import warnings

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pydeseq2")

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.auditor import ScientificAuditor
from eacbp.capabilities.advanced_statistics import (
    PyDESeq2LeaveOneDonorOutCapability, PyDESeq2PseudobulkCapability,
)
from eacbp.capabilities.fit_diagnostics import (
    FitReceipt, canonical_sha256, fit_receipts_uri,
)
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def _registry(tmp_path, n_pairs=4):
    rng = np.random.default_rng(1907)
    rows, obs = [], []
    for index in range(n_pairs):
        donor = f"d{index + 1}"
        for condition in ("A", "B"):
            center = [90 if condition == "A" else 10, 10 if condition == "A" else 90, *([35] * 10)]
            for cell in range(4):
                rows.append(rng.poisson(center).astype(np.int64))
                obs.append({"cell_id": f"{donor}_{condition}_{cell}", "condition": condition, "donor_id": donor})
    counts = np.asarray(rows)
    data = SCData(X=counts.astype(np.float32), obs=pd.DataFrame(obs),
                  var=pd.DataFrame({"gene_name": [f"g{i}" for i in range(counts.shape[1])]}),
                  layers={"counts": counts})
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    registry.register("adata://diagnostics/input/v1", data, ArtifactType.ANNDATA, "diagnostics", "input", "input")
    return registry


def _contract(*, loo=False, receipt=True):
    primary = "table://diagnostics/loo/v1" if loo else "table://diagnostics/deg/v1"
    params = {"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True, "n_cpus": 1}
    outputs = [primary]
    if receipt:
        params["fit_receipts_uri"] = fit_receipts_uri(primary)
        outputs.append(params["fit_receipts_uri"])
    return TaskContract(task_id="diagnostic_loo" if loo else "diagnostic_deg",
                        capability="donor_sensitivity" if loo else "deg",
                        method="pydeseq2_leave_one_donor_out_v1" if loo else "pydeseq2_pseudobulk_v1",
                        input_artifacts=["adata://diagnostics/input/v1"],
                        expected_outputs=outputs, parameters=params,
                        validation_requirements=["advanced_statistics_integrity"])


def _assert_admitted(contract, result, registry):
    report = ScientificAuditor().audit_task(contract, result, registry)
    assert report.overall_passed, [(check.check_name, check.message) for check in report.checks if not check.passed]
    return report


def test_real_receipt_is_registered_admitted_and_does_not_change_statistics(tmp_path):
    registry = _registry(tmp_path / "new")
    contract = _contract()
    contract.parameters["fit_type"] = "mean"
    staged = TaskArtifactTransaction(registry)
    try:
        result = PyDESeq2PseudobulkCapability().execute(contract, staged)
        report = _assert_admitted(contract, result, staged)
        diagnostics = staged.load_payload(contract.parameters["fit_receipts_uri"])
        fit = diagnostics["fits"][0]
        assert fit["status"] == "completed"
        assert fit["actual_dispersion_trend"] == "mean"
        assert fit["dispersion_trend_source"] == "dds.uns.disp_function_type"
        assert fit["counts"]["frame"]["index"] == fit["metadata"]["frame"]["index"] == fit["design_matrix"]["frame"]["index"]
        assert fit["effective_options"]["fit_type"] == "mean"
        assert fit["effective_options"]["alpha"] == .05
        assert fit["software_versions"]["pydeseq2"] != "unavailable"
        assert any(check.check_name == "advanced_statistics_fit_receipts" and check.passed for check in report.checks)
        # Registered diagnostics and table are both part of the compute commit.
        registry.commit_task(staged, "diagnostic_commit", result)
        assert registry.exists(contract.parameters["fit_receipts_uri"])
        registry.begin_audit("diagnostic_commit", contract, result, auditor_name=report.auditor_name)
        report = _assert_admitted(contract, result, registry)
        registry.record_audit("diagnostic_commit", contract, result, report, auditor_name=report.auditor_name)
        admitted = registry.get_audited_many("diagnostic_commit", contract)
        assert {metadata.uri for metadata, _ in admitted} == set(contract.expected_outputs)
        enabled = registry.load_payload(contract.expected_outputs[0])
    finally:
        staged.close()

    legacy_registry = _registry(tmp_path / "legacy")
    legacy_contract = _contract(receipt=False)
    legacy_contract.parameters["fit_type"] = "mean"
    legacy = PyDESeq2PseudobulkCapability().execute(legacy_contract, legacy_registry)
    legacy_report = _assert_admitted(legacy_contract, legacy, legacy_registry)
    assert legacy.metrics["fit_receipts_status"] == "unavailable_legacy_task"
    assert any("unavailable" in check.message for check in legacy_report.checks if check.check_name == "advanced_statistics_fit_receipts")
    pd.testing.assert_frame_equal(enabled, legacy_registry.load_payload(legacy_contract.expected_outputs[0]))


def test_real_loo_warning_belongs_only_to_designated_fold_and_failed_fold_is_retained(tmp_path, monkeypatch):
    from pydeseq2.dds import DeseqDataSet
    original = DeseqDataSet.deseq2

    def controlled(self, *args, **kwargs):
        observed = set(self.obs.donor.astype(str))
        if "d2" not in observed:
            warnings.warn("specified donor d2 diagnostic warning", UserWarning)
        if "d3" not in observed:
            raise RuntimeError("specified donor d3 diagnostic failure")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(DeseqDataSet, "deseq2", controlled)
    registry = _registry(tmp_path)
    contract = _contract(loo=True)
    result = PyDESeq2LeaveOneDonorOutCapability().execute(contract, registry)
    _assert_admitted(contract, result, registry)
    fits = registry.load_payload(contract.parameters["fit_receipts_uri"])["fits"]
    assert len(fits) == 5
    assert len({fit["fit_id"] for fit in fits}) == 5
    assert [fit["left_out_donor"] for fit in fits if any("specified donor d2" in item["message"] for item in fit["warnings"])] == ["d2"]
    failed = next(fit for fit in fits if fit["left_out_donor"] == "d3")
    assert failed["status"] == "failed"
    assert failed["error"] == {"type": "RuntimeError", "message": "specified donor d3 diagnostic failure"}
    assert failed["design_matrix"] is not None
    assert result.metrics["n_failed_fits"] == 1
    assert result.metrics["scientific_robustness_claim_supported"] is False
    for fit in fits:
        assert len(fit["sample_ids"]) == (8 if fit["left_out_donor"] is None else 6)
        assert all(not value.startswith(str(fit["left_out_donor"]) + "__") for value in fit["sample_ids"]) if fit["left_out_donor"] else True


@pytest.mark.parametrize("tamper", ["counts", "design", "fold", "coverage", "options"])
def test_independent_audit_rejects_forged_rehashed_receipts(tmp_path, monkeypatch, tamper):
    registry = _registry(tmp_path)
    contract = _contract()
    result = PyDESeq2PseudobulkCapability().execute(contract, registry)
    _assert_admitted(contract, result, registry)
    receipt_uri = contract.parameters["fit_receipts_uri"]
    forged = deepcopy(registry.load_payload(receipt_uri))
    fit = forged["fits"][0]
    if tamper == "counts":
        fit["counts"]["frame"]["data"][0][0] += 1
        fit["counts"]["sha256"] = canonical_sha256(fit["counts"]["frame"])
        fit["counts_sha256"] = fit["counts"]["sha256"]
    elif tamper == "design":
        fit["design_matrix"]["frame"]["data"][0][0] += 1
        fit["design_matrix"]["sha256"] = canonical_sha256(fit["design_matrix"]["frame"])
    elif tamper == "fold":
        fit["left_out_donor"] = "d2"
    elif tamper == "coverage":
        forged["fits"] = []
        result.metrics["n_fit_receipts"] = 0
    elif tamper == "options":
        fit["effective_options"]["alpha"] = .9
    identity = {key: fit[key] for key in ("task_id", "method", "left_out_donor", "counts_sha256", "metadata_sha256", "design_formula", "contrast")}
    fit["fit_id"] = canonical_sha256(identity)
    fit["receipt_sha256"] = canonical_sha256({key: value for key, value in fit.items() if key != "receipt_sha256"})
    original_get = registry.get

    def get(uri):
        metadata, payload = original_get(uri)
        return (metadata, forged) if uri == receipt_uri else (metadata, payload)

    monkeypatch.setattr(registry, "get", get)
    report = ScientificAuditor().audit_task(contract, result, registry)
    check = next(check for check in report.checks if check.check_name == "advanced_statistics_fit_receipts")
    assert not check.passed and not report.overall_passed
    assert "hash differs" not in check.message  # Rehashing cannot disguise input/design/coverage inconsistency.


def test_underpowered_loo_records_no_fits_and_remains_an_explicit_skip(tmp_path):
    registry = _registry(tmp_path, n_pairs=2)
    contract = _contract(loo=True)
    result = PyDESeq2LeaveOneDonorOutCapability().execute(contract, registry)
    _assert_admitted(contract, result, registry)
    payload = registry.load_payload(contract.parameters["fit_receipts_uri"])
    assert payload["fits"] == []
    assert payload["fit_status"] == "not_attempted"
    assert payload["not_attempted_reason"]


def test_warning_as_error_policy_is_preserved():
    frame = pd.DataFrame([[1]], index=["sample"], columns=["gene"])
    receipt = FitReceipt(task_id="warning_policy", method="test", counts=frame,
                         metadata=pd.DataFrame({"condition": ["A"]}, index=frame.index),
                         design_formula="~condition", contrast={}, params={}, left_out_donor=None)
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        with pytest.raises(UserWarning, match="caller policy"):
            with receipt:
                warnings.warn("caller policy", UserWarning)
    assert receipt.record["status"] == "failed"
    assert receipt.record["error"]["type"] == "UserWarning"
