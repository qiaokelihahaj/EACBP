"""Biological replication must follow individuals, never the library count."""

import json

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.computational import ComputationalValidator
from eacbp.capabilities.advanced_statistics import AdvancedStatisticsInputError, _prepare_pseudobulk
from eacbp.capabilities.advanced_communication import LianaCommunicationCapability, LianaCommunicationInputError
from eacbp.capabilities.biological_units import resolve_donor_role
from eacbp.capabilities.deg import DifferentialAbundanceCapability, DifferentialExpressionCapability
from eacbp.capabilities.qc import DatasetAuditCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import ScientificResultStatus, TaskContract


def _data(*, n_mice=1, paired=False, biological_col="mouse_id"):
    rows = []
    counts = []
    for condition_idx, condition in enumerate(("control", "treated")):
        for mouse_idx in range(n_mice):
            mouse = f"m{mouse_idx}" if paired else f"{condition}_m{mouse_idx}"
            for library_idx in range(4):
                for cell_idx in range(2):
                    rows.append({
                        "condition": condition, biological_col: mouse,
                        "sample_id": f"{mouse}_{condition}_s{library_idx}",
                        "library_id": f"{mouse}_{condition}_l{library_idx}",
                        "cell_type": "state_a" if cell_idx == 0 else "state_b",
                    })
                    counts.append([3 + mouse_idx + library_idx, 6 + condition_idx * 5 + cell_idx])
    counts = np.asarray(counts, dtype=np.int64)
    return SCData(counts.copy(), pd.DataFrame(rows), pd.DataFrame({"gene_name": ["g1", "g2"]}), layers={"counts": counts})


def _registry(tmp_path, data):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://biological_units/input/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "biological_units", "input", "input")
    return registry, uri


def _contract(uri, capability, **params):
    return TaskContract(task_id="biological_unit_test", capability=capability, input_artifacts=[uri],
                        parameters={"condition_a": "treated", "condition_b": "control", **params})


def test_four_libraries_from_one_mouse_do_not_satisfy_replication(tmp_path):
    registry, uri = _registry(tmp_path, _data())
    audited = DatasetAuditCapability().execute(_contract(uri, "dataset_audit"), registry)
    assert audited.metrics["donor_replicates"] == {"control": 1, "treated": 1}
    assert audited.metrics["biological_replication_sufficient"] is False
    assert audited.metrics["donor_col"] == "mouse_id"
    assert audited.metrics["donor_role_source"] == "observed_biological_alias"
    with pytest.raises(AdvancedStatisticsInputError, match="at least 2 donors"):
        _prepare_pseudobulk(_data(), {})
    with pytest.raises(ValueError, match="not estimable.*two biological donors"):
        DifferentialAbundanceCapability().execute(_contract(uri, "differential_abundance"), registry)


@pytest.mark.parametrize("identifier", ["sample", "sample_id", "library", "library_id", "replicate"])
def test_technical_identifiers_never_automatically_become_donors(tmp_path, identifier):
    data = _data(n_mice=2)
    data.obs = data.obs[["condition", "cell_type"]].assign(**{identifier: [f"technical_{i}" for i in range(data.n_obs)]})
    registry, uri = _registry(tmp_path, data)
    audited = DatasetAuditCapability().execute(_contract(uri, "dataset_audit"), registry)
    assert audited.metrics["donor_col"] is None
    assert audited.metrics["biological_replication_sufficient"] is False
    result = DifferentialExpressionCapability().execute(_contract(uri, "deg"), registry)
    table = registry.load_payload(result.output_artifacts[0])
    assert result.scientific_result.status == ScientificResultStatus.NOT_ESTIMABLE
    assert table["log2_fold_change"].notna().all()
    assert table["p_value"].isna().all() and table["fdr_q_value"].isna().all()
    assert not any("t_test" in operation for operation in result.executed_operations)
    with pytest.raises(AdvancedStatisticsInputError, match="not estimable.*sample/library"):
        _prepare_pseudobulk(data, {})


@pytest.mark.parametrize("biological_col", ["donor", "donor_id", "subject", "subject_id", "patient", "patient_id", "mouse", "mouse_id"])
def test_individual_aliases_preserve_real_independent_replication(tmp_path, biological_col):
    data = _data(n_mice=2, biological_col=biological_col)
    registry, uri = _registry(tmp_path, data)
    result = DifferentialExpressionCapability().execute(_contract(uri, "deg"), registry)
    table = registry.load_payload(result.output_artifacts[0])
    assert result.metrics["is_pseudobulk"] is True
    assert result.metrics["donor_col"] == biological_col
    assert table["n_donors_condition_a"].eq(2).all()
    assert table["n_donors_condition_b"].eq(2).all()
    assert table["p_value"].notna().all()
    assert result.metrics["donor_independence_established"] is False


def test_paired_individuals_aggregate_all_libraries_into_four_pseudobulks():
    data = _data(n_mice=2, paired=True)
    prepared = _prepare_pseudobulk(data, {"paired": True})
    assert len(prepared.metadata) == 4
    assert prepared.donor_ids_a == prepared.donor_ids_b == ["m0", "m1"]
    for sample, metadata in prepared.metadata.iterrows():
        selected = ((data.obs["mouse_id"] == metadata["donor"]) & (data.obs["condition"] == metadata["condition"]))
        np.testing.assert_array_equal(prepared.counts.loc[sample].to_numpy(), data.layers["counts"][selected].sum(axis=0))


@pytest.mark.parametrize("explicit", ["sample_id", "library_id"])
def test_explicit_technical_mapping_cannot_split_a_known_mouse(tmp_path, explicit):
    data = _data()
    registry, uri = _registry(tmp_path, data)
    with pytest.raises(ValueError, match="not independent biological replicates"):
        DatasetAuditCapability().execute(_contract(uri, "dataset_audit", donor_col=explicit), registry)
    with pytest.raises(ValueError, match="not independent biological replicates"):
        DifferentialExpressionCapability().execute(_contract(uri, "deg", donor_col=explicit), registry)
    with pytest.raises(AdvancedStatisticsInputError, match="not independent biological replicates"):
        _prepare_pseudobulk(data, {"donor_col": explicit})


def test_specimens_across_conditions_cannot_hide_the_same_paired_individual():
    data = _data(n_mice=2, paired=True)
    data.obs["sample_id"] = data.obs["mouse_id"] + "_" + data.obs["condition"]
    with pytest.raises(AdvancedStatisticsInputError, match="not independent biological replicates"):
        _prepare_pseudobulk(data, {"donor_col": "sample_id"})
    assert _prepare_pseudobulk(data, {"paired": True}).paired is True


@pytest.mark.parametrize("alias", ["donor", "donor_id"])
def test_explicit_biological_alias_cannot_disguise_library_ids_as_donors(tmp_path, alias):
    data = _data()
    data.obs[alias] = data.obs["library_id"]
    registry, uri = _registry(tmp_path, data)
    with pytest.raises(ValueError, match="splits one biological individual.*mouse_id"):
        DatasetAuditCapability().execute(_contract(uri, "dataset_audit", donor_col=alias), registry)
    with pytest.raises(AdvancedStatisticsInputError, match="not independent biological replicates"):
        _prepare_pseudobulk(data, {"donor_col": alias})


def test_explicit_one_to_one_individual_aliases_can_use_different_names():
    data = _data(n_mice=2, paired=True)
    data.obs["donor"] = "author_donor_" + data.obs["mouse_id"]
    prepared = _prepare_pseudobulk(data, {"donor_col": "donor", "paired": True})
    assert prepared.donor_ids_a == prepared.donor_ids_b == ["author_donor_m0", "author_donor_m1"]
    assert len(prepared.metadata) == 4
    assert prepared.donor_role.source == "user_explicit"


def test_absent_donor_role_survives_csv_reload_without_numeric_nan(tmp_path):
    data = _data()
    data.obs = data.obs.drop(columns=["mouse_id"])
    registry, uri = _registry(tmp_path, data)
    contract = _contract(uri, "dataset_audit")
    result = DatasetAuditCapability().execute(contract, registry)
    reopened = ArtifactRegistry(str(tmp_path / "artifacts"))
    table = reopened.load_payload(result.output_artifacts[0])
    assert table["donor_col"].eq("unavailable").all()
    assert result.metrics["donor_col"] is None
    audited = ComputationalValidator().audit(contract, result, reopened)
    assert all(check.passed for check in audited.checks), audited.checks


def test_missing_donor_audit_is_admitted_in_native_lifecycle(tmp_path, monkeypatch):
    import anndata as ad
    from eacbp.application import autonomous
    from eacbp.orchestrator import loop

    original = loop.build_study_tasks
    monkeypatch.setattr(loop, "build_study_tasks", lambda *args, **kwargs:
                        [task for task in original(*args, **kwargs) if task.capability == "dataset_audit"])
    data = _data()
    data.obs = data.obs.drop(columns=["mouse_id"])
    source = tmp_path / "source.h5ad"
    ad.AnnData(np.column_stack([data.X, np.ones(data.n_obs, dtype=np.int64)]),
               obs=data.obs, var=pd.DataFrame({"gene_name": ["g1", "g2", "g3"]})).write_h5ad(source)
    run_dir = tmp_path / "run"
    result = autonomous.analyze_study(question="Discover cell populations", data=source,
                                     species="human", tissue="blood", run_dir=run_dir)
    assert result["status"] == "success", result
    assert result["analysis_request_status"] == "partial"
    saved = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))
    assert saved["status"] == "success" and saved["analysis_request_status"] == "partial"


def _communication_task(tmp_path, data, **params):
    resource = tmp_path / "lr.csv"
    pd.DataFrame({"ligand": ["g1"], "receptor": ["g2"]}).to_csv(resource, index=False)
    registry, uri = _registry(tmp_path, data)
    task = _contract(uri, "liana_communication", species="human", lr_resource_version="test-v1",
                     lr_resource_source="local test fixture", lr_resource_path=str(resource),
                     min_cells=2, **params)
    return registry, task


@pytest.mark.parametrize("case", ["sample_only", "explicit_split", "blank"])
def test_communication_rejects_invalid_donor_roles_before_loading_liana(tmp_path, monkeypatch, case):
    data = _data()
    params = {}
    if case == "sample_only":
        data.obs = data.obs.drop(columns=["mouse_id"])
    elif case == "explicit_split":
        params["donor_col"] = "library_id"
    else:
        data.obs.loc[0, "mouse_id"] = " "
    registry, task = _communication_task(tmp_path, data, **params)
    def reject_import():
        raise AssertionError("Invalid biological roles must fail before loading the statistical method")
    monkeypatch.setattr("eacbp.capabilities.advanced_communication._load_liana", reject_import)
    with pytest.raises(LianaCommunicationInputError):
        LianaCommunicationCapability().execute(task, registry)


def test_communication_strata_and_explicit_role_receipt_do_not_count_libraries(tmp_path, monkeypatch):
    from types import SimpleNamespace
    data = _data()
    data.obs["specimen_identity"] = data.obs.pop("mouse_id")
    registry, task = _communication_task(tmp_path, data, donor_col="specimen_identity")
    def ranks(*args, **kwargs):
        return pd.DataFrame({"source": ["state_a"], "target": ["state_b"],
                             "ligand": ["g1"], "receptor": ["g2"],
                             "magnitude_rank": [0.2], "specificity_rank": [0.3]})
    # This test covers role routing and receipts, not LIANA numerical inference.
    monkeypatch.setattr("eacbp.capabilities.advanced_communication._load_liana",
                        lambda: SimpleNamespace(mt=SimpleNamespace(rank_aggregate=ranks)))
    result = LianaCommunicationCapability().execute(task, registry)
    assert result.metrics["donors_by_condition"] == {"control": ["control_m0"], "treated": ["treated_m0"]}
    assert result.metrics["comparison_status"] == "insufficient_donors"
    assert result.metrics["donor_role_source"] == "user_explicit"
    assert result.metrics["donor_independence_established"] is False
    assert registry.get_metadata(result.output_artifacts[0]).parameters["donor_role_source"] == "user_explicit"


def test_explicit_user_mapping_retains_role_source_and_independence_limit(tmp_path):
    data = _data(n_mice=2)
    data.obs["sample_id"] = data.obs.pop("mouse_id")
    registry, uri = _registry(tmp_path, data)
    result = DifferentialExpressionCapability().execute(_contract(uri, "deg", donor_col="sample_id"), registry)
    table = registry.load_payload(result.output_artifacts[0])
    assert result.metrics["donor_role_source"] == "user_explicit"
    assert "sampling design" in result.metrics["donor_role_limitation"]
    assert result.metrics["donor_independence_established"] is False
    assert table["donor_role_source"].eq("user_explicit").all()
    stored = registry.get_metadata(result.output_artifacts[0])
    assert stored.parameters["donor_role_source"] == "user_explicit"
    prepared = _prepare_pseudobulk(data, {"donor_col": "sample_id"})
    assert prepared.donor_role.source == "user_explicit"


def test_missing_explicit_column_and_ambiguous_individual_roles_are_rejected():
    data = _data(n_mice=2)
    with pytest.raises(ValueError, match="not present"):
        resolve_donor_role(data.obs, {"donor_col": "absent"})
    data.obs["donor"] = [f"distinct_{i}" for i in range(data.n_obs)]
    with pytest.raises(ValueError, match="Ambiguous donor metadata"):
        resolve_donor_role(data.obs, {})


@pytest.mark.parametrize("missing", [None, "", "  ", "nan"])
def test_incomplete_donor_identity_is_descriptive_and_never_a_partial_bulk_fit(tmp_path, missing):
    data = _data(n_mice=2)
    data.obs.loc[0, "mouse_id"] = missing
    registry, uri = _registry(tmp_path, data)
    result = DifferentialExpressionCapability().execute(_contract(uri, "deg"), registry)
    assert result.metrics["is_pseudobulk"] is False
    assert result.scientific_result.status == ScientificResultStatus.NOT_ESTIMABLE
    assert registry.load_payload(result.output_artifacts[0])["p_value"].isna().all()
    with pytest.raises(AdvancedStatisticsInputError, match="missing or empty"):
        _prepare_pseudobulk(data, {})
