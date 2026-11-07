"""Independent reconstruction must not inflate donors with technical libraries."""

from types import SimpleNamespace
import hashlib

import numpy as np
import pandas as pd
import pytest

from eacbp.auditor.advanced_statistics import AdvancedStatisticsValidator, _reconstruct_input
from eacbp.auditor.advanced_communication import LianaValidator
from eacbp.auditor.requirements import run_payload_requirements
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


METHODS = (
    "pydeseq2_pseudobulk_v1",
    "decoupler_ulm_v2",
    "pydeseq2_leave_one_donor_out_v1",
)


def _data(*, biological_column=None):
    observations = []
    for donor in ("m1", "m2", "m3"):
        for condition in ("control", "treated"):
            for library in range(4):
                identity = f"{donor}_{condition}_library{library}"
                row = {"condition": condition, "sample_id": identity, "library_id": identity}
                if biological_column:
                    row[biological_column] = donor
                observations.append(row)
    counts = np.tile(np.array([2, 3, 5], dtype=np.int64), (len(observations), 1))
    return SCData(
        X=counts.astype(float), obs=pd.DataFrame(observations),
        var=pd.DataFrame(index=["g1", "g2", "g3"]), layers={"counts": counts},
    )


def _reconstruct(data, *, method=METHODS[0], **extra):
    return _reconstruct_input(
        data, {"condition_a": "treated", "condition_b": "control", **extra}, {}, method,
    )


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("technical_column", ["sample_id", "sample", "library_id"])
def test_technical_identities_never_silently_supply_donors(method, technical_column):
    data = _data()
    data.obs = data.obs[["condition", "sample_id"]].rename(columns={"sample_id": technical_column})
    with pytest.raises(ValueError, match="Biological donor identity is unavailable"):
        _reconstruct(data, method=method)


@pytest.mark.parametrize("column", ["mouse_id", "mouse", "subject_id", "subject", "patient_id", "patient"])
def test_biological_aliases_keep_four_libraries_as_one_donor_condition_unit(column):
    data = _data(biological_column=column)
    details = _reconstruct(data, paired=True)
    assert details.donor_col == column
    assert details.donors_a == details.donors_b == ["m1", "m2", "m3"]
    assert len(details.metadata) == 6
    assert details.design_rank == 4 and details.residual_df == 2


@pytest.mark.parametrize("method", METHODS)
def test_explicit_library_mapping_cannot_split_observed_biological_individual(method):
    with pytest.raises(ValueError, match="not independent biological replicates"):
        _reconstruct(_data(biological_column="mouse_id"), method=method, donor_col="sample_id")


def test_conflicting_biological_aliases_require_explicit_resolution():
    data = _data(biological_column="mouse_id")
    data.obs["donor"] = "different_" + data.obs["mouse_id"]
    with pytest.raises(ValueError, match="Ambiguous donor metadata"):
        _reconstruct(data, paired=True)
    assert _reconstruct(data, paired=True, donor_col="mouse_id").donors_a == ["m1", "m2", "m3"]


@pytest.mark.parametrize("missing", ["", "  ", "null", "NaN", None])
def test_incomplete_donor_identity_is_rejected(missing):
    data = _data(biological_column="donor")
    data.obs.loc[0, "donor"] = missing
    with pytest.raises(ValueError, match="missing values"):
        _reconstruct(data, paired=True)


def test_forged_successful_statistics_cannot_bypass_input_donor_audit():
    data = _data()
    input_uri, output_uri = "adata://units/input/v1", "table://units/deg/v1"
    contract = TaskContract(
        task_id="unit_audit", capability="deg", method=METHODS[0],
        input_artifacts=[input_uri], expected_outputs=[output_uri],
        parameters={"condition_a": "treated", "condition_b": "control"},
    )
    result = TaskResult(
        task_id=contract.task_id, capability=contract.capability, method_used=METHODS[0],
        status=TaskStatus.SUCCESS, input_artifacts=[input_uri], output_artifacts=[output_uri],
        metrics={"n_donors_condition_a": 12, "n_donors_condition_b": 12, "biological_replication_sufficient": True},
    )
    payloads = {input_uri: data, output_uri: pd.DataFrame({"gene": ["g1", "g2", "g3"]})}
    registry = SimpleNamespace(get=lambda uri: (None, payloads[uri]))
    report = AdvancedStatisticsValidator().audit(contract, result, registry)
    check = next(check for check in report.checks if check.check_name == "advanced_statistics_input_counts_metadata")
    assert not check.passed and not report.overall_passed
    assert "Biological donor identity is unavailable" in check.message


def _replicate_check(data, recorded_minimum, *, recorded_sufficient=False, **parameters):
    input_uri, output_uri = "adata://units/input/v1", "table://units/audit/v1"
    contract = TaskContract(
        task_id="replication", capability="dataset_audit", input_artifacts=[input_uri],
        expected_outputs=[output_uri], validation_requirements=["replicate_check"], parameters=parameters,
    )
    result = TaskResult(
        task_id=contract.task_id, capability=contract.capability, method_used="sc_audit_v1",
        status=TaskStatus.SUCCESS, input_artifacts=[input_uri], output_artifacts=[output_uri],
    )
    audit_table = pd.DataFrame([{"min_replicates_per_condition": recorded_minimum,
                                 "biological_replication_sufficient": recorded_sufficient}])
    payloads = {
        input_uri: (SimpleNamespace(type=ArtifactType.ANNDATA), data),
        output_uri: (SimpleNamespace(type=ArtifactType.TABLE), audit_table),
    }
    checks = run_payload_requirements(contract, result, SimpleNamespace(get=lambda uri: payloads[uri]), [])
    return next(check for check in checks if check.check_name == "replicate_check")


@pytest.mark.parametrize("technical_column", ["sample_id", "sample", "library_id"])
def test_descriptive_dataset_audit_without_donors_remains_admissible(technical_column):
    data = _data()
    data.obs = data.obs[["condition", "sample_id"]].rename(columns={"sample_id": technical_column})
    check = _replicate_check(data, 0)
    assert check.passed
    assert check.metrics["min_replicates"] == 0
    assert check.metrics["replication_assessed"] is False
    assert check.metrics["donor_col"] is None
    assert check.metrics["donor_independence_established"] is False


@pytest.mark.parametrize("reported_count,reported_sufficient", [(12, True), (1, False), (0, True)])
def test_dataset_audit_rejects_technical_counts_or_false_replication_claim(reported_count, reported_sufficient):
    check = _replicate_check(_data(), reported_count, recorded_sufficient=reported_sufficient)
    assert not check.passed
    assert check.metrics["min_replicates"] == 0


@pytest.mark.parametrize("column", ["mouse", "subject", "patient_id"])
def test_dataset_audit_recounts_biological_alias_without_inflating_libraries(column):
    check = _replicate_check(_data(biological_column=column), 3, recorded_sufficient=True)
    assert check.passed and check.metrics["replication_assessed"]
    assert check.metrics["donor_col"] == column
    assert check.metrics["min_replicates"] == 3
    assert check.metrics["donor_independence_established"] is False


def test_dataset_audit_rejects_explicit_specimen_mapping_that_splits_mouse():
    check = _replicate_check(
        _data(biological_column="mouse_id"), 12, recorded_sufficient=True, donor_col="sample_id",
    )
    assert not check.passed
    assert "not independent biological replicates" in check.message


@pytest.mark.parametrize("donor_column,expected_pass", [("sample_id", False), ("mouse_id", True)])
def test_liana_auditor_rejects_library_donors_but_keeps_real_donor_ranks(tmp_path, donor_column, expected_pass):
    data = _data(biological_column="mouse_id")
    data.obs["cell_type"] = "T"
    resource_path = tmp_path / "lr.csv"
    pd.DataFrame({"ligand": ["g1"], "receptor": ["g2"]}).to_csv(resource_path, index=False)
    digest = hashlib.sha256(resource_path.read_bytes()).hexdigest()
    parameters = {
        "species": "human", "donor_col": donor_column,
        "condition_col": "condition", "cell_type_col": "cell_type", "spatial_key": None,
        "resource": {
            "resource_path": str(resource_path.resolve()), "resource_sha256": digest,
            "resource_version": "fixture-v1", "resource_source": "test fixture", "species": "human",
            "resource_species_values": [], "resource_species_columns": [],
        },
        "rank_aggregate": {"return_all_lrs": True},
    }
    rows = []
    for (donor, condition), size in data.obs.groupby([donor_column, "condition"]).size().items():
        rows.append({
            "donor_id": donor, "condition": condition, "donor_condition": f"{donor}::{condition}",
            "n_cells_group": size, "source": "T", "target": "T", "ligand": "g1", "receptor": "g2",
            "magnitude_rank": .2, "specificity_rank": .3, "species": "human",
            "resource_version": "fixture-v1", "resource_source": "test fixture", "resource_sha256": digest,
            "rank_is_fdr": False, "spatial_used": False,
        })
    input_uri, output_uri = "adata://units/input/v1", "table://units/ranks/v1"
    payloads = {
        input_uri: (SimpleNamespace(type=ArtifactType.ANNDATA), data),
        output_uri: (SimpleNamespace(type=ArtifactType.TABLE, parameters=parameters, summary_metrics={}), pd.DataFrame(rows)),
    }
    contract = TaskContract(
        task_id="liana_units", capability="liana_communication", input_artifacts=[input_uri],
        expected_outputs=[output_uri], parameters={"species": "human", "donor_col": donor_column},
    )
    result = TaskResult(
        task_id=contract.task_id, capability=contract.capability, method_used="liana_rank_aggregate_v1",
        status=TaskStatus.SUCCESS, input_artifacts=[input_uri], output_artifacts=[output_uri],
    )
    report = LianaValidator().audit(contract, result, SimpleNamespace(get=lambda uri: payloads[uri]))
    assert report.overall_passed is expected_pass, [check.message for check in report.checks]
    if not expected_pass:
        assert "not independent biological replicates" in report.checks[0].message
