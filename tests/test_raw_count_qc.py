"""Raw H5AD inputs must not silently bypass mitochondrial QC."""

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.computational import ComputationalValidator
from eacbp.capabilities.qc import QCCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def _run(tmp_path, data, parameters=None):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://raw_qc/input/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "raw_qc", "input", "test")
    contract = TaskContract(
        task_id="qc", capability="qc", input_artifacts=[uri],
        parameters={"min_genes": 1, "max_mito_pct": 20, **(parameters or {})},
    )
    result = QCCapability().execute(contract, registry)
    return contract, result, registry


def _execute(tmp_path, data, parameters=None):
    _, result, registry = _run(tmp_path, data, parameters)
    return result, registry.load_scdata(result.output_artifacts[0])


@pytest.mark.parametrize("as_sparse", [False, True])
def test_mito_percent_computed_from_raw_counts(tmp_path, as_sparse):
    counts = np.array([[8, 2], [1, 9], [0, 0]])
    matrix = sparse.csr_matrix(counts) if as_sparse else counts
    data = SCData(matrix, pd.DataFrame(index=["high", "low", "zero"]),
                  pd.DataFrame(index=["MT-CO1", "GENE"]))
    result, output = _execute(tmp_path, data)
    assert list(output.obs.index) == ["low"]
    assert output.obs.percent_mito.tolist() == [10.0]
    assert result.metrics["mitochondrial_filtering_applied"] is True
    assert result.metrics["n_mitochondrial_genes"] == 1
    assert "percent_mito" not in data.obs


def test_counts_layer_drives_qc_with_normalized_x(tmp_path):
    data = SCData(np.ones((2, 2)), pd.DataFrame(index=["high", "low"]),
                  pd.DataFrame({"mt": [True, False]}, index=["ENSG1", "ENSG2"]),
                  layers={"counts": sparse.csr_matrix([[8, 2], [1, 9]])})
    result, output = _execute(tmp_path, data)
    assert list(output.obs.index) == ["low"]
    assert result.metrics["counts_source"] == "layers.counts"
    assert result.metrics["mitochondrial_metric_source"] == "var.mt"


def test_unknown_mitochondrial_features_not_reported_as_filtered(tmp_path):
    data = SCData(np.ones((2, 2)), pd.DataFrame(index=["a", "b"]),
                  pd.DataFrame(index=["ENSG1", "ENSG2"]))
    result, output = _execute(tmp_path, data)
    assert output.n_obs == 2
    assert "percent_mito" not in output.obs
    assert result.metrics["mitochondrial_filtering_applied"] is False
    assert "mitochondrial_filtering" not in result.executed_operations
    table = output.uns["qc"]["cell_disposition"]
    assert np.isnan(table["mito_pct"]).all()
    assert not table["mito_pct_available"].any()
    assert not table["high_mito_flag"].any()


def test_invalid_precomputed_mito_percent_rejected(tmp_path):
    data = SCData(np.ones((2, 2)), pd.DataFrame({"percent_mito": [np.nan, 5]}),
                  pd.DataFrame(index=["MT-CO1", "GENE"]))
    with pytest.raises(ValueError, match="finite percentages"):
        _execute(tmp_path, data)


def _disposition_input(as_sparse=False):
    values = np.array([[1, 9, 2], [0, 1, 0], [8, 1, 1], [4, 0, 0], [0, 0, 0]])
    return SCData(sparse.csr_matrix(values) if as_sparse else values,
                  pd.DataFrame(index=["retained", "low_rna", "high_mito", "both", "zero"]),
                  pd.DataFrame(index=["MT-GENE", "GENE1", "GENE2"]))


@pytest.mark.parametrize("as_sparse", [False, True])
def test_all_input_cell_dispositions_persist_and_independent_audit_passes(tmp_path, as_sparse):
    contract, result, registry = _run(tmp_path, _disposition_input(as_sparse), {"min_genes": 2})
    output = registry.load_scdata(result.output_artifacts[0])
    qc = output.uns["qc"]
    table = qc["cell_disposition"]
    assert list(table["ids"]) == ["retained", "low_rna", "high_mito", "both", "zero"]
    assert list(table["reason"]) == ["retained", "low_gene_count", "high_mito_percentage",
                                    "low_gene_count_and_high_mito_percentage", "low_gene_count"]
    assert table["retained"].tolist() == [True, False, False, False, False]
    assert table["n_genes"].tolist() == [3, 1, 3, 1, 0]
    assert table["total_counts"].tolist() == [12, 1, 10, 4, 0]
    assert table["mito_pct_available"].tolist() == [True, True, True, True, False]
    assert np.isnan(table["mito_pct"][-1])
    assert list(output.obs.index) == ["retained"]
    assert result.metrics["cell_disposition_counts"] == {
        "retained": 1, "low_gene_count": 2, "high_mito_percentage": 1,
        "low_gene_count_and_high_mito_percentage": 1,
    }
    assert "does not independently establish cellular damage" in qc["interpretation_warning"]
    assert "not a damage label" in qc["interpretation_warning"]
    assert qc["input_sha256"] == registry.get_metadata(contract.input_artifacts[0]).sha256_hash
    report = ComputationalValidator().audit(contract, result, registry)
    check = next(check for check in report.checks if check.check_name == "qc_cell_disposition_integrity")
    assert report.overall_passed
    assert check.metrics["all_input_cells_assessed"] is True
    assert check.metrics["initial_cells"] == 5


def test_qc_observability_preserves_threshold_boundaries(tmp_path):
    data = SCData(np.array([[2, 8, 0], [3, 7, 0], [0, 1, 0]]),
                  pd.DataFrame(index=["exact", "above", "below"]),
                  pd.DataFrame(index=["MT-GENE", "GENE1", "GENE2"]))
    contract, result, registry = _run(tmp_path, data, {"min_genes": 2, "max_mito_pct": 20})
    output = registry.load_scdata(result.output_artifacts[0])
    assert list(output.obs.index) == ["exact"]
    assert output.uns["qc"]["cell_disposition"]["low_gene_flag"].tolist() == [False, False, True]
    assert output.uns["qc"]["cell_disposition"]["high_mito_flag"].tolist() == [False, True, False]
    assert ComputationalValidator().audit(contract, result, registry).overall_passed


def test_zero_count_cell_mask_unchanged_when_zero_gene_threshold_requested(tmp_path):
    data = SCData(np.zeros((1, 2)), pd.DataFrame(index=["zero"]),
                  pd.DataFrame(index=["MT-GENE", "GENE"]))
    contract, result, registry = _run(tmp_path, data, {"min_genes": 0})
    output = registry.load_scdata(result.output_artifacts[0])
    assert list(output.obs.index) == ["zero"]  # Historical requested mask.
    assert output.obs.percent_mito.tolist() == [0.0]
    table = output.uns["qc"]["cell_disposition"]
    assert table["retained"].tolist() == [True]
    assert table["mito_pct_available"].tolist() == [False]
    assert np.isnan(table["mito_pct"][0])  # No scientific zero fraction invented.
    assert ComputationalValidator().audit(contract, result, registry).overall_passed


@pytest.mark.parametrize("metric_source", ["counts_layer", "precomputed", "unavailable"])
def test_independent_qc_audit_accepts_supported_metric_sources(tmp_path, metric_source):
    if metric_source == "counts_layer":
        data = SCData(np.ones((2, 2)), pd.DataFrame(index=["a", "b"]),
                      pd.DataFrame({"mt": [True, False]}, index=["ENSG1", "ENSG2"]),
                      layers={"counts": sparse.csr_matrix([[8, 2], [1, 9]])})
    elif metric_source == "precomputed":
        data = SCData(np.ones((2, 2)), pd.DataFrame({"percent_mito": [25.0, 10.0]}, index=["a", "b"]),
                      pd.DataFrame(index=["GENE1", "GENE2"]))
    else:
        data = SCData(np.ones((2, 2)), pd.DataFrame(index=["a", "b"]),
                      pd.DataFrame(index=["GENE1", "GENE2"]))
    contract, result, registry = _run(tmp_path, data)
    report = ComputationalValidator().audit(contract, result, registry)
    assert report.overall_passed, report


@pytest.mark.parametrize("tampering", ["missing_cell", "damage_label", "counts", "retained_flag",
                                       "retained_id", "input_hash", "missing_table", "fake_zero_mito"])
def test_independent_qc_audit_rejects_tampered_dispositions(tmp_path, tampering):
    contract, result, registry = _run(tmp_path, _disposition_input(), {"min_genes": 2})
    original_meta = registry.get_metadata(result.output_artifacts[0])
    changed = registry.load_scdata(result.output_artifacts[0])
    qc = changed.uns["qc"]
    table = qc["cell_disposition"]
    if tampering == "missing_cell":
        qc["cell_disposition"] = {key: np.asarray(value)[:-1] for key, value in table.items()}
    elif tampering == "damage_label":
        table["reason"][1] = "damaged"
    elif tampering == "counts":
        qc["cell_disposition_counts"]["low_gene_count"] = 0
    elif tampering == "retained_flag":
        table["retained"][1] = True
    elif tampering == "retained_id":
        changed.obs.index = pd.Index(["low_rna"])
    elif tampering == "input_hash":
        qc["input_sha256"] = "sha256:" + "0" * 64
    elif tampering == "missing_table":
        del qc["cell_disposition"]
    elif tampering == "fake_zero_mito":
        table["mito_pct"][-1] = 0.0
    # Register internally coherent bytes so failure tests independent scientific
    # validation, rather than merely the registry's file-tampering guard.
    changed_uri = "adata://raw_qc/changed_qc/v1"
    registry.register(changed_uri, changed, ArtifactType.ANNDATA, "raw_qc", "qc", "filter_cells_qc",
                      parent_uris=contract.input_artifacts, parameters=original_meta.parameters,
                      summary_metrics=original_meta.summary_metrics)
    changed_result = result.model_copy(update={"output_artifacts": [changed_uri]})
    report = ComputationalValidator().audit(contract, changed_result, registry)
    check = next(check for check in report.checks if check.check_name == "qc_cell_disposition_integrity")
    assert not report.overall_passed
    assert not check.passed
    assert check.metrics["all_input_cells_assessed"] is False


def test_legacy_qc_artifact_disposition_missing_is_explicitly_unassessed(tmp_path):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://raw_qc/legacy_qc/v1"
    data = SCData(np.ones((1, 1)), pd.DataFrame(index=["a"]), pd.DataFrame(index=["GENE"]))
    registry.register(uri, data, ArtifactType.ANNDATA, "raw_qc", "qc", "filter_cells_qc")
    from eacbp.schemas.task import TaskResult, TaskStatus
    result = TaskResult(task_id="qc", capability="qc", status=TaskStatus.SUCCESS,
                        method_used="sc_qc_v1", output_artifacts=[uri])
    report = ComputationalValidator().audit(TaskContract(task_id="qc", capability="qc"), result, registry)
    check = next(check for check in report.checks if check.check_name == "qc_cell_disposition_integrity")
    assert report.overall_passed
    assert check.metrics == {"cell_disposition_available": False, "all_input_cells_assessed": False}
