"""Candidate counts remain exact, isolated and independently auditable."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from pydantic import ValidationError
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.qc_review import (
    CELL_COLUMNS, FEATURE_COLUMNS, INTERPRETATION_WARNING, QCReviewCapability,
    QCReviewParameters, REASON_NAMES, _plan_qc_review,
)
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskStatus


def _data(as_sparse=False):
    counts = np.array([[1, 9, 2], [0, 1, 0], [8, 1, 1], [4, 0, 0], [0, 0, 0]])
    if as_sparse:
        counts = sparse.csr_matrix(counts)
    return SCData(
        counts, pd.DataFrame({"cell_type": ["a", "b", "c", "a", "b"]},
                             index=["retained", "low_rna", "high_mito", "both", "zero"]),
        pd.DataFrame(index=["MT-GENE", "GENE1", "GENE2"]),
        obsm={"X_pca": np.ones((5, 2))}, uns={"old_reference": {"answer": "secret"}},
    )


def _run(tmp_path, data=None, parameters=None):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    in_uri = "adata://review_test/input/v1"
    registry.register(in_uri, data or _data(), ArtifactType.ANNDATA, "review_test", "input", "test")
    contract = TaskContract(
        task_id="review", capability="qc_review", method="raw_count_qc_review_v1",
        input_artifacts=[in_uri], parameters={"min_genes": 2, "max_mito_pct": 20.0, **(parameters or {})},
        expected_outputs=["json://review_test/qc_review/v1", "adata://review_test/qc_candidates/v1",
                          "table://review_test/qc_review_cells/v1", "table://review_test/qc_review_features/v1"],
    )
    result = QCReviewCapability().execute(contract, registry)
    return contract, result, registry


@pytest.mark.parametrize("as_sparse", [False, True])
def test_candidate_counts_full_cell_table_and_all_feature_groups(tmp_path, as_sparse):
    original = _data(as_sparse)
    contract, result, registry = _run(tmp_path, original)
    summary = registry.load_payload(result.output_artifacts[0])
    candidate = registry.load_scdata(result.output_artifacts[1])
    cells = registry.load_payload(result.output_artifacts[2])
    features = registry.load_payload(result.output_artifacts[3])
    assert result.status == TaskStatus.SUCCESS
    assert result.metrics == summary
    assert summary["candidate_cells"] == 4 and summary["retained_cells"] == 1
    assert summary["restoration_performed"] is False
    assert summary["damage_diagnosis"] is False and summary["reference_labels_used"] is False
    assert summary["interpretation_warning"] == INTERPRETATION_WARNING
    assert tuple(cells.columns) == CELL_COLUMNS and tuple(features.columns) == FEATURE_COLUMNS
    assert cells.cell_id.tolist() == ["retained", "low_rna", "high_mito", "both", "zero"]
    assert cells.reason.tolist() == ["retained", "low_gene_count", "high_mito_percentage",
                                   "low_gene_count_and_high_mito_percentage", "low_gene_count"]
    assert cells.retained.tolist() == [True, False, False, False, False]
    assert not cells.loc[4, "mito_pct_available"] and np.isnan(cells.loc[4, "mito_pct"])
    assert candidate.obs.index.tolist() == cells.loc[~cells.retained, "cell_id"].tolist()
    expected = original.X[1:]
    expected = expected.toarray() if sparse.issparse(expected) else expected
    np.testing.assert_array_equal(candidate.X.toarray(), expected)
    assert sparse.isspmatrix_csr(candidate.X)
    assert candidate.var.equals(original.var)
    assert candidate.obs.qc_review_candidate.all() and not candidate.obs.qc_passed.any()
    assert set(candidate.uns) == {"qc_review"} and candidate.uns["qc_review"] == summary
    assert set(candidate.layers) == {"counts"}
    assert not candidate.obsm and not candidate.obsp and candidate.raw is None
    assert not candidate.varm and not candidate.varp
    assert features.shape == (12, len(FEATURE_COLUMNS))
    assert features.groupby("reason").size().to_dict() == {reason: 3 for reason in REASON_NAMES}
    low = features[(features.reason == "low_gene_count") & (features.gene_index == 1)].iloc[0]
    assert low.n_cells == 2 and low.n_detected == 1 and low.total_counts == 1
    assert low.detection_fraction == 0.5 and low.mean_raw_count == 0.5
    assert all(registry.get_metadata(uri).parent_uris == contract.input_artifacts for uri in result.output_artifacts)
    assert original.uns == {"old_reference": {"answer": "secret"}}
    assert "qc_review_candidate" not in original.obs


def test_layers_counts_preferred_and_large_integer_counts_preserved(tmp_path):
    large = 2**24 + 1
    counts = np.array([[large, 0], [1, 2]], dtype=np.int64)
    data = SCData(np.ones((2, 2)), pd.DataFrame(index=["candidate", "retained"]),
                  pd.DataFrame(index=["GENE1", "GENE2"]), layers={"counts": counts})
    _, result, registry = _run(tmp_path, data)
    candidate = registry.load_scdata(result.output_artifacts[1])
    assert result.metrics["counts_source"] == "layers.counts"
    assert candidate.X.dtype == counts.dtype
    assert candidate.X[0, 0] == large
    np.testing.assert_array_equal(candidate.layers["counts"], counts[:1])


@pytest.mark.parametrize("min_genes,expected_candidates", [(0, 0), (10, 2)])
def test_empty_or_all_candidate_sets_are_descriptive_not_restored(tmp_path, min_genes, expected_candidates):
    data = SCData(np.array([[1, 2], [3, 4]]), pd.DataFrame(index=["a", "b"]),
                  pd.DataFrame(index=["G1", "G2"]))
    _, result, registry = _run(tmp_path, data, {"min_genes": min_genes})
    candidate = registry.load_scdata(result.output_artifacts[1])
    features = registry.load_payload(result.output_artifacts[3])
    assert candidate.shape == (expected_candidates, 2)
    assert result.metrics["restoration_performed"] is False
    empty = features[features.n_cells == 0]
    assert not empty.empty and empty.status.eq("no_cells").all()
    assert empty.n_detected.eq(0).all() and empty.total_counts.eq(0).all()
    assert empty.detection_fraction.isna().all() and empty.mean_raw_count.isna().all()
    assert result.metrics["candidate_cells"] + result.metrics["retained_cells"] == 2


@pytest.mark.parametrize("mito_source", ["precomputed", "annotation", "unavailable"])
def test_metric_availability_paths_and_qc_boundaries(tmp_path, mito_source):
    obs = pd.DataFrame(index=["exact", "high", "low"])
    var = pd.DataFrame(index=["GENE1", "GENE2"])
    values = np.array([[2, 8], [3, 7], [0, 1]])
    if mito_source == "precomputed":
        obs["percent_mito"] = [20.0, 30.0, 0.0]
    if mito_source == "annotation":
        var["mt"] = [True, False]
    _, result, registry = _run(tmp_path, SCData(values, obs, var))
    cells = registry.load_payload(result.output_artifacts[2])
    assert cells.low_gene_flag.tolist() == [False, False, True]
    if mito_source == "unavailable":
        assert not cells.high_mito_flag.any() and not cells.mito_pct_available.any()
        assert cells.mito_pct.isna().all()
        assert result.metrics["mitochondrial_metric_source"] == "unavailable"
    else:
        assert cells.retained.tolist() == [True, False, False]
        assert cells.high_mito_flag.tolist() == [False, True, False]
        assert cells.mito_pct_available.all()


def test_duplicate_gene_symbols_keep_distinct_gene_indices(tmp_path):
    data = SCData(np.array([[1, 2]]), pd.DataFrame(index=["a"]),
                  pd.DataFrame({"gene_name": ["DUP", "DUP"]}, index=["id1", "id2"]))
    _, result, registry = _run(tmp_path, data, {"min_genes": 0})
    table = registry.load_payload(result.output_artifacts[3])
    retained = table[table.reason == "retained"]
    assert retained.gene.tolist() == ["DUP", "DUP"]
    assert retained.gene_index.tolist() == [0, 1] and retained.gene_id.tolist() == ["id1", "id2"]
    assert retained.total_counts.tolist() == [1, 2]


@pytest.mark.parametrize("counts", [np.array([[0.5, 1.0]]), np.array([[-1, 1]]),
                                    np.array([[np.nan, 1]]), np.array([[np.inf, 1]])])
def test_invalid_raw_counts_rejected_before_outputs_are_registered(tmp_path, counts):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://review_test/input/v1"
    data = SCData(np.ones((1, 2)), pd.DataFrame(index=["a"]),
                  pd.DataFrame(index=["G1", "G2"]), layers={"counts": counts})
    registry.register(uri, data, ArtifactType.ANNDATA, "review_test", "input", "test")
    with pytest.raises(ValueError, match="finite nonnegative integer"):
        QCReviewCapability().execute(TaskContract(task_id="review", capability="qc_review", input_artifacts=[uri]), registry)
    assert registry.list_artifacts(study_id="review_test") == [registry.get_metadata(uri)]


def test_integer_reduction_overflow_rejected_before_outputs(tmp_path):
    counts = np.array([[np.iinfo(np.int64).max, 1]], dtype=np.int64)
    data = SCData(np.ones((1, 2)), pd.DataFrame(index=["a"]),
                  pd.DataFrame(index=["G1", "G2"]), layers={"counts": counts})
    with pytest.raises(ValueError, match="could overflow"):
        _run(tmp_path, data)


def test_cell_identifier_string_collision_rejected(tmp_path):
    data = SCData(np.ones((2, 2)), pd.DataFrame(index=pd.Index([1, "1"], dtype=object)),
                  pd.DataFrame(index=["G1", "G2"]))
    # Use direct producer input to exercise identity checks before H5AD's own
    # conversion to strings can make an originally mixed index ambiguous.
    registry = SimpleNamespace(get_scdata=lambda uri: (SimpleNamespace(sha256_hash="sha256:" + "0" * 64), data))
    with pytest.raises(ValueError, match="collide"):
        QCReviewCapability().execute(TaskContract(task_id="review", capability="qc_review", input_artifacts=["adata://s/raw/v1"]), registry)


@pytest.mark.parametrize("option", [{"restore_cells": True}, {"reference_label": "cell_type"},
                                    {"marker_genes": ["PF4"]}, {"min_genes": 1.5},
                                    {"max_mito_pct": np.inf}])
def test_unknown_inferential_options_and_invalid_thresholds_rejected(option):
    with pytest.raises(ValidationError):
        QCReviewParameters.model_validate(option)


def test_plan_inherits_actual_qc_source_thresholds_and_dependencies():
    qc = TaskContract(task_id="qc", capability="qc", method="sc_qc_v1",
                      input_artifacts=["adata://s/raw/v1"], depends_on=["audit"],
                      parameters={"min_genes": 17, "max_mito_pct": 25.0})
    task = _plan_qc_review(manifest=SimpleNamespace(study_id="s"), parameters={}, tasks=[qc])
    assert task.parameters == {"min_genes": 17, "max_mito_pct": 25.0}
    assert task.input_artifacts == qc.input_artifacts and task.depends_on == ["audit"]
    assert "qc" not in task.depends_on
    assert task.expected_outputs == ["json://s/qc_review/v1", "adata://s/qc_candidates/v1",
                                     "table://s/qc_review_cells/v1", "table://s/qc_review_features/v1"]
    assert QCReviewCapability().get_descriptor().scope == "shared"
    with pytest.raises(ValueError, match="must match"):
        _plan_qc_review(manifest=SimpleNamespace(study_id="s"), parameters={"min_genes": 18}, tasks=[qc])


def test_no_reference_label_dependence_in_selection_or_feature_counts(tmp_path):
    a = _data(True)
    b = _data(True)
    b.obs["cell_type"] = ["new", "x", "arbitrary", "new", "x"]
    _, result_a, registry_a = _run(tmp_path / "a", a)
    _, result_b, registry_b = _run(tmp_path / "b", b)
    for index in (2, 3):
        pd.testing.assert_frame_equal(registry_a.load_payload(result_a.output_artifacts[index]),
                                      registry_b.load_payload(result_b.output_artifacts[index]))
    assert registry_a.load_scdata(result_a.output_artifacts[1]).obs.index.equals(
        registry_b.load_scdata(result_b.output_artifacts[1]).obs.index)


@pytest.mark.parametrize("min_genes", [0, 2, 10])
def test_real_producer_outputs_pass_independent_auditor_including_empty_candidates(tmp_path, min_genes):
    from eacbp.auditor.qc_review import QCReviewValidator

    data = SCData(np.array([[1, 2], [3, 4]]), pd.DataFrame({"cell_type": ["a", "b"]}, index=["a", "b"]),
                  pd.DataFrame(index=["G1", "G2"]))
    contract, result, registry = _run(tmp_path, data, {"min_genes": min_genes})
    report = QCReviewValidator().audit(contract, result, registry)
    assert report.overall_passed, [check.message for check in report.checks]
    assert report.checks[0].metrics["outputs_verified"] == 4


def test_mixed_reasons_and_raw_count_layer_pass_independent_audit(tmp_path):
    from eacbp.auditor.qc_review import QCReviewValidator

    data = _data(True)
    data.layers["counts"] = data.X.copy()
    data.X = np.ones(data.shape, dtype=np.float32)
    contract, result, registry = _run(tmp_path, data)
    report = QCReviewValidator().audit(contract, result, registry)
    assert report.overall_passed, [check.message for check in report.checks]


def test_review_evidence_hook_returns_no_claims_without_loading_payloads():
    from eacbp.capabilities.registry import CapabilityRegistry
    from eacbp.schemas.task import TaskResult

    capabilities = CapabilityRegistry()
    capabilities.register(QCReviewCapability())
    contract = TaskContract(task_id="review", capability="qc_review", method="raw_count_qc_review_v1")
    result = TaskResult(task_id="review", capability="qc_review", method_used="raw_count_qc_review_v1",
                        status=TaskStatus.SUCCESS)
    def no_payload_reads(*args, **kwargs):
        raise AssertionError("Review-only evidence hook must not read artifact payloads")
    registry = SimpleNamespace(get=no_payload_reads, load_payload=no_payload_reads)
    report = SimpleNamespace(overall_passed=True, stop_rule_triggered=False)
    assert capabilities.extract_evidence(contract, result, report, registry) == []
