"""Raw-count QC review is audited independently after durable serialization."""

from copy import deepcopy
from pathlib import Path
import hashlib

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor import ScientificAuditor
from eacbp.auditor.qc_review import QCReviewValidator
from eacbp.capabilities.qc_review import QCReviewCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def _data(*, as_sparse=False):
    counts = np.array([[3, 2, 0, 0], [1, 0, 0, 0], [1, 1, 0, 10],
                       [0, 0, 0, 4], [0, 0, 0, 0]], dtype=np.int64)
    matrix = sparse.csr_matrix(counts) if as_sparse else counts
    return SCData(
        X=sparse.csr_matrix(counts * -2), layers={"counts": matrix},
        obs=pd.DataFrame({"donor_id": ["d1", "d1", "d2", "d2", "d3"],
                          "cell_type": ["untrusted_A", "untrusted_B", "untrusted_A", "untrusted_B", "untrusted_C"]},
                         index=[f"cell{i}" for i in range(5)]),
        var=pd.DataFrame({"gene_name": ["duplicate", "duplicate", "G3", "MT-X"]},
                         index=["feature_a", "feature_b", "feature_c", "feature_d"]),
        obsm={"X_umap": np.zeros((5, 2))}, uns={"old_clustering": {"answer": "ignored"}},
    )


def _run(tmp_path, data=None, parameters=None, *, source_type=ArtifactType.ANNDATA):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://review/raw/v1"
    registry.register(uri, data or _data(), source_type, "review", "ingest", "h5ad_import")
    outputs = ["json://review/qc_review/v1", "adata://review/qc_candidates/v1",
               "table://review/qc_review_cells/v1", "table://review/qc_review_features/v1"]
    contract = TaskContract(
        task_id="review_task", capability="qc_review", method="raw_count_qc_review_v1",
        input_artifacts=[uri], expected_outputs=outputs,
        parameters=parameters or {"min_genes": 2, "max_mito_pct": 40.0},
        allowed_operations=["summarize_qc_candidates"], validation_requirements=["qc_review_integrity"],
    )
    source_meta = registry.get_metadata(uri)
    source_hash = hashlib.sha256(Path(source_meta.storage_path).read_bytes()).hexdigest()
    result = QCReviewCapability().execute(contract, registry)
    assert registry.get_metadata(uri).sha256_hash == "sha256:" + source_hash
    assert hashlib.sha256(Path(source_meta.storage_path).read_bytes()).hexdigest() == source_hash
    return registry, contract, result


@pytest.mark.parametrize("as_sparse", [False, True])
@pytest.mark.parametrize("mito_source", ["gene_symbol_prefix", "var.mt", "obs.percent_mito", "unavailable"])
def test_four_outputs_are_independently_verified_after_reopen(tmp_path, as_sparse, mito_source):
    data = _data(as_sparse=as_sparse)
    if mito_source == "var.mt":
        data.var["mt"] = [False, False, False, True]
    elif mito_source == "obs.percent_mito":
        # Overrides the raw-derived percentage, including the zero-count row.
        data.obs["percent_mito"] = [0.0, 100.0, 0.0, 100.0, 10.0]
    elif mito_source == "unavailable":
        data.var["gene_name"] = ["g1", "g2", "g3", "g4"]
    registry, contract, result = _run(tmp_path, data)
    reopened = ArtifactRegistry(registry.storage.base_dir)
    report = QCReviewValidator().audit(contract, result, reopened)
    assert report.overall_passed, report
    assert report.checks[0].metrics["outputs_verified"] == 4
    assert result.metrics["mitochondrial_metric_source"] == mito_source
    assert ScientificAuditor().audit_task(contract, result, reopened).overall_passed


@pytest.mark.parametrize("candidate_count,min_genes,max_mito", [(0, 0, 100.0), (5, 10, 100.0)])
def test_empty_or_all_candidate_branch_is_valid_and_not_restored(tmp_path, candidate_count, min_genes, max_mito):
    registry, contract, result = _run(tmp_path, parameters={"min_genes": min_genes, "max_mito_pct": max_mito})
    assert result.metrics["candidate_cells"] == candidate_count
    assert result.metrics["restoration_performed"] is False
    assert QCReviewValidator().audit(contract, result, registry).overall_passed
    # The independent audit must run even without descriptor preparation.
    contract.validation_requirements = []
    report = ScientificAuditor().audit_task(contract, result, registry)
    assert report.overall_passed, report
    assert any(check.check_name == "qc_review_integrity" and check.passed for check in report.checks)


def _reregister(registry, contract, result, tmp_path, corruption):
    """Publish tampered artifacts normally so all payload hashes remain valid."""
    replacement = ArtifactRegistry(tmp_path / "rehashed")
    source_meta, source = registry.get_scdata(contract.input_artifacts[0])
    if corruption in {"input_fraction", "input_nan", "input_inf", "input_negative", "input_overflow"}:
        matrix = source.layers["counts"].astype(float).copy() if corruption in {
            "input_fraction", "input_nan", "input_inf"} else source.layers["counts"].copy()
        matrix[0, 0] = {"input_fraction": 0.5, "input_nan": np.nan, "input_inf": np.inf,
                        "input_negative": -1, "input_overflow": np.iinfo(np.int64).max // 2 + 1}[corruption]
        source.layers["counts"] = matrix
    elif corruption == "input_id_collision":
        source.obs.index = ["same", "same", *source.obs.index[2:]]
    replacement.register(source_meta.uri, source, source_meta.type, source_meta.study_id,
                         source_meta.created_by_task, source_meta.operation,
                         parameters=source_meta.parameters, summary_metrics=source_meta.summary_metrics)
    records = [registry.get(uri) for uri in result.output_artifacts]
    summary = deepcopy(records[0][1])
    summary["source_sha256"] = replacement.get_metadata(source_meta.uri).sha256_hash
    candidate = SCData.from_payload(records[1][1]).copy()
    cells, features = records[2][1].copy(), records[3][1].copy()
    if corruption == "source_sha":
        summary["source_sha256"] = "0" * 64
    elif corruption == "damage_diagnosis":
        summary["damage_diagnosis"] = True
    elif corruption == "summary_extra":
        summary["fabricated_p_value"] = 0.001
    candidate.uns = {"qc_review": deepcopy(summary)}
    if corruption == "candidate_counts":
        candidate.X.data[0] += 1
        candidate.layers["counts"] = candidate.X.copy()
    elif corruption == "candidate_order":
        candidate = candidate.subset_obs(np.arange(candidate.n_obs)[::-1])
    elif corruption == "candidate_id":
        candidate.obs.index = ["forged", *candidate.obs.index[1:]]
    elif corruption == "candidate_annotation":
        candidate.obs["cell_type"] = "fabricated_damage_identity"
    elif corruption == "candidate_restored_flag":
        candidate.obs["qc_passed"] = True
    elif corruption == "candidate_embedding":
        candidate.obsm["X_umap"] = np.zeros((candidate.n_obs, 2))
    elif corruption == "cell_reason":
        cells.loc[cells.index[0], "reason"] = "low_gene_count"
    elif corruption == "cell_count":
        cells.loc[cells.index[0], "total_counts"] += 1
    elif corruption == "cell_missing":
        cells = cells.iloc[:-1]
    elif corruption == "cell_mito_availability":
        cells.loc[cells.index[-1], "mito_pct_available"] = True
    elif corruption == "cell_infinite_mito":
        cells.loc[cells.index[0], "mito_pct"] = np.inf
    elif corruption == "feature_detection":
        features.loc[features.index[0], "n_detected"] += 1
    elif corruption == "feature_counts":
        features.loc[features.index[0], "total_counts"] += 1
    elif corruption == "feature_mean":
        features.loc[features.index[0], "mean_raw_count"] += 1
    elif corruption == "feature_missing":
        features = features.iloc[:-1]
    elif corruption == "feature_duplicate":
        features.iloc[-1] = features.iloc[0]
    elif corruption == "feature_gene_id":
        features.loc[features.index[0], "gene_id"] = "forged"
    elif corruption == "feature_inference":
        features["fdr_q_value"] = 0.001
    elif corruption == "feature_empty_status":
        features.loc[features["n_cells"] == 0, "status"] = "estimated_supported"
    payloads = [summary, candidate, cells, features]
    for index, ((meta, _), payload) in enumerate(zip(records, payloads)):
        parameters = deepcopy(meta.parameters)
        metadata_summary = deepcopy(summary)
        parents = meta.parent_uris
        if corruption == "artifact_threshold" and index == 3:
            parameters["min_genes"] += 1
        if corruption == "artifact_parent" and index == 3:
            parents = []
        if corruption == "artifact_extra_inference" and index == 3:
            metadata_summary["fdr_q_value"] = 0.001
        replacement.register(meta.uri, payload, meta.type, meta.study_id, meta.created_by_task,
                             meta.operation, parent_uris=parents, parameters=parameters, summary_metrics=metadata_summary)
    changed = result.model_copy(deep=True)
    changed.metrics = deepcopy(summary)
    if corruption == "omit_output":
        changed.output_artifacts = changed.output_artifacts[:-1]
    if corruption == "result_only_restoration":
        changed.metrics["restoration_performed"] = True
    if corruption == "operation_restore":
        changed.executed_operations = ["restore_cells"]
    replacement.verify_many(changed.output_artifacts)
    return replacement, changed


@pytest.mark.parametrize("corruption", [
    "source_sha", "damage_diagnosis", "summary_extra", "candidate_counts", "candidate_order", "candidate_id",
    "candidate_annotation", "candidate_restored_flag", "candidate_embedding", "cell_reason", "cell_count",
    "cell_missing", "cell_mito_availability", "cell_infinite_mito", "feature_detection", "feature_counts",
    "feature_mean", "feature_missing", "feature_duplicate", "feature_gene_id", "feature_inference",
    "artifact_threshold", "artifact_parent", "artifact_extra_inference", "omit_output", "result_only_restoration",
    "operation_restore",
])
def test_rehashed_tampering_is_rejected_for_every_output(tmp_path, corruption):
    registry, contract, result = _run(tmp_path)
    replacement, changed = _reregister(registry, contract, result, tmp_path, corruption)
    report = QCReviewValidator().audit(contract, changed, replacement)
    assert not report.overall_passed, (corruption, report)
    assert report.checks[0].check_name == "qc_review_integrity" and not report.checks[0].passed
    contract.validation_requirements = []
    assert not ScientificAuditor().audit_task(contract, changed, replacement).overall_passed


def test_rehashed_empty_group_cannot_claim_an_estimated_result(tmp_path):
    registry, contract, result = _run(tmp_path, parameters={"min_genes": 0, "max_mito_pct": 100.0})
    replacement, changed = _reregister(registry, contract, result, tmp_path, "feature_empty_status")
    assert not QCReviewValidator().audit(contract, changed, replacement).overall_passed


def test_rehashed_unmodified_control_still_passes_independent_audit(tmp_path):
    registry, contract, result = _run(tmp_path)
    replacement, changed = _reregister(registry, contract, result, tmp_path, "unmodified")
    assert QCReviewValidator().audit(contract, changed, replacement).overall_passed
    assert ScientificAuditor().audit_task(contract, changed, replacement).overall_passed


@pytest.mark.parametrize("corruption", ["input_fraction", "input_nan", "input_inf", "input_negative",
                                         "input_overflow", "input_id_collision"])
def test_independent_reconstruction_rejects_invalid_rehashed_raw_input(tmp_path, corruption):
    registry, contract, result = _run(tmp_path)
    replacement, changed = _reregister(registry, contract, result, tmp_path, corruption)
    assert not QCReviewValidator().audit(contract, changed, replacement).overall_passed


def test_large_integer_candidate_count_is_preserved_exactly(tmp_path):
    data = _data(as_sparse=True)
    data.layers["counts"][1, 0] = 2 ** 24 + 1
    registry, contract, result = _run(tmp_path, data)
    candidate = registry.get_scdata(result.output_artifacts[1])[1]
    assert candidate.X.dtype == np.int64
    assert candidate.X[0, 0] == 2 ** 24 + 1
    assert QCReviewValidator().audit(contract, result, registry).overall_passed


def test_spatial_input_is_described_without_carrying_coordinates_to_candidates(tmp_path):
    registry, contract, result = _run(tmp_path, _data(as_sparse=True),
                                       source_type=ArtifactType.SPATIAL_DATA)
    assert QCReviewValidator().audit(contract, result, registry).overall_passed
    assert not registry.get_scdata(result.output_artifacts[1])[1].obsm


def test_audit_uses_registered_source_without_reexecuting_producer(tmp_path, monkeypatch):
    registry, contract, result = _run(tmp_path)
    def prohibit(*args, **kwargs):
        raise AssertionError("Independent audit must not execute the producer")
    monkeypatch.setattr(QCReviewCapability, "execute", prohibit)
    assert QCReviewValidator().audit(contract, result, registry).overall_passed


def test_registered_origin_markers_are_preserved_as_metadata_extras(tmp_path):
    data = _data()
    data.uns.update({"is_simulated": False, "data_origin": "local test fixture"})
    registry, contract, result = _run(tmp_path, data)
    assert QCReviewValidator().audit(contract, result, registry).overall_passed
    for uri in result.output_artifacts:
        metadata = registry.get_metadata(uri)
        assert metadata.summary_metrics["is_simulated"] is False
        assert metadata.summary_metrics["data_origin"] == "local test fixture"


def test_real_executor_routing_metadata_does_not_change_scientific_summary(tmp_path):
    from eacbp.capabilities import create_default_capability_registry
    from eacbp.orchestrator.execution import TaskExecutor

    registry = ArtifactRegistry(tmp_path / "artifacts")
    uri = "adata://review/raw/v1"
    registry.register(uri, _data(), ArtifactType.ANNDATA, "review", "ingest", "h5ad_import")
    contract = TaskContract(
        task_id="review_task", capability="qc_review", method="raw_count_qc_review_v1",
        input_artifacts=[uri], parameters={"min_genes": 2, "max_mito_pct": 40.0, "study_id": "review"},
        expected_outputs=["json://review/qc_review/v1", "adata://review/qc_candidates/v1",
                          "table://review/qc_review_cells/v1", "table://review/qc_review_features/v1"],
        allowed_operations=["summarize_qc_candidates"], validation_requirements=["qc_review_integrity"],
    )
    outcome = TaskExecutor(registry, create_default_capability_registry()).execute(contract)
    try:
        assert outcome.result.metrics["execution_attempts"] == 1
        report = QCReviewValidator().audit(contract, outcome.result, outcome.staged)
        assert report.overall_passed, report
        assert ScientificAuditor().audit_task(contract, outcome.result, outcome.staged).overall_passed
        summary = outcome.staged.get(outcome.result.output_artifacts[0])[1]
        assert "execution_attempts" not in summary
        assert {key: outcome.result.metrics[key] for key in summary} == summary
    finally:
        outcome.staged.close()


@pytest.mark.parametrize("attempts", [True, 0, -1, 5, 1.5, "1"])
def test_executor_metadata_cannot_bypass_retry_contract_or_summary_schema(tmp_path, attempts):
    registry, contract, result = _run(tmp_path)
    result.metrics["execution_attempts"] = attempts
    assert not QCReviewValidator().audit(contract, result, registry).overall_passed
