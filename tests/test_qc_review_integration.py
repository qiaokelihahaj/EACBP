"""Candidate preservation is a separate, admitted branch of native QC."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.sc_data import SCData
from eacbp.evidence.graph import EvidenceGraph
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.orchestrator.planning import build_study_tasks
from eacbp.report.markdown_report import ScientificReportGenerator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskResult, TaskStatus


def manifest(targets=()):
    return StudyManifest(study_id="review_flow", biological_design=BiologicalDesign(
        species="human", tissue="blood", target_cell_types=list(targets)),
        data={"raw_artifact_uri": "adata://review_flow/raw/v1"})


def config(min_genes=2):
    return {"method_profile": "baseline", "autonomous_analysis": {"intent": "exploration"},
            "analysis_extensions": {"qc_review": True},
            "capability_parameters": {"qc": {"min_genes": min_genes, "max_mito_pct": 100.0},
                                      "normalization": {"n_top_genes": 15},
                                      "integration": {"n_components": 5}}}


def source(registry):
    values = np.random.default_rng(19).poisson(3, (60, 30)).astype(np.int64)
    values[:3] = 0
    values[3:30, :5] += 20
    values[30:, 5:10] += 20
    data = SCData(sparse.csr_matrix(values), pd.DataFrame(index=[f"c{i}" for i in range(60)]),
                  pd.DataFrame(index=[f"g{i}" for i in range(30)]))
    registry.register(manifest().data.raw_artifact_uri, data, ArtifactType.ANNDATA,
                      "review_flow", "input", "import")
    return values


@pytest.mark.parametrize("targets", [(), ("T cells", "B cells")])
def test_review_is_shared_before_qc_and_never_feeds_statistics(targets):
    tasks = build_study_tasks(manifest(targets), config(), create_default_capability_registry())
    review = [task for task in tasks if task.capability == "qc_review"]
    assert len(review) == 1
    review = review[0]
    qc = next(task for task in tasks if task.capability == "qc")
    assert review.task_id in qc.depends_on
    assert qc.task_id not in review.depends_on
    assert review.input_artifacts == qc.input_artifacts
    assert {key: review.parameters[key] for key in ("min_genes", "max_mito_pct")} == {
        key: qc.parameters[key] for key in ("min_genes", "max_mito_pct")}
    assert tasks.index(review) < tasks.index(qc)
    assert all(not set(review.expected_outputs).intersection(task.input_artifacts)
               for task in tasks if task.task_id != review.task_id)


@pytest.mark.parametrize("min_genes,candidates,main_success", [(0, 0, True), (2, 3, True), (31, 60, False)])
def test_native_candidate_admission_survives_empty_and_all_filtered_inputs(tmp_path, min_genes, candidates, main_success):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    counts = source(registry)
    before = registry.get_metadata(manifest().data.raw_artifact_uri).sha256_hash
    orchestrator = ScientificOrchestrator(registry)
    summary = orchestrator.run_study(manifest(), config(min_genes))
    assert (summary["status"] == "success") is main_success, summary["failures"]
    review = next(task for task in orchestrator.task_history if task.capability == "qc_review")
    assert review.status == TaskStatus.SUCCESS
    receipt = next(record for record in registry.list_audit_records() if record.task_id == review.task_id)
    assert receipt.overall_passed
    outputs = registry.get_audited_many(receipt.signature, receipt.contract)
    assert len(outputs) == 4
    candidate = SCData.from_payload(outputs[1][1])
    assert candidate.n_obs == candidates and candidate.n_vars == 30
    expected = counts if candidates == 60 else counts[:candidates]
    assert np.array_equal(candidate.layers["counts"].toarray(), expected)
    if not main_success:
        qc = next(task for task in orchestrator.task_history if task.capability == "qc")
        assert qc.status == TaskStatus.SCIENTIFIC_FAILURE
        assert all(task.status == TaskStatus.BLOCKED for task in orchestrator.task_history
                   if task.capability in {"normalization", "integration", "clustering"})
    else:
        for task in orchestrator.task_history:
            if task.capability in {"normalization", "integration", "clustering"}:
                data = registry.get_scdata(task.output_artifacts[0])[1]
                assert not set(data.obs.index).intersection(candidate.obs.index)
    report = ScientificReportGenerator(manifest(), orchestrator.evidence_graph, registry,
                                        orchestrator.task_history, orchestrator.audit_reports)._qc_review_details()
    text = "\n".join(report)
    assert "Excluded-cell QC review" in text and "review unavailable" not in text
    assert f"excluded candidates: {candidates}" in text
    assert "not cell-type assignments" in text
    assert registry.get_metadata(manifest().data.raw_artifact_uri).sha256_hash == before


def test_candidate_review_resumes_without_recomputation(tmp_path, monkeypatch):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    source(registry)
    first = ScientificOrchestrator(registry)
    assert first.run_study(manifest(), config())["status"] == "success"
    review = next(task for task in first.task_history if task.capability == "qc_review")
    hashes = {uri: registry.get_metadata(uri).sha256_hash for uri in review.output_artifacts}
    from eacbp.capabilities.qc_review import QCReviewCapability
    def forbidden(*args, **kwargs):
        raise AssertionError("Admitted candidate review must resume")
    monkeypatch.setattr(QCReviewCapability, "execute", forbidden)
    second = ScientificOrchestrator(registry)
    resumed = config()
    resumed["resume"] = True
    assert second.run_study(manifest(), resumed)["status"] == "success"
    assert all(registry.get_metadata(uri).sha256_hash == digest for uri, digest in hashes.items())


def test_review_disabled_preserves_primary_scientific_artifacts(tmp_path):
    tables = []
    for enabled in (False, True):
        registry = ArtifactRegistry(tmp_path / str(enabled))
        source(registry)
        state = config()
        state["analysis_extensions"]["qc_review"] = enabled
        orchestrator = ScientificOrchestrator(registry)
        assert orchestrator.run_study(manifest(), state)["status"] == "success"
        outputs = {}
        for task in orchestrator.task_history:
            if task.capability in {"qc", "normalization", "integration", "clustering"}:
                data = registry.get_scdata(task.output_artifacts[0])[1]
                matrix = data.X.toarray() if sparse.issparse(data.X) else np.asarray(data.X)
                layers = {name: value.toarray() if sparse.issparse(value) else np.asarray(value)
                          for name, value in data.layers.items()}
                outputs[task.capability] = (data.obs.copy(), matrix, deepcopy(data.obsm), layers)
        tables.append(outputs)
    for capability in tables[0]:
        old, new = tables[0][capability], tables[1][capability]
        pd.testing.assert_frame_equal(old[0], new[0])
        assert np.array_equal(old[1], new[1])
        assert old[2].keys() == new[2].keys()
        assert all(np.array_equal(old[2][name], new[2][name]) for name in old[2])
        assert old[3].keys() == new[3].keys()
        assert all(np.array_equal(old[3][name], new[3][name]) for name in old[3])


def test_report_does_not_present_a_successful_unadmitted_review(tmp_path):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    task = TaskResult(task_id="review", capability="qc_review", method_used="raw_count_qc_review_v1",
                      status=TaskStatus.SUCCESS, output_artifacts=["json://review_flow/qc_review/v1"])
    text = "\n".join(ScientificReportGenerator(manifest(), EvidenceGraph(), registry, [task])._qc_review_details())
    assert "review unavailable" in text
    assert "excluded candidates:" not in text
