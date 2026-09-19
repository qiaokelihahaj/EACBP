"""Real Scrublet + local CellTypist integration through the study orchestrator."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")
pytest.importorskip("celltypist")


def test_real_qc_extensions_preserve_input_and_emit_audited_neutral_evidence(
    tmp_path, monkeypatch
):
    """Run both real QC extensions through independent audit and reporting."""

    import celltypist

    from eacbp.artifact.registry import ArtifactRegistry
    from eacbp.capabilities.sc_data import SCData
    from eacbp.orchestrator.loop import ScientificOrchestrator
    from eacbp.report.markdown_report import ScientificReportGenerator
    from eacbp.schemas.artifact import ArtifactType
    from eacbp.schemas.evidence import EvidencePolarity, EvidenceType
    from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest

    rng = np.random.default_rng(21)
    n_cells, n_genes = 64, 30
    counts = rng.poisson(3, (n_cells, n_genes)).astype(np.int32)
    counts[:32, :10] += 10
    counts[32:, 10:20] += 10
    counts += 1
    normalized = np.log1p(counts / counts.sum(axis=1)[:, None] * 10000)
    genes = [f"g{i}" for i in range(n_genes)]

    # Leave four query cells outside the SGD training rows so the real
    # probability-match mode retains Unassigned observations at a high
    # threshold.  This exercises uncertainty rather than manufacturing labels.
    model = celltypist.train(
        X=normalized[:60],
        labels=["SyntheticA"] * 30 + ["SyntheticB"] * 30,
        genes=genes,
        n_jobs=1,
        max_iter=100,
        check_expression=True,
        use_SGD=True,
        random_state=6,
    )
    model_path = tmp_path / "synthetic_sgd_model.pkl"
    model.write(str(model_path))

    # A local model path is part of the contract.  Fail if CellTypist tries
    # named-model discovery/download while this integration test is running.
    def reject_remote_discovery(*args, **kwargs):
        raise AssertionError("QC integration must use the supplied local model")

    monkeypatch.setattr(celltypist.models, "get_all_models", reject_remote_discovery)

    cell_ids = [f"cell_{i}" for i in range(n_cells)]
    obs = pd.DataFrame(
        {
            "cell_id": cell_ids,
            "batch": ["library_a"] * 32 + ["library_b"] * 32,
            # Keep conditions in disjoint donors for the baseline Welch branch.
            "donor": [f"donor_{i // 8}" for i in range(n_cells)],
            "condition": ["A"] * 32 + ["B"] * 32,
            "percent_mito": np.zeros(n_cells),
        },
        index=cell_ids,
    )
    var = pd.DataFrame({"gene_name": genes}, index=genes)
    source = SCData(counts, obs, var, layers={"counts": counts.copy()})

    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    study_id = "real_qc_flow"
    raw_uri = f"adata://{study_id}/raw/v1"
    registry.register(
        raw_uri,
        source,
        ArtifactType.ANNDATA,
        study_id,
        "input",
        "input",
        summary_metrics={"is_simulated": True},
    )
    manifest = StudyManifest(
        study_id=study_id,
        biological_design=BiologicalDesign(species="synthetic", tissue="synthetic"),
        data=DataSpec(raw_artifact_uri=raw_uri),
    )
    state = {
        "method_profile": "baseline",
        "analysis_extensions": {
            "doublet_detection": {
                "batch_key": "batch",
                "n_prin_comps": 8,
                "random_seed": 13,
                "filter_doublets": False,
            },
            "cell_annotation": {
                "model_path": str(model_path),
                "model_species": "synthetic",
                "model_tissue": "synthetic",
                "mode": "prob match",
                "p_thres": 0.9999,
                "use_as_cell_type": True,
            },
        },
    }

    orchestrator = ScientificOrchestrator(registry)
    result = orchestrator.run_study(manifest, state)
    assert result["status"] == "success", result["failures"]

    doublet_task = next(
        task for task in orchestrator.task_history if task.capability == "doublet_detection"
    )
    annotation_task = next(
        task for task in orchestrator.task_history if task.capability == "cell_annotation"
    )
    assert doublet_task.method_used == "scanpy_scrublet_v1"
    assert annotation_task.method_used == "celltypist_local_v1"
    assert doublet_task.metrics["n_cells_before"] == n_cells
    assert doublet_task.metrics["n_cells_after"] == n_cells
    assert annotation_task.metrics["n_unknown_annotations"] > 0
    assert annotation_task.metrics["n_label_conflicts"] > 0

    _, doublet_output = registry.get(doublet_task.output_artifacts[0])
    _, clustered_output = registry.get(f"adata://{study_id}/annotated/v4")
    _, annotation_output = registry.get(annotation_task.output_artifacts[0])
    for output in (doublet_output, annotation_output):
        np.testing.assert_array_equal(output.layers["counts"], counts)
        assert output.obs.index.tolist() == cell_ids
        assert output.obs["cell_id"].tolist() == cell_ids
        assert output.var.index.tolist() == genes

    # The clustering capability always emits cell_type.  With conflicting
    # model labels, use_as_cell_type=True retains that existing column and
    # keeps the reference labels in the separate CellTypist column.
    pd.testing.assert_series_equal(
        annotation_output.obs["cell_type"], clustered_output.obs["cell_type"]
    )
    pd.testing.assert_series_equal(
        annotation_output.obs["cell_type_existing"],
        clustered_output.obs["cell_type"],
        check_names=False,
    )
    assert annotation_output.obs["cell_type_celltypist_unknown"].any()
    assert annotation_output.obs["cell_type_celltypist_conflict"].any()

    audits = {
        report.target_task_id: report
        for report in orchestrator.audit_reports
        if report.target_task_id in {doublet_task.task_id, annotation_task.task_id}
    }
    assert set(audits) == {doublet_task.task_id, annotation_task.task_id}
    assert all(report.overall_passed for report in audits.values())
    assert all(
        check.passed
        for report in audits.values()
        for check in report.checks
        if check.check_name == "advanced_qc_integrity"
    )

    evidence = {
        node.source_task_id: node
        for node in orchestrator.evidence_graph.evidence_nodes.values()
        if node.source_task_id in {doublet_task.task_id, annotation_task.task_id}
    }
    assert set(evidence) == {doublet_task.task_id, annotation_task.task_id}
    assert evidence[doublet_task.task_id].type == EvidenceType.QC_METRICS
    assert evidence[annotation_task.task_id].type == EvidenceType.CELL_ANNOTATION
    for node in evidence.values():
        assert node.audit_passed is True
        assert node.score == 0
        assert node.polarity == EvidencePolarity.NEUTRAL
        assert node.data_origin_uris == [raw_uri]

    report = ScientificReportGenerator(
        manifest,
        orchestrator.evidence_graph,
        registry,
        orchestrator.task_history,
        orchestrator.audit_reports,
    ).generate_markdown()
    assert (
        f"| `n_unknown_annotations` | {annotation_task.metrics['n_unknown_annotations']} |"
        in report
    )
    assert (
        f"| `n_label_conflicts` | {annotation_task.metrics['n_label_conflicts']} |"
        in report
    )
    assert (
        f"| `n_doublets_marked` | {doublet_task.metrics['n_doublets_marked']} |"
        in report
    )
    assert "Unknown and conflicting CellTypist labels are retained as uncertainty" in report
