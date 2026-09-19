"""Real local LIANA smoke and donor comparison regressions."""
import hashlib
import numpy as np
import pandas as pd
import pytest
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.advanced_communication import LianaValidator, _safe_rank_comparison_pvalue
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.advanced_communication import (
    LianaCommunicationCapability,
    LianaCommunicationInputError,
)
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def _communication_fixture(*, genes=None, spatial=None):
    genes = genes or ["L1", "R1", "G1"]
    counts = np.ones((4, len(genes)), dtype=float)
    obs = pd.DataFrame({
        "donor_id": ["d1"] * 4,
        "condition": ["A"] * 4,
        "cell_type": ["T", "T", "B", "B"],
    })
    obsm = {"spatial": spatial} if spatial is not None else {}
    return SCData(counts, obs, pd.DataFrame({"gene_name": genes}), obsm=obsm)


def _communication_task(tmp_path, data, resource_frame, *, species="human", spatial=False):
    resource = tmp_path / "lr.csv"
    resource_frame.to_csv(resource, index=False)
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://communication/fixture-input/v1"
    registry.register(input_uri, data, ArtifactType.ANNDATA, "communication", "input", "input")
    parameters = {
        "lr_resource_path": str(resource),
        "species": species,
        "lr_resource_version": "fixture-v1",
        "lr_resource_source": "test fixture",
        "min_cells": 2,
        "n_perms": 0,
    }
    if spatial:
        parameters.update({"spatial_key": "spatial", "use_spatial": True})
    task = TaskContract(
        task_id="liana_fixture",
        capability="liana_communication",
        input_artifacts=[input_uri],
        expected_outputs=["table://communication/fixture-ranks/v1"],
        parameters=parameters,
    )
    return registry, task


def _fake_rank_aggregate(*args, **kwargs):
    return pd.DataFrame({
        "source": ["T"],
        "target": ["B"],
        "ligand": ["L1"],
        "receptor": ["R1"],
        "magnitude_rank": [0.2],
        "specificity_rank": [0.3],
    })


@pytest.mark.parametrize("pipeline", [False, True])
def test_real_liana_local_resource_per_donor(tmp_path, pipeline):
    pytest.importorskip("liana")
    rng = np.random.default_rng(4)
    genes = ["L1", "L2", "R1", "R2", "G1", "G2"]
    counts = rng.poisson(3, (96, 6)) + 1
    obs = pd.DataFrame({"donor_id": np.repeat(["d1", "d2", "d3", "d4"], 24),
                        "condition": ["A"] * 48 + ["B"] * 48,
                        "cell_type": (["T"] * 12 + ["B"] * 12) * 4})
    data = SCData(np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4),
                  obs, pd.DataFrame({"gene_name": genes}), layers={"counts": counts})
    resource = tmp_path / "lr.csv"
    pd.DataFrame({"ligand": ["L1", "L2"], "receptor": ["R1", "R2"]}).to_csv(resource, index=False)
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://communication/normalized/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "communication", "input", "input")
    task = TaskContract(task_id="liana", capability="liana_communication", input_artifacts=[uri],
                        expected_outputs=["table://communication/ranks/v1", "table://communication/comparison/v1"],
                        parameters={"lr_resource_path": str(resource), "species": "human",
                                    "lr_resource_version": "synthetic-v1", "lr_resource_source": "test fixture",
                                    "n_perms": 10, "min_cells": 5})
    result = LianaCommunicationCapability().execute(task, registry)
    _, table = registry.get(result.output_artifacts[0])
    assert set(table.donor_id) == {"d1", "d2", "d3", "d4"}
    assert table.groupby("donor_id").size().to_dict() == {"d1": 8, "d2": 8, "d3": 8, "d4": 8}
    assert not table.rank_is_fdr.any()
    assert table[["magnitude_rank", "specificity_rank"]].apply(lambda column: column.between(0, 1).all()).all()
    from eacbp.auditor import ScientificAuditor
    task.validation_requirements = ["liana_integrity"]
    audit = ScientificAuditor().audit_task(task, result, registry)
    assert audit.overall_passed, audit.model_dump()
    from eacbp.auditor.advanced_communication import LianaValidator
    from unittest.mock import patch
    original_get = registry.get
    def corrupted_comparison(uri):
        meta, payload = original_get(uri)
        if uri == result.output_artifacts[1]:
            payload = payload.copy()
            payload.loc[0, "comparison_fdr_magnitude"] = -1
        return meta, payload
    with patch.object(registry, "get", side_effect=corrupted_comparison):
        assert not LianaValidator().audit(task, result, registry).overall_passed
    if pipeline:
        from eacbp.orchestrator.loop import ScientificOrchestrator
        from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec
        from eacbp.report.markdown_report import ScientificReportGenerator
        # Use independently supplied labels, preserved by preprocessing.
        data.obs["reference_type"] = data.obs.cell_type
        # Separate input URI: registered payloads are immutable.
        pipeline_uri = "adata://communication/pipeline_input/v1"
        registry.register(pipeline_uri, data, ArtifactType.ANNDATA, "communication", "pipeline_input", "input")
        manifest = StudyManifest(study_id="communication", biological_design=BiologicalDesign(species="human", tissue="test"), data=DataSpec(raw_artifact_uri=pipeline_uri))
        state = {"method_profile": "baseline", "analysis_extensions": {
            "liana_communication": {**task.parameters, "cell_type_col": "reference_type"}},
            "capability_parameters": {"qc": {"min_genes": 1, "min_cells": 1}}}
        orch = ScientificOrchestrator(registry)
        outcome = orch.run_study(manifest, state)
        assert outcome["status"] == "success", outcome["failures"]
        report_text = ScientificReportGenerator(manifest, orch.evidence_graph, registry, orch.task_history, orch.audit_reports).generate_markdown()
        assert "LIANA ranks are not FDR" in report_text
        assert "liana_condition_comparison" in report_text
        resumed = orch.run_study(manifest, {**state, "resume": True})
        assert resumed["status"] == "success", resumed["failures"]


def test_paired_comparison_means_exclude_unpaired_donors():
    rows = pd.DataFrame({"source": "T", "target": "B", "ligand": "L", "receptor": "R",
                         "donor_id": ["d1", "d2", "d1", "d2", "extra"],
                         "condition": ["A", "A", "B", "B", "A"],
                         "magnitude_rank": [.1, .2, .3, .5, .99],
                         "specificity_rank": [.1, .2, .3, .5, .99]})
    table, metrics = LianaCommunicationCapability._comparison_table(rows, {
        "condition_a": "A", "condition_b": "B", "paired": True, "min_donors": 2,
        "donors_a": {"d1", "d2", "extra"}, "donors_b": {"d1", "d2"}})
    assert table.iloc[0].mean_magnitude_rank_a == pytest.approx(.15)
    assert table.iloc[0].n_donors_a == 2
    assert metrics["comparison_fdr_applied"]
    assert table.iloc[0].comparison_fdr_magnitude == pytest.approx(table.iloc[0].comparison_p_value_magnitude)


def test_rank_comparison_rejects_degenerate_variance_as_not_estimable():
    """Auditor math agrees with the executor for zero-variance donor ranks."""

    assert np.isnan(_safe_rank_comparison_pvalue([0.1, 0.1], [0.2, 0.2], paired=False))
    assert np.isnan(_safe_rank_comparison_pvalue([0.1, 0.2], [0.2, 0.3], paired=True))
    assert np.isfinite(_safe_rank_comparison_pvalue([0.1, 0.2], [0.2, 0.4], paired=False))


def test_auditor_accepts_paired_not_estimated_comparison(tmp_path):
    """A constant paired difference is retained and audited as not_estimated."""

    data = SCData(
        np.ones((8, 2), dtype=float),
        pd.DataFrame({
            "donor_id": ["d1", "d1", "d1", "d1", "d2", "d2", "d2", "d2"],
            "condition": ["A", "A", "B", "B", "A", "A", "B", "B"],
            "cell_type": ["T", "B"] * 4,
        }),
        pd.DataFrame({"gene_name": ["L1", "R1"]}),
    )
    resource = tmp_path / "lr.csv"
    pd.DataFrame({"ligand": ["L1"], "receptor": ["R1"]}).to_csv(resource, index=False)
    digest = hashlib.sha256(resource.read_bytes()).hexdigest()
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://communication/degenerate-input/v1"
    rank_uri = "table://communication/degenerate-ranks/v1"
    comparison_uri = "table://communication/degenerate-comparison/v1"
    registry.register(input_uri, data, ArtifactType.ANNDATA, "communication", "input", "input")
    params = {
        "species": "human",
        "donor_col": "donor_id",
        "condition_col": "condition",
        "cell_type_col": "cell_type",
        "condition_a": "A",
        "condition_b": "B",
        "paired": True,
        "min_donors_for_comparison": 2,
        "spatial_key": None,
        "resource": {
            "resource_path": str(resource.resolve()),
            "resource_sha256": digest,
            "resource_version": "fixture-v1",
            "resource_source": "test fixture",
            "species": "human",
            "resource_species": None,
            "resource_species_columns": [],
            "resource_species_values": [],
        },
        "rank_aggregate": {"return_all_lrs": True},
    }
    ranks = pd.DataFrame({
        "donor_id": ["d1", "d1", "d2", "d2"],
        "condition": ["A", "B", "A", "B"],
        "donor_condition": ["d1::A", "d1::B", "d2::A", "d2::B"],
        "n_cells_group": [2, 2, 2, 2],
        "source": ["T"] * 4,
        "target": ["B"] * 4,
        "ligand": ["L1"] * 4,
        "receptor": ["R1"] * 4,
        "magnitude_rank": [0.1, 0.2, 0.1, 0.2],
        "specificity_rank": [0.3, 0.4, 0.3, 0.4],
        "species": ["human"] * 4,
        "resource_version": ["fixture-v1"] * 4,
        "resource_source": ["test fixture"] * 4,
        "resource_sha256": [digest] * 4,
        "rank_is_fdr": [False] * 4,
        "spatial_used": [False] * 4,
    })
    comparison = pd.DataFrame({
        "source": ["T"], "target": ["B"], "ligand": ["L1"], "receptor": ["R1"],
        "condition_a": ["A"], "condition_b": ["B"], "n_donors_a": [2], "n_donors_b": [2],
        "comparison_test": ["paired_t_test"], "comparison_status": ["not_estimated"],
        "mean_magnitude_rank_a": [0.1], "mean_magnitude_rank_b": [0.2],
        "mean_specificity_rank_a": [0.3], "mean_specificity_rank_b": [0.4],
        "comparison_p_value_magnitude": [np.nan], "comparison_p_value_specificity": [np.nan],
        "comparison_fdr_magnitude": [np.nan], "comparison_fdr_specificity": [np.nan],
    })
    registry.register(rank_uri, ranks, ArtifactType.TABLE, "communication", "degenerate", "liana_rank_aggregate_donor_condition", parameters=params)
    registry.register(comparison_uri, comparison, ArtifactType.TABLE, "communication", "degenerate", "liana_rank_condition_comparison_donor_level", parameters=params)
    contract = TaskContract(
        task_id="degenerate",
        capability="liana_communication",
        input_artifacts=[input_uri],
        expected_outputs=[rank_uri, comparison_uri],
        parameters={"species": "human"},
    )
    result = TaskResult(
        task_id="degenerate", status=TaskStatus.SUCCESS,
        capability="liana_communication", method_used="liana_rank_aggregate_v1",
        input_artifacts=[input_uri], output_artifacts=[rank_uri, comparison_uri],
    )
    report = LianaValidator().audit(contract, result, registry)
    assert report.overall_passed, report.model_dump()


def test_missing_gene_pairs_are_recorded_and_independently_audited(tmp_path, monkeypatch):
    data = _communication_fixture()
    resource = pd.DataFrame({
        "ligand": ["L1", "L2"],
        "receptor": ["R1", "R2"],
    })
    registry, task = _communication_task(tmp_path, data, resource)
    monkeypatch.setattr(
        "eacbp.capabilities.advanced_communication._load_liana",
        lambda: type("LI", (), {"mt": type("MT", (), {"rank_aggregate": _fake_rank_aggregate})()})(),
    )

    result = LianaCommunicationCapability().execute(task, registry)
    missing = result.metrics["group_missing_non_evaluable"]
    assert missing == {
        "d1::A": [
            {"ligand": "L2", "receptor": "R2", "missing_genes": ["L2", "R2"]}
        ]
    }
    _, metadata_payload = registry.get(result.output_artifacts[0])
    assert set(metadata_payload.columns) >= {"ligand", "receptor"}
    audit = LianaValidator().audit(task, result, registry)
    assert audit.overall_passed, audit.model_dump()


def test_resource_species_mismatch_is_rejected_before_inference(tmp_path):
    data = _communication_fixture()
    resource = pd.DataFrame({
        "ligand": ["L1"],
        "receptor": ["R1"],
        "species": ["mouse"],
    })
    registry, task = _communication_task(tmp_path, data, resource)

    with pytest.raises(LianaCommunicationInputError, match="do not match requested species"):
        LianaCommunicationCapability().execute(task, registry)


def test_invalid_spatial_coordinates_are_rejected(tmp_path, monkeypatch):
    data = _communication_fixture(spatial=np.ones((4, 1), dtype=float))
    resource = pd.DataFrame({"ligand": ["L1"], "receptor": ["R1"]})
    registry, task = _communication_task(tmp_path, data, resource, spatial=True)
    monkeypatch.setattr(
        "eacbp.capabilities.advanced_communication._load_liana",
        lambda: type("LI", (), {"mt": type("MT", (), {"rank_aggregate": _fake_rank_aggregate})()})(),
    )

    with pytest.raises(LianaCommunicationInputError, match="spatial coordinates"):
        LianaCommunicationCapability().execute(task, registry)


def test_real_liana_missing_gene_fixture_records_non_evaluable_pair(tmp_path):
    pytest.importorskip("liana")
    rng = np.random.default_rng(14)
    genes = ["L1", "R1", "G1", "G2"]
    counts = rng.poisson(3, (24, len(genes))).astype(float) + 1
    obs = pd.DataFrame({
        "donor_id": ["d1"] * 24,
        "condition": ["A"] * 24,
        "cell_type": (["T"] * 12) + (["B"] * 12),
    })
    data = SCData(
        np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4),
        obs,
        pd.DataFrame({"gene_name": genes}),
        layers={"counts": counts},
    )
    resource = tmp_path / "lr_missing_gene.csv"
    pd.DataFrame({
        "ligand": ["L1", "Lmissing"],
        "receptor": ["R1", "Rmissing"],
    }).to_csv(resource, index=False)
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    input_uri = "adata://communication/real-missing-input/v1"
    registry.register(input_uri, data, ArtifactType.ANNDATA, "communication", "input", "input")
    task = TaskContract(
        task_id="real_liana_missing_gene",
        capability="liana_communication",
        input_artifacts=[input_uri],
        expected_outputs=["table://communication/real-missing-ranks/v1"],
        parameters={
            "lr_resource_path": str(resource),
            "species": "human",
            "lr_resource_version": "fixture-v1",
            "lr_resource_source": "test fixture",
            "min_cells": 5,
            "n_perms": 10,
        },
    )
    result = LianaCommunicationCapability().execute(task, registry)
    assert result.metrics["group_missing_non_evaluable"]["d1::A"] == [
        {"ligand": "Lmissing", "receptor": "Rmissing", "missing_genes": ["Lmissing", "Rmissing"]}
    ]
    audit = LianaValidator().audit(task, result, registry)
    assert audit.overall_passed, audit.model_dump()
