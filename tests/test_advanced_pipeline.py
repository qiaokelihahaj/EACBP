import numpy as np
import pandas as pd
import pytest
pytest.importorskip("pydeseq2")
pytest.importorskip("decoupler")
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec


@pytest.mark.parametrize("with_functional", [False, True])
def test_advanced_statistics_run_through_audit_and_report(tmp_path, with_functional):
    from eacbp.report.markdown_report import ScientificReportGenerator
    rng = np.random.default_rng(14)
    counts = rng.negative_binomial(10, .4, size=(60, 20))
    counts[30:, :5] *= 3
    obs = pd.DataFrame({"donor": np.repeat([f"d{i}" for i in range(6)], 10),
                        "condition": ["A"] * 30 + ["B"] * 30})
    data = SCData(counts, obs, pd.DataFrame({"gene_name": [f"g{i}" for i in range(20)]}), layers={"counts": counts})
    reg = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "adata://advanced_pipeline/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "advanced_pipeline", "input", "input", summary_metrics={"is_simulated": True})
    manifest = StudyManifest(study_id="advanced_pipeline", biological_design=BiologicalDesign(species="human", tissue="test"), data=DataSpec(raw_artifact_uri=uri))
    state = {"method_profile": "baseline", "advanced_analysis": True,
             "analysis_extensions": {"donor_sensitivity": {"min_loo_donors": 4}},
             "capability_parameters": {"deg": {"condition_a": "A", "condition_b": "B", "n_cpus": 1}}}
    if with_functional:
        state["analysis_extensions"]["functional_activity"] = {
            "network": {"source": ["ToyPathway"] * 20, "target": [f"g{i}" for i in range(20)],
                        "weight": [1.] * 5 + [-1.] * 15},
            "network_source": "synthetic_test_fixture", "network_version": "1", "tmin": 5,
        }
    orch = ScientificOrchestrator(reg)
    result = orch.run_study(manifest, state)
    assert result["status"] == "success", result["failures"]
    assert next(t for t in orch.task_history if t.capability == "deg").method_used == "pydeseq2_pseudobulk_v1"
    sensitivity = next(t for t in orch.task_history if t.capability == "donor_sensitivity")
    assert sensitivity.metrics["skipped"]
    if with_functional:
        functional = next(t for t in orch.task_history if t.capability == "functional_activity")
        assert functional.method_used == "decoupler_ulm_v2"
    report = ScientificReportGenerator(manifest, orch.evidence_graph, reg, orch.task_history, orch.audit_reports).generate_markdown()
    assert "Advanced analysis results" in report and "skip_reason" in report
    assert "SIMULATED DATA" in report
    resumed = orch.run_study(manifest, {**state, "resume": True})
    assert resumed["status"] == "success", resumed["failures"]
