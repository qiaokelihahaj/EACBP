"""Execute the actual optional scientific libraries on controlled small matrices."""
import numpy as np
import pandas as pd
import pytest
pytest.importorskip("scanpy")
pytest.importorskip("harmonypy")
pytest.importorskip("leidenalg")
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.real_methods import HarmonyIntegrationCapability, LeidenClusteringCapability, DPTTrajectoryCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def setup_data(tmp_path):
    rng = np.random.default_rng(4)
    t = np.linspace(0, 3, 80)
    x = np.abs(t[:, None] * np.linspace(.2, 2, 12)[None, :] + rng.normal(0, .15, (80, 12)))
    data = SCData(x, pd.DataFrame({"cell_id": [f"c{i}" for i in range(80)], "batch": ["a", "b"] * 40,
                                  "cell_type_ground_truth": ["DO_NOT_USE"] * 80}),
                  pd.DataFrame({"gene_name": [f"g{i}" for i in range(12)]}))
    reg = ArtifactRegistry(str(tmp_path))
    uri = "adata://real/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "real", "input", "input")
    return reg, uri


def test_actual_harmony_leiden_umap_preserve_cells_and_provenance(tmp_path):
    reg, uri = setup_data(tmp_path)
    harmony = HarmonyIntegrationCapability()
    task = TaskContract(task_id="h", capability="integration", input_artifacts=[uri], parameters={"random_seed": 7, "n_components": 5})
    result = harmony.execute(task, reg)
    meta, data = reg.get(result.output_artifacts[0])
    assert result.method_used == "harmonypy_v1"
    assert "harmonypy" in meta.software_versions
    assert data.obsm["X_pca_harmony"].shape == (80, 5)
    assert not np.allclose(data.obsm["X_pca"], data.obsm["X_pca_uncorrected"])
    task = TaskContract(task_id="l", capability="clustering", input_artifacts=result.output_artifacts, parameters={"n_neighbors": 10})
    clustered = LeidenClusteringCapability().execute(task, reg)
    _, data = reg.get(clustered.output_artifacts[0])
    assert data.n_obs == 80 and data.obsp["connectivities"].shape == (80, 80)
    assert data.obsm["X_umap"].shape == (80, 2) and np.isfinite(data.obsm["X_umap"]).all()
    assert data.obs["cluster"].nunique() > 1
    assert "DO_NOT_USE" not in data.obs["cell_type"].tolist()


def test_actual_diffusion_pseudotime_requires_root_and_saves_per_cell_output(tmp_path):
    reg, uri = setup_data(tmp_path)
    cap = DPTTrajectoryCapability()
    task = TaskContract(task_id="d", capability="trajectory_inference", input_artifacts=[uri], parameters={"n_neighbors": 12})
    with pytest.raises(ValueError, match="root_cell_id"):
        cap.execute(task, reg)
    task.parameters["root_cell_id"] = "c0"
    result = cap.execute(task, reg)
    _, data = reg.get(result.output_artifacts[1])
    pt = data.obs["dpt_pseudotime"].to_numpy()
    assert pt[0] == 0 and np.isfinite(pt).all() and pt[-1] > .5
    _, table = reg.get(result.output_artifacts[0])
    assert table.fdr_q_value.between(0, 1).all()
    assert result.metrics["stability_evaluated"] is False


def test_standard_profile_executes_real_libraries_through_audit_and_resume(tmp_path):
    from eacbp.orchestrator.loop import ScientificOrchestrator
    from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec
    reg, uri = setup_data(tmp_path)
    manifest = StudyManifest(study_id="real", biological_design=BiologicalDesign(species="human", tissue="test"), data=DataSpec(raw_artifact_uri=uri))
    orch = ScientificOrchestrator(artifact_registry=reg)
    state = {"method_profile": "standard", "capability_parameters": {"trajectory_inference": {"root_cell_id": "c0", "run_paga": True}}}
    result = orch.run_study(manifest, state)
    assert result["status"] == "success", result["failures"]
    methods = {t.method_used for t in orch.task_history}
    assert {"harmonypy_v1", "scanpy_leiden_umap_v1", "scanpy_dpt_v1"} <= methods
    _, data = reg.get("adata://real/diffusion_pseudotime/v1")
    assert "paga" in data.uns and data.uns["paga"]["connectivities"].shape[0] >= 2
    resumed = orch.run_study(manifest, {**state, "resume": True})
    assert resumed["status"] == "success", resumed["failures"]


def test_dpt_to_cellrank_is_connected_in_study_plan(tmp_path):
    pytest.importorskip("cellrank")
    from eacbp.orchestrator.loop import ScientificOrchestrator
    from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec
    reg, uri = setup_data(tmp_path)
    manifest = StudyManifest(study_id="real", biological_design=BiologicalDesign(species="human", tissue="test"), data=DataSpec(raw_artifact_uri=uri))
    state = {"method_profile": "standard", "capability_parameters": {"trajectory_inference": {"root_cell_id": "c0"}},
             "cellrank_terminal_states": {"end1": ["c78"], "end2": ["c79"]}}
    result = ScientificOrchestrator(artifact_registry=reg).run_study(manifest, state)
    assert result["status"] == "success", result["failures"]
    _, probabilities = reg.get("table://real/fate_probabilities/v1")
    assert np.allclose(probabilities.sum(axis=1), 1, atol=1e-3)
