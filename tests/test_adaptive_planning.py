import numpy as np
import pandas as pd
import pytest
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.subset import SubsetCapability
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.schemas.study import StudyManifest, BiologicalDesign, DataSpec
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def test_missing_target_never_substitutes_majority(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    data = SCData(np.ones((10, 4)), pd.DataFrame({"cell_type": ["Neurons"] * 10}), pd.DataFrame(index=range(4)))
    uri = "adata://s/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "s", "input", "input")
    task = TaskContract(task_id="subset", capability="subset_cells", input_artifacts=[uri], parameters={"cell_type": "Microglia"})
    with pytest.raises(ValueError, match="refusing to substitute"):
        SubsetCapability().execute(task, reg)
    assert len(reg.list_artifacts()) == 1


def test_metadata_audit_prunes_unsupported_statistics_and_consumers(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    data = SCData(np.random.default_rng(8).poisson(4, (60, 30)),
                  pd.DataFrame(index=range(60)), pd.DataFrame({"gene_name": [f"g{i}" for i in range(30)]}))
    uri = "adata://s/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "s", "input", "input")
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="human", tissue="kidney"), data=DataSpec(raw_artifact_uri=uri))
    orch = ScientificOrchestrator(artifact_registry=reg)
    result = orch.run_study(manifest, {"include_knowledge": True, "method_profile": "baseline"})
    assert result["status"] == "success", result["failures"]
    omitted = {item["task_id"] for item in result["planning_decisions"]}
    assert {"task_007_abundance", "task_008_deg", "task_016_knowledge"} <= omitted
    assert not any(t.capability == "subset_cells" for t in orch.task_history)
    assert not any("Microglia" in c["statement"] or "Alzheimer" in c["statement"] for c in result["claims"])
    resumed = orch.run_study(manifest, {"include_knowledge": True, "method_profile": "baseline", "resume": True})
    assert resumed["status"] == "success", resumed["failures"]
    assert resumed["planning_decisions"] == result["planning_decisions"]


def test_invalid_explicit_contrast_is_not_silently_omitted():
    tasks = [TaskContract(task_id="deg", capability="deg", parameters={"condition_a": "missing", "condition_b": "control"})]
    with pytest.raises(ValueError, match="Explicit contrast"):
        ComputationalDAGPlanner.adapt_after_audit(tasks, {"observed_conditions": ["control"], "condition_metadata_complete": True})
