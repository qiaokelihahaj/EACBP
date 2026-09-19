import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix
pytest.importorskip("cellrank")
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.real_methods import CellRankFateCapability
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def test_actual_cellrank_matches_known_absorption_probabilities(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    transition = csr_matrix([[0, 0, .25, .75], [0, 0, .6, .4], [0, 0, 1, 0], [0, 0, 0, 1]])
    data = SCData(np.ones((4, 3)), pd.DataFrame({"cell_id": ["a", "b", "left", "right"]}),
                  pd.DataFrame(index=range(3)), obsp={"transition_matrix": transition})
    uri = "adata://fate/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "fate", "input", "input")
    task = TaskContract(task_id="f", capability="fate_mapping", input_artifacts=[uri],
                        parameters={"kernel": "precomputed", "terminal_states": {"L": ["left"], "R": ["right"]}})
    result = CellRankFateCapability().execute(task, reg)
    metadata, table = reg.get(result.output_artifacts[0])
    assert np.allclose(table.loc["a", ["L", "R"]], [.25, .75])
    assert np.allclose(table.loc["b", ["L", "R"]], [.6, .4])
    assert "cellrank" in metadata.software_versions
    assert result.method_used == "cellrank_fate_v1"
