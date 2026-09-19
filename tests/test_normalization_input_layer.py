import numpy as np
import pandas as pd
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.normalization import NormalizationCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract


def test_corrected_layer_normalization_preserves_original_counts(tmp_path):
    raw = np.array([[100, 10], [50, 50]])
    corrected = np.array([[20, 10], [40, 40]])
    data = SCData(raw, pd.DataFrame(index=["a", "b"]), pd.DataFrame(index=["x", "y"]),
                  layers={"counts": raw, "corrected_counts": corrected})
    reg = ArtifactRegistry(str(tmp_path))
    uri = "adata://s/raw/v1"
    reg.register(uri, data, ArtifactType.ANNDATA, "s", "import", "import")
    task = TaskContract(task_id="norm", capability="normalization", input_artifacts=[uri],
                        parameters={"input_layer": "corrected_counts", "target_sum": 100})
    result = NormalizationCapability().execute(task, reg)
    _, out = reg.get(result.output_artifacts[0])
    if isinstance(out, dict):
        out = SCData.from_dict(out)
    np.testing.assert_array_equal(out.layers["counts"], raw)
    np.testing.assert_array_equal(out.layers["corrected_counts"], corrected)
    np.testing.assert_allclose(out.X, np.log1p(corrected / corrected.sum(axis=1)[:, None] * 100), rtol=1e-6)
    assert out.uns["normalization"]["input_layer"] == "corrected_counts"
