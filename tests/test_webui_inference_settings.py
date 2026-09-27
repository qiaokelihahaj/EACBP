from types import SimpleNamespace

import pandas as pd

from eacbp.webui.results import _deg_view


def test_deg_preview_preserves_actual_contrast_and_threshold():
    metadata = SimpleNamespace(
        uri="table://inference/deg/v1",
        parameters={"condition_a": "treated", "condition_b": "control"},
        summary_metrics={},
    )
    task = SimpleNamespace(
        task_id="deg", method_used="pydeseq2_pseudobulk_v1",
        metrics={"condition_a": "control", "condition_b": "treated",
                 "contrast_label": "control versus treated", "alpha": 0.1,
                 "effect_definition": "PyDESeq2 log2 fold change for control versus treated"},
    )
    table = pd.DataFrame({"gene": ["g"], "log2_fold_change": [-1.0], "fdr_q_value": [0.08]})
    view = _deg_view(metadata, table, task)
    assert view["condition_a"] == "control"
    assert view["contrast_label"] == "control versus treated"
    assert view["effect_definition"] == task.metrics["effect_definition"]
    assert view["alpha"] == 0.1
