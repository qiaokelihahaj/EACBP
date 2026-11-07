"""Marker representation must respect the study's explicit reference."""

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.biological import BiologicalValidator
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


@pytest.mark.parametrize(
    "label,genes,parameters,expected",
    [
        ("B cells", ["MS4A1", "CD79A"], {"marker_reference": {"B cells": ["MS4A1"]}}, True),
        ("B cells", ["OTHER"], {"marker_reference": {"B cells": ["MS4A1"]}}, False),
        ("B cells", ["MS4A1"], {"marker_reference": {"B cells": "MS4A1"}}, False),
        ("B cells", ["MS4A1"], {"marker_reference": {}}, False),
        ("Microglia", ["Cx3cr1"], {}, True),
        ("Microglia", ["OTHER"], {}, False),
    ],
)
def test_marker_reference_and_missing_markers(tmp_path, label, genes, parameters, expected):
    registry = ArtifactRegistry(tmp_path / "artifacts")
    data = SCData(
        X=np.ones((3, len(genes))),
        obs=pd.DataFrame({"cell_type": [label] * 3}),
        var=pd.DataFrame(index=pd.Index(genes)),
    )
    uri = "adata://public_marker_test/clustered/v1"
    registry.register(uri, data, ArtifactType.ANNDATA, "public_marker_test", "cluster", "test")
    contract = TaskContract(task_id="cluster", capability="clustering", parameters=parameters)
    result = TaskResult(
        task_id="cluster", capability="clustering", status=TaskStatus.SUCCESS,
        method_used="test", output_artifacts=[uri],
    )
    check = BiologicalValidator().audit(contract, result, registry).checks[0]
    assert check.passed is expected
    assert "annotation accuracy" in check.message
    if expected:
        assert check.metrics["matched_markers"][label]
