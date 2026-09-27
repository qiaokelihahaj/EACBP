"""Report summaries must use the threshold of the audited analysis."""

from types import SimpleNamespace

import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.evidence.graph import EvidenceGraph
from eacbp.report.markdown_report import ScientificReportGenerator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import TaskStatus


@pytest.mark.parametrize("alpha, expected", [(0.1, 2), (0.01, 0), (0.05, 1)])
def test_report_uses_recorded_alpha_and_effect_definition(tmp_path, alpha, expected):
    registry = ArtifactRegistry(str(tmp_path / "artifacts"))
    uri = "table://report_inference/deg/v1"
    table = pd.DataFrame({
        "gene": ["g1", "g2", "g3"],
        "fdr_q_value": [0.02, 0.08, float("nan")],
        "effect_definition": ["control versus treated"] * 3,
        "alpha": [alpha] * 3,
    })
    registry.register(uri, table, ArtifactType.TABLE, "report_inference", "deg", "test")
    task = SimpleNamespace(
        task_id="deg", capability="deg", method_used="pydeseq2_pseudobulk_v1",
        status=TaskStatus.SUCCESS, metrics={"alpha": alpha}, output_artifacts=[uri],
    )
    manifest = StudyManifest(
        study_id="report_inference",
        biological_design=BiologicalDesign(species="human", tissue="test"),
        data=DataSpec(raw_artifact_uri="adata://report_inference/raw/v1"),
    )
    audit = SimpleNamespace(target_task_id="deg", overall_passed=True, stop_rule_triggered=False)
    report = ScientificReportGenerator(manifest, EvidenceGraph(), registry, [task], [audit])
    text = "\n".join(report._advanced_analysis_details())
    assert f"FDR < {alpha:g}: {expected}; unestimated FDR: 1." in text
    assert "control versus treated" in text
