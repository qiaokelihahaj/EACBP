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


def test_report_distinguishes_question_completion_from_successful_tasks(tmp_path):
    manifest = StudyManifest(study_id="question", biological_design=BiologicalDesign(species="human", tissue="test"))
    report = ScientificReportGenerator(
        manifest, EvidenceGraph(), ArtifactRegistry(tmp_path / "artifacts"), [],
        run_summary={"status": "success", "analysis_request_status": "partial",
                     "autonomy_kind": "deterministic_rule_driven", "policy_version": "observed_single_cell_v1",
                     "scientific_gaps": ["Missing donor contrast"],
                     "analysis_limitations": ["Direction is model inferred"]},
    )
    text = report.generate_markdown()
    assert "Requested analysis: `partial`" in text
    assert "Missing donor contrast" in text and "Direction is model inferred" in text
    assert "does not certify agreement with an external reproduction benchmark" in text


def test_legacy_report_has_no_invented_question_completion(tmp_path):
    manifest = StudyManifest(study_id="legacy", biological_design=BiologicalDesign(species="human", tissue="test"))
    report = ScientificReportGenerator(manifest, EvidenceGraph(), ArtifactRegistry(tmp_path / "artifacts"), [])
    assert "Scientific question completion" not in report.generate_markdown()


def test_qc_report_only_explains_admitted_dispositions(tmp_path):
    manifest = StudyManifest(study_id="qc_report", biological_design=BiologicalDesign(species="human", tissue="test"))
    task = SimpleNamespace(
        task_id="qc", capability="qc", status=TaskStatus.SUCCESS,
        metrics={"cell_disposition_schema_version": 1,
                 "cell_disposition_counts": {"retained": 8, "low_gene_count": 2},
                 "interpretation_warning": "Low RNA alone is not independent evidence of damage."},
        output_artifacts=["adata://qc_report/qc/v1"],
    )
    report = ScientificReportGenerator(manifest, EvidenceGraph(), ArtifactRegistry(tmp_path / "artifacts"), [task])
    assert report._qc_disposition_details() == []
    report.audit_reports = [SimpleNamespace(target_task_id="qc", overall_passed=True, stop_rule_triggered=False)]
    text = "\n".join(report._qc_disposition_details())
    assert "| low_gene_count | 2 |" in text
    assert "adata://qc_report/qc/v1" in text
    assert "not independent evidence of damage" in text
