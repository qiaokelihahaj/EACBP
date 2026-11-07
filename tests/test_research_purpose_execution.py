"""Real application lifecycle for purpose-driven runs (standard CI lane)."""
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from eacbp.application.study_service import run_study, resume_study, report_study
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.webui.service import WebService


def input_file(tmp_path):
    obs = pd.DataFrame({"condition": ["treated"] * 40 + ["control"] * 40,
                        "donor_id": [f"d{i // 10}" for i in range(80)]},
                       index=[f"cell{i}" for i in range(80)])
    data = ad.AnnData(np.random.default_rng(25).poisson(4, (80, 40)).astype(np.int64),
                     obs=obs, var=pd.DataFrame(index=[f"g{i}" for i in range(40)]))
    path = tmp_path / "synthetic.h5ad"
    data.write_h5ad(path)
    return path


def test_overview_does_not_require_a_comparison_direction(tmp_path):
    path = input_file(tmp_path)
    service = WebService(tmp_path, tmp_path / "runs")
    preview = service.preview({"data": str(path), "species": "human", "tissue": "synthetic", "research_purpose": "overview",
                               "method_profile": "baseline", "min_genes": 1, "condition_col": "condition", "donor_col": "donor_id"})
    assert preview["can_run"], preview["design"]
    assert not {"deg", "differential_abundance", "trajectory_inference"}.intersection(t["capability"] for t in preview["plan"]["tasks"])
    preview = service.preview({"data": str(path), "species": "human", "tissue": "synthetic", "research_purpose": "overview",
                               "method_profile": "baseline", "min_genes": 1, "condition_col": "condition", "donor_col": "donor_id",
                               "config_overrides": {"capability_parameters": {"deg": {"donor_col": "unused_missing_column"}}}})
    assert preview["can_run"], preview["design"]


def test_purpose_completion_persists_through_resume_and_report(tmp_path):
    path = input_file(tmp_path)
    manifest = StudyManifest(study_id="purpose", biological_design=BiologicalDesign(species="human", tissue="synthetic"))
    run_dir = tmp_path / "run"
    config = {"research_purpose": "overview", "method_profile": "baseline", "capability_parameters": {"qc": {"min_genes": 1}}}
    summary = run_study(manifest=manifest, config=config, data=path, run_dir=run_dir)
    assert summary["status"] == "success", summary
    assert summary["analysis_request_status"] == "complete", summary
    saved = json.loads((run_dir / "run_config.json").read_text())
    assert saved["config"]["analysis_request"]["requirement_contract"]["source_sha256"] == saved["source"]["sha256"]
    assert json.loads((run_dir / "config.json").read_text()) == saved["config"]
    resumed = resume_study(run_dir=run_dir)
    assert resumed["requirement_contract_id"] == summary["requirement_contract_id"]
    assert resumed["analysis_request_status"] == "complete"
    report_study(run_dir=run_dir)
    report = (run_dir / "report.md").read_text(encoding="utf-8")
    assert "Requested analysis: `complete`" in report
    reuse = WebService(tmp_path, tmp_path).reuse_configuration("run")
    assert reuse["form"]["research_purpose"] == "overview"
    assert reuse["research_purpose_saved"] is True
    assert "analysis_request" not in reuse["form"]["config_overrides"]
    assert "analysis_request" in reuse["original"]["config"]


@pytest.mark.parametrize("purpose", ["deg", "abundance"])
def test_baseline_inference_can_complete_without_positive_findings(tmp_path, purpose):
    path = input_file(tmp_path)
    data = ad.read_h5ad(path)
    data.X[40:] = data.X[:40]
    data.write_h5ad(path)
    manifest = StudyManifest(study_id="baseline", biological_design=BiologicalDesign(species="human", tissue="synthetic"))
    summary = run_study(manifest=manifest,
                        config={"research_purpose": purpose, "method_profile": "baseline",
                                "capability_parameters": {"qc": {"min_genes": 1}}},
                        data=path, run_dir=tmp_path / "run")
    assert summary["status"] == "success", summary
    assert summary["analysis_request_status"] == "complete", summary
    scientific = next(check for check in summary["analysis_completion_checks"]
                      if check["dimension"] == ("deg" if purpose == "deg" else "differential_abundance"))
    assert scientific["status"] == "passed"
    assert scientific["scientific_statuses"] == ["estimated_inconclusive"]
