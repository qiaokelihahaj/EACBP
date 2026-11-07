"""Sensitivity orchestration preserves native evidence and failed branches."""
from copy import deepcopy
import json

import numpy as np
import pandas as pd
import pytest

from eacbp.application import autonomous, study_service
from eacbp.application import qc_sensitivity_execution as execution
from eacbp.application.qc_sensitivity import freeze_qc_sensitivity_plan
from eacbp.capabilities.sc_data import SCData
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.schemas.study import StudyManifest, BiologicalDesign


def test_native_plan_fit_receipts_are_unique_and_declared_after_target_expansion():
    manifest = StudyManifest(study_id="fit_plan", biological_design=BiologicalDesign(
        species="human", tissue="blood", conditions=["ctrl", "stim"], target_cell_types=["T", "B"]))
    tasks = ComputationalDAGPlanner.build_study_plan(manifest, {"advanced_analysis": True})
    fits = [task for task in tasks if task.capability in {"deg", "donor_sensitivity"}]
    assert len(fits) == 4
    assert len({task.parameters["fit_receipts_uri"] for task in fits}) == 4
    for task in fits:
        assert task.expected_outputs[-1] == task.parameters["fit_receipts_uri"]
        assert task.parameters["target_branch"] in task.expected_outputs[-1]
        assert task.expected_outputs[0].startswith("table://")
    assert all("fit_receipts_uri" not in task.parameters for task in tasks if task.capability == "functional_activity")


@pytest.mark.parametrize("failed", [False, True])
def test_sensitivity_keeps_primary_files_and_records_failed_branch(tmp_path, monkeypatch, failed):
    payload = SCData(np.array([[1] * count + [0] * (20 - count) for count in (1, 3, 5, 7, 8, 9, 10, 10, 11, 12, 15, 20)]),
                     pd.DataFrame(index=[f"c{i}" for i in range(12)]), pd.DataFrame(index=[f"g{i}" for i in range(20)]))
    source = tmp_path / "raw.h5ad"
    source.write_bytes(b"fixed input identity")
    digest = study_service._hash_file(source)
    observation = autonomous.observe_input(payload, question="Explore populations")
    config = {"autonomous_analysis": {"question": "Explore populations", "source_sha256": digest},
              "capability_parameters": {"qc": {key: observation["qc"][key] for key in ("min_genes", "max_mito_pct")}},
              "method_overrides": {"clustering": "unchanged"}}
    plan = freeze_qc_sensitivity_plan(payload=payload, observations=observation, source_sha256=digest, baseline_config=config)
    config["autonomous_analysis"]["qc_sensitivity_plan"] = plan
    run_dir = tmp_path / "primary"
    run_dir.mkdir()
    for name in ("snapshot.json", "summary.json", "report.md", "config.json"):
        (run_dir / name).write_text("original " + name, encoding="utf-8")
    before = {path.name: path.read_bytes() for path in run_dir.iterdir()}
    original_config = deepcopy(config)
    calls = []

    def run(**kwargs):
        calls.append(kwargs)
        if failed and kwargs["run_dir"].name == "stricter_2mad":
            raise RuntimeError("a genuine downstream failure")
        return {"status": "success"}

    def collect(path, **kwargs):
        return {"qc_admitted": True, "provenance": {"report": str(path / "report.md"), "snapshot_sha256": "a" * 64}}

    def compare(frozen, *, baseline, branches, observed_groups):
        assert frozen == plan
        return {"status": "partial" if len(branches) < 2 else "completed", "plan_id": frozen["plan_id"],
                "interpretation": frozen["interpretation"], "branches": []}

    monkeypatch.setattr(study_service, "run_study", run)
    monkeypatch.setattr(execution, "collect_qc_run", collect)
    monkeypatch.setattr(execution, "compare_qc_sensitivity_runs", compare)
    manifest = StudyManifest(study_id="original", biological_design=BiologicalDesign(species="human", tissue="blood"))
    summary = {"status": "success", "analysis_request_status": "partial"}
    result = execution.execute_qc_sensitivity(plan=plan, manifest=manifest, config=config, data=source,
                                               run_dir=run_dir, baseline_summary=summary)
    assert result["status"] == "success" and result["analysis_request_status"] == "partial"
    assert result["qc_sensitivity"]["status"] == ("partial" if failed else "completed")
    assert config == original_config and manifest.study_id == "original"
    assert all((run_dir / name).read_bytes() == content for name, content in before.items())
    assert len(calls) == 2
    for call in calls:
        assert call["config"]["method_overrides"] == config["method_overrides"]
        assert call["data"] == source
        assert call["manifest"].study_id.startswith("original_qc_")
    report = json.loads((run_dir / "qc_sensitivity" / "comparison.json").read_text(encoding="utf-8"))
    assert bool(report["failures"]) is failed
    if failed:
        assert "genuine downstream failure" in report["failures"][0]["reason"]
