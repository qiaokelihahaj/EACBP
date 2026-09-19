"""End-to-end coverage for the legacy standard-study entry point."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from eacbp import cli
from eacbp.orchestrator.execution import TaskExecutor
from eacbp.orchestrator.dag import ComputationalDAGPlanner


def _source(path: Path) -> Path:
    data = ad.AnnData(
        np.random.default_rng(17).poisson(4, (16, 6)).astype(float),
        obs=pd.DataFrame({"condition": ["case"] * 8 + ["control"] * 8}, index=[f"c{i}" for i in range(16)]),
        var=pd.DataFrame(index=[f"g{i}" for i in range(6)]),
    )
    data.write_h5ad(path)
    return path


@pytest.fixture
def audit_only(monkeypatch):
    original = ComputationalDAGPlanner.build_study_plan

    def plan(manifest, config):
        return [task for task in original(manifest, config) if task.capability == "dataset_audit"]

    monkeypatch.setattr(ComputationalDAGPlanner, "build_study_plan", plan)


def _script_module():
    path = Path(__file__).parents[1] / "scripts" / "run_standard_study.py"
    spec = importlib.util.spec_from_file_location("run_standard_study_under_test", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_script_resume_reuses_snapshot(tmp_path, capsys, audit_only):
    source = _source(tmp_path / "source.h5ad")
    output_dir = tmp_path / "r"
    script = _script_module()

    assert script.main([
        "--data", str(source),
        "--species", "human",
        "--tissue", "brain",
        "--study-id", "s",
        "--output-dir", str(output_dir),
    ]) == 0
    report_line = next(line for line in capsys.readouterr().out.splitlines() if line.startswith("Report: "))
    report = Path(report_line.removeprefix("Report: "))
    run_dir = report.parent
    assert report.is_file()
    assert (run_dir / "snapshot.json").is_file()
    saved = json.loads((run_dir / "run_config.json").read_text(encoding="utf-8"))
    assert saved["import_completed"] is True

    source.write_bytes(b"changed after import; resume must use the registry")
    with patch.object(cli, "_import_h5ad", side_effect=AssertionError("resume must not re-import")), \
         patch.object(TaskExecutor, "execute", side_effect=AssertionError("resume must reuse checkpoint")):
        resumed = cli.resume_study(run_dir=run_dir)
    assert resumed["status"] == "success", resumed
    assert len(list((run_dir / "snapshots").glob("*.json"))) == 2
