"""Focused coverage for typed runtime state and orchestration boundaries."""

from unittest.mock import patch

import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.auditor.base import ValidationReport
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.runtime import ExecutionState, RunConfig
from eacbp.schemas.study import BiologicalDesign, StudyManifest
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def test_execution_state_keeps_configuration_separate_from_observations():
    state = ExecutionState.from_input({
        "method_profile": "baseline",
        "mode": "demo",
        "n_cells": 10,
        "resume": True,
    })
    state.record_metrics({
        "method_profile": "standard",
        "mode": "real",
        "method_overrides": {"deg": "unsafe"},
        "analysis_extensions": {"malicious": True},
        "n_cells": 2,
        "observed_conditions": ["control", "treated"],
    })

    assert state.resume_requested is True
    assert state.current_state["method_profile"] == "baseline"
    assert state.current_state["mode"] == "demo"
    assert state.current_state["n_cells"] == 10
    assert "method_overrides" not in state.current_state
    assert "analysis_extensions" not in state.current_state
    assert state.current_state["observed_conditions"] == ["control", "treated"]


def test_run_config_rejects_coerced_control_values():
    with pytest.raises(ValueError):
        RunConfig.from_mapping({"method_profile": "baseline", "include_knowledge": "false"})


@pytest.mark.parametrize("value", ["false", "true", 0, 1])
def test_legacy_resume_rejects_non_boolean_controls(value):
    with pytest.raises(ValueError):
        ExecutionState.from_input({"resume": value})


def test_target_observations_do_not_contaminate_siblings():
    state = ExecutionState.from_input({"method_profile": "baseline"})
    state.record_metrics({"n_cells": 100, "batches": ["shared"]})
    state.record_metrics({"n_cells": 12, "method_overrides": {"deg": "rogue"}}, branch="T")
    state.record_metrics({"n_cells": 24}, branch="B")
    assert state.current_state["n_cells"] == 100
    assert state.state_for_branch("T")["n_cells"] == 12
    assert state.state_for_branch("B")["n_cells"] == 24
    assert state.state_for_branch("T")["batches"] == ["shared"]
    assert "method_overrides" not in state.state_for_branch("T")


def test_orchestrator_does_not_let_result_metrics_replace_user_options(tmp_path):
    registry = ArtifactRegistry(str(tmp_path))
    manifest = StudyManifest(
        study_id="runtime_state",
        biological_design=BiologicalDesign(species="human", tissue="test"),
    )
    uri = "table://runtime_state/audit/v1"
    task = TaskContract(
        task_id="audit",
        capability="dataset_audit",
        expected_outputs=[uri],
    )

    def execute(contract, staged):
        staged.register(uri, pd.DataFrame({"n_cells": [4]}), ArtifactType.TABLE,
                        "runtime_state", "audit", "audit")
        return TaskResult(
            task_id=contract.task_id,
            capability=contract.capability,
            method_used="sc_audit_v1",
            status=TaskStatus.SUCCESS,
            output_artifacts=[uri],
            metrics={
                "method_profile": "standard",
                "mode": "real",
                "n_cells": 4,
                "method_overrides": {"deg": "unsafe"},
            },
        )

    report = ValidationReport(auditor_name="test", target_task_id="audit")
    orchestrator = ScientificOrchestrator(artifact_registry=registry)
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[task]), \
         patch.object(orchestrator.capability_registry, "execute_contract", side_effect=execute), \
         patch.object(orchestrator.auditor, "audit_task", return_value=report):
        summary = orchestrator.run_study(
            manifest,
            {"method_profile": "baseline", "mode": "demo", "n_cells": 9},
        )

    assert summary["status"] == "success", summary["failures"]
    assert orchestrator.current_state["method_profile"] == "baseline"
    assert orchestrator.current_state["mode"] == "demo"
    assert orchestrator.current_state["n_cells"] == 9
    assert "method_overrides" not in orchestrator.current_state
