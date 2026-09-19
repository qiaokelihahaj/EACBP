from pathlib import Path

import pytest

from eacbp.orchestrator.advanced_plan import extend_plan
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import TaskContract


def _manifest(study_id="background_plan"):
    return StudyManifest(
        study_id=study_id,
        biological_design=BiologicalDesign(species="human", tissue="test"),
        data=DataSpec(raw_artifact_uri=f"adata://{study_id}/raw/v1"),
    )


def _base_tasks(study_id="background_plan", normalization_parameters=None):
    raw = f"adata://{study_id}/raw/v1"
    qc = f"adata://{study_id}/qc/v1"
    normalized = f"adata://{study_id}/normalized/v2"
    return [
        TaskContract(
            task_id="task_001_audit",
            capability="dataset_audit",
            input_artifacts=[raw],
            expected_outputs=[f"table://{study_id}/dataset_audit/v1"],
        ),
        TaskContract(
            task_id="task_002_qc",
            capability="qc",
            input_artifacts=[raw],
            expected_outputs=[qc],
        ),
        TaskContract(
            task_id="task_003_norm",
            capability="normalization",
            input_artifacts=[qc],
            parameters=dict(normalization_parameters or {}),
            expected_outputs=[normalized],
        ),
        TaskContract(
            task_id="deg",
            capability="deg",
            input_artifacts=[normalized],
        ),
    ]


def _cellbender_parameters(tmp_path):
    raw_input = tmp_path / "unfiltered.h5"
    executable = tmp_path / "cellbender"
    output = tmp_path / "corrected.h5"
    raw_input.write_bytes(b"unfiltered droplet input")
    executable.write_bytes(b"cellbender executable")
    return {
        "unfiltered_input_path": str(raw_input),
        "executable": str(executable),
        "output_path": str(output),
        "run_cwd": str(tmp_path),
    }


def test_background_removal_is_available_and_precedes_qc(tmp_path):
    study_id = "background_plan"
    manifest = _manifest(study_id)
    state = {
        "analysis_extensions": {"background_removal": _cellbender_parameters(tmp_path)},
    }
    planned = extend_plan(_base_tasks(study_id), manifest, state)

    background = next(task for task in planned if task.capability == "background_removal")
    qc = next(task for task in planned if task.capability == "qc")
    normalization = next(task for task in planned if task.capability == "normalization")
    corrected = f"adata://{study_id}/background_corrected/v1"

    assert background.input_artifacts == [f"adata://{study_id}/raw/v1"]
    assert background.expected_outputs == [
        corrected,
        f"json://{study_id}/cellbender_report/v1",
    ]
    assert qc.input_artifacts == [corrected]
    assert normalization.parameters["input_layer"] == "corrected_counts"


def test_background_dag_connects_external_output_to_qc_and_report(tmp_path):
    study_id = "background_dag"
    manifest = _manifest(study_id)
    state = {
        "analysis_extensions": {"background_removal": _cellbender_parameters(tmp_path)},
    }
    planned = ComputationalDAGPlanner.build_study_plan(manifest, state)
    background = next(task for task in planned if task.capability == "background_removal")
    qc = next(task for task in planned if task.capability == "qc")
    normalization = next(task for task in planned if task.capability == "normalization")

    assert background.expected_outputs[1] == f"json://{study_id}/cellbender_report/v1"
    assert background.task_id in qc.depends_on
    assert qc.expected_outputs[0] in normalization.input_artifacts
    ordered = ComputationalDAGPlanner.order_tasks(planned)
    assert [task.task_id for task in ordered].index(background.task_id) < [task.task_id for task in ordered].index(qc.task_id)


def test_background_removal_missing_paths_fail_during_planning(tmp_path):
    with pytest.raises(ValueError, match="unfiltered_input_path.*executable.*output_path"):
        extend_plan(
            _base_tasks(),
            _manifest(),
            {"analysis_extensions": {"background_removal": True}},
        )


def test_background_removal_rejects_conflicting_normalization_layer(tmp_path):
    parameters = _cellbender_parameters(tmp_path)
    with pytest.raises(ValueError, match="corrected_counts"):
        extend_plan(
            _base_tasks(normalization_parameters={"input_layer": "counts"}),
            _manifest(),
            {"analysis_extensions": {"background_removal": parameters}},
        )


def test_background_removal_rejects_invalid_output_path(tmp_path):
    parameters = _cellbender_parameters(tmp_path)
    parameters["output_path"] = str(tmp_path / "corrected.h5ad")
    with pytest.raises(ValueError, match=r"output_path.*\.h5"):
        extend_plan(
            _base_tasks(),
            _manifest(),
            {"analysis_extensions": {"background_removal": parameters}},
        )
