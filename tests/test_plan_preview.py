from copy import deepcopy
from unittest.mock import patch

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.planning import build_study_tasks, preview_study_plan, resolve_task_contract
from eacbp.schemas.study import BiologicalDesign, DataSpec, StudyManifest
from eacbp.schemas.task import TaskContract


def manifest():
    return StudyManifest(study_id='preview', biological_design=BiologicalDesign(
        species='human', tissue='brain', target_cell_types=['Microglia', 'Neurons']),
        data=DataSpec(raw_artifact_uri='adata://preview/raw/v1'))


def test_preview_resolves_branches_without_storage_or_execution(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    config = {'method_profile': 'standard', 'target_parameters': {
        'Microglia': {'trajectory_inference': {'root_cell_id': 'm1'}}}}
    before = deepcopy(config)
    with patch.object(ArtifactRegistry, '__init__', side_effect=AssertionError('must not open storage')), \
         patch.object(CapabilityRegistry, 'execute_contract', side_effect=AssertionError('must not compute')):
        preview = preview_study_plan(manifest(), config)
    assert preview['valid'], preview
    assert config == before
    assert list(tmp_path.iterdir()) == []
    roots = [task for task in preview['tasks'] if task['capability'] == 'trajectory_inference']
    assert len(roots) == 2
    assert all(task['method'] == 'scanpy_dpt_v1' for task in roots)
    assert preview['missing_parameters'] == [{
        'task_id': 'task_009_trajectory_neurons', 'parameter': 'root_cell_id',
        'effect': 'trajectory_and_consumers_omitted_after_dataset_audit'}]
    assert preview['external_inputs'] == ['adata://preview/raw/v1']


def test_preview_matches_runtime_contract_resolution():
    study = manifest()
    state = {'method_profile': 'baseline', 'method_overrides': {'integration': 'no_correction_v1'}}
    registry = create_default_capability_registry()
    tasks = build_study_tasks(study, state, registry)
    for task in tasks:
        resolve_task_contract(task, study, state, registry)
    preview = preview_study_plan(study, state, registry)
    assert preview['valid']
    assert [{key: value for key, value in task.items() if key not in {'target_branch', 'target_cell_type'}}
            for task in preview['tasks']] == [task.model_dump(mode='json') for task in tasks]


def test_preview_reports_static_error_without_executing():
    preview = preview_study_plan(manifest(), {'method_overrides': {'deg': 'not_registered'}})
    assert not preview['valid']
    assert {error['task_id'] for error in preview['errors']} == {
        'task_008_deg_microglia', 'task_008_deg_neurons'}
    preview = preview_study_plan(manifest(), {'analysis_extensions': {'background_removal': True}})
    assert not preview['valid']
    assert 'unfiltered_input_path' in preview['errors'][0]['message']


def test_preview_rejects_colliding_outputs():
    tasks = [TaskContract(task_id=name, capability='dataset_audit', expected_outputs=['json://preview/same/v1'])
             for name in ['a', 'b']]
    with patch.object(ComputationalDAGPlanner, 'build_study_plan', return_value=tasks):
        preview = preview_study_plan(manifest())
    assert not preview['valid']
    assert 'same artifact URI' in preview['errors'][0]['message']


def test_preview_marks_missing_fastq_resources_as_blocking():
    study = manifest()
    study.data.has_raw_fastq = True
    study.data.raw_artifact_uri = 'fastq://preview/raw/v1'
    preview = preview_study_plan(study, {'quant_tool': 'kb_python_v1'})
    assert not preview['valid']
    assert {'index_path', 't2g_path'} <= {item['parameter'] for item in preview['missing_parameters']}
