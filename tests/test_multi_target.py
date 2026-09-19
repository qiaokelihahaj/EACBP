from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.orchestrator.dag import ComputationalDAGPlanner, target_branch_slug_map
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest, BiologicalDesign


def manifest(targets):
    return StudyManifest(study_id='s', biological_design=BiologicalDesign(
        species='human', tissue='brain', target_cell_types=targets))


def test_plan_shares_preprocessing_and_namespaces_all_consumers():
    study = manifest(['Microglia', 'Neurons'])
    tasks = ComputationalDAGPlanner.build_study_plan(study, {
        'include_knowledge': True, 'analysis_extensions': {'donor_sensitivity': True}})
    assert sum(t.capability == 'qc' for t in tasks) == 1
    outputs = [uri for t in tasks for uri in t.expected_outputs]
    assert len(outputs) == len(set(outputs))
    by_id = {t.task_id: t for t in tasks}
    for capability in ('subset_cells', 'deg', 'trajectory_inference', 'donor_sensitivity', 'knowledge_retrieval'):
        assert sum(t.capability == capability for t in tasks) == 2
    for task in tasks:
        branch = task.parameters.get('target_branch')
        if not branch:
            continue
        if task.capability == 'subset_cells':
            assert task.parameters['cell_type'] == task.parameters['target_cell_type']
        for parent in task.depends_on:
            assert by_id[parent].parameters.get('target_branch') in (None, branch)
    ComputationalDAGPlanner.order_tasks(tasks)


def test_slug_collisions_and_per_target_roots():
    targets = ['T cell', 'T-cell', '星形胶质细胞']
    slugs = target_branch_slug_map(targets)
    assert len(set(slugs.values())) == 3
    assert slugs == target_branch_slug_map(list(reversed(targets)))
    with pytest.raises(ValueError, match='Duplicate'):
        target_branch_slug_map(['T', 'T'])
    tasks = ComputationalDAGPlanner.build_study_plan(manifest(['Microglia', 'Neurons']), {
        'method_profile': 'standard', 'capability_parameters': {'trajectory_inference': {'root_cell_id': 'global'}},
        'target_parameters': {'Microglia': {'trajectory_inference': {'root_cell_id': 'microglia_root'}}}})
    roots = {t.parameters['target_cell_type']: t.parameters.get('root_cell_id')
             for t in tasks if t.capability == 'trajectory_inference'}
    assert roots == {'Microglia': 'microglia_root', 'Neurons': None}


def test_cellrank_requires_separate_terminal_states_for_every_branch():
    study = manifest(['Microglia', 'Neurons'])
    state = {'method_profile': 'standard', 'target_parameters': {
        target: {'trajectory_inference': {'root_cell_id': target + '_root'},
                 'fate_mapping': {'terminal_states': {'terminal': [target + '_terminal']}}}
        for target in study.biological_design.target_cell_types}}
    tasks = ComputationalDAGPlanner.build_study_plan(study, state)
    fates = [task for task in tasks if task.capability == 'fate_mapping']
    assert len(fates) == 2
    for task in fates:
        target = task.parameters['target_cell_type']
        assert task.parameters['terminal_states'] == {'terminal': [target + '_terminal']}
        assert task.input_artifacts[0].endswith('/' + task.parameters['target_branch'] + '/v1')
    del state['target_parameters']['Neurons']['fate_mapping']
    with pytest.raises(ValueError, match='terminal_states'):
        ComputationalDAGPlanner.build_study_plan(study, state)


def test_subset_override_cannot_mislabel_a_branch():
    with pytest.raises(ValueError, match='match its target branch'):
        ComputationalDAGPlanner.build_study_plan(manifest(['Microglia', 'Neurons']), {
            'target_parameters': {'Neurons': {'subset_cells': {'cell_type': 'Microglia'}}}})


@pytest.mark.parametrize('missing_target', [False, True])
def test_real_branch_execution_resume_and_failure_isolation(tmp_path, missing_target):
    study = manifest(['Microglia', 'Absent' if missing_target else 'Neurons'])
    tasks = [t for t in ComputationalDAGPlanner.build_study_plan(study, {})
             if t.capability in {'subset_cells', 'deg'}]
    task_ids = {t.task_id for t in tasks}
    for task in tasks:
        task.depends_on = [parent for parent in task.depends_on if parent in task_ids]
    obs = pd.DataFrame([
        {'cell_type': target, 'condition': condition, 'donor': f'{condition}_{donor}'}
        for target in ['Microglia', 'Neurons'] for condition in ['case', 'control']
        for donor in range(4) for cell in range(8)
    ], index=[f'cell_{i}' for i in range(128)])
    counts = np.random.default_rng(6).poisson(5, (len(obs), 25)).astype(float)
    counts[obs['condition'].to_numpy() == 'case', :5] *= 5
    data = SCData(np.log1p(counts), obs, pd.DataFrame({'gene_name': [f'g{i}' for i in range(25)]}),
                  layers={'counts': counts})
    registry = ArtifactRegistry(str(tmp_path))
    registry.register('adata://s/annotated/v4', data, ArtifactType.ANNDATA, 's', 'input', 'fixture')
    orchestrator = ScientificOrchestrator(registry)
    with patch.object(ComputationalDAGPlanner, 'build_study_plan', side_effect=lambda *_: [t.model_copy(deep=True) for t in tasks]):
        first = orchestrator.run_study(study)
        assert first['status'] == ('failed' if missing_target else 'success'), first
        nodes = list(orchestrator.evidence_graph.evidence_nodes.values())
        assert nodes
        assert {node.biological_context['target_cell_type'] for node in nodes} == (
            {'Microglia'} if missing_target else {'Microglia', 'Neurons'})
        successful = [t for t in orchestrator.task_history if t.status.value == 'success']
        assert len(successful) == (2 if missing_target else 4)
        if not missing_target:
            resumed = orchestrator.run_study(study, {'resume': True})
            assert resumed['status'] == 'success', resumed
            assert resumed['evidence_nodes_count'] == first['evidence_nodes_count']
    assert registry.get('adata://s/annotated/v4')[1].n_obs == 128
