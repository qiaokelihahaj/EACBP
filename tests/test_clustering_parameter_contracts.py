"""Effective Leiden options reject ignored or misspelled scientific settings."""
import pytest
from pydantic import ValidationError

from eacbp.capabilities.parameters import LeidenParameters, NoCorrectionParameters, HarmonyParameters
from eacbp.orchestrator.planning import preview_study_plan
from eacbp.schemas.study import StudyManifest, BiologicalDesign


@pytest.mark.parametrize("options", [
    {"resoluton": 1.0}, {"n_pcs": 10}, {"k_clusters": 4},
    {"resolution": True}, {"resolution": 0.0}, {"resolution": float('nan')},
    {"n_neighbors": 1}, {"n_neighbors": True}, {"n_components": 0},
])
def test_leiden_rejects_ignored_or_invalid_options(options):
    with pytest.raises(ValidationError):
        LeidenParameters.model_validate(options)


def test_integration_contracts_preserve_supported_actual_options():
    assert NoCorrectionParameters(n_components=30, batch_col=None).n_components == 30
    assert HarmonyParameters(n_components=30, batch_col="library", max_iter_harmony=12).max_iter_harmony == 12
    with pytest.raises(ValidationError):
        NoCorrectionParameters(n_pcs=30)
    with pytest.raises(ValidationError):
        NoCorrectionParameters(zero_center=True)


def study():
    return StudyManifest(study_id="cluster_contracts", biological_design=BiologicalDesign(
        species="human", tissue="pancreas", target_cell_types=[]))


def test_standard_preview_removes_only_the_injected_kmeans_default():
    result = preview_study_plan(study(), {"method_profile": "standard",
        "method_overrides": {"integration": "no_correction_v1"},
        "capability_parameters": {"clustering": {"resolution": 0.5, "n_neighbors": 15}}})
    assert result['valid'], result['errors']
    task = next(task for task in result['tasks'] if task['capability'] == 'clustering')
    assert task['method'] == 'scanpy_leiden_umap_v1'
    assert 'k_clusters' not in task['parameters']
    assert task['parameters']['resolution'] == 0.5
    assert 'leiden_execution_integrity' in task['validation_requirements']
    invalid = preview_study_plan(study(), {'method_profile': 'standard',
        'capability_parameters': {'clustering': {'k_clusters': 4}}})
    assert not invalid['valid']
    assert any('k_clusters' in error['message'] for error in invalid['errors'])


def test_baseline_kmeans_contract_retains_its_own_algorithm_option():
    result = preview_study_plan(study(), {'method_profile': 'baseline',
        'capability_parameters': {'clustering': {'k_clusters': 6}}})
    assert result['valid'], result['errors']
    task = next(task for task in result['tasks'] if task['capability'] == 'clustering')
    assert task['method'] == 'kmeans_marker_embedding_v1'
    assert task['parameters']['k_clusters'] == 6
