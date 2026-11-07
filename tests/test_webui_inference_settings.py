from types import SimpleNamespace

import pandas as pd
import pytest

from eacbp.webui.results import _deg_view


def test_deg_preview_preserves_actual_contrast_and_threshold():
    metadata = SimpleNamespace(
        uri="table://inference/deg/v1",
        parameters={"condition_a": "treated", "condition_b": "control"},
        summary_metrics={},
    )
    task = SimpleNamespace(
        task_id="deg", method_used="pydeseq2_pseudobulk_v1",
        metrics={"condition_a": "control", "condition_b": "treated",
                 "contrast_label": "control versus treated", "alpha": 0.1,
                 "effect_definition": "PyDESeq2 log2 fold change for control versus treated"},
    )
    table = pd.DataFrame({"gene": ["g"], "log2_fold_change": [-1.0], "fdr_q_value": [0.08]})
    view = _deg_view(metadata, table, task)
    assert view["condition_a"] == "control"
    assert view["contrast_label"] == "control versus treated"
    assert view["effect_definition"] == task.metrics["effect_definition"]
    assert view["alpha"] == 0.1


@pytest.mark.parametrize('capability,method,missing_package', [
    ('rna_velocity', 'scvelo_stochastic_v1', 'scvelo'),
    ('liana_communication', 'liana_rank_aggregate_v1', 'liana'),
    ('doublet_detection', 'scanpy_scrublet_v1', 'scikit-image'),
])
def test_preflight_checks_optional_packages_from_method_descriptors(tmp_path, monkeypatch,
                                                                  capability, method, missing_package):
    from importlib.metadata import PackageNotFoundError
    from eacbp.webui.service import WebService

    # Isolate package preflight from scientific execution and matrix reading.
    source = tmp_path / 'input.h5ad'
    source.write_bytes(b'preflight fixture; not a scientific matrix')
    service = WebService(tmp_path, tmp_path / 'runs')
    plan = {'valid': True, 'tasks': [{'capability': capability, 'method': method}]}
    monkeypatch.setattr('eacbp.orchestrator.planning.preview_study_plan', lambda *args, **kwargs: plan)
    monkeypatch.setattr(service, 'design', lambda *args, **kwargs: {
        'overview': {'samples': {'count': 0}}, 'can_run': True,
    })

    def installed_version(name):
        if name == missing_package:
            raise PackageNotFoundError(name)
        return 'test-installed-version'

    monkeypatch.setattr('eacbp.webui.service.version', installed_version)
    result = service.preview({'data': str(source), 'species': 'human', 'tissue': 'test'})
    assert missing_package in service.settings()['packages']
    assert result['missing_packages'] == [missing_package]
    assert result['can_run'] is False
    assert not list((tmp_path / 'runs').glob('*/run_config.json'))

    monkeypatch.setattr('eacbp.webui.service.version', lambda name: 'test-installed-version')
    result = service.preview({'data': str(source), 'species': 'human', 'tissue': 'test'})
    assert result['missing_packages'] == []
    assert result['can_run'] is True
