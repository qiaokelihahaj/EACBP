import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from eacbp import cli
from eacbp.application import study_service
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.events import read_run_events
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import BiologicalDesign, StudyManifest


def study():
    return StudyManifest(study_id='cli_test', biological_design=BiologicalDesign(species='human', tissue='brain'))


def input_file(tmp_path):
    import anndata as ad
    obs = pd.DataFrame({
        'condition': ['case'] * 12 + ['control'] * 12,
        'donor': [f'd{i // 3}' for i in range(24)],
    }, index=[f'cell{i}' for i in range(24)])
    data = ad.AnnData(np.random.default_rng(5).poisson(5, (24, 8)).astype(float), obs=obs,
                      var=pd.DataFrame(index=[f'g{i}' for i in range(8)]))
    path = tmp_path / 'input.h5ad'
    data.write_h5ad(path)
    return path


@pytest.fixture
def audit_only(monkeypatch):
    original = ComputationalDAGPlanner.build_study_plan
    def plan(manifest, config, **kwargs):
        return [t for t in original(manifest, config, **kwargs) if t.capability == 'dataset_audit']
    monkeypatch.setattr(ComputationalDAGPlanner, 'build_study_plan', plan)


def test_cli_plan_is_read_only_and_reports_errors(tmp_path, monkeypatch, capsys):
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(study().model_dump_json())
    monkeypatch.chdir(tmp_path)
    before = set(tmp_path.iterdir())
    with patch.object(study_service, '_import_h5ad', side_effect=AssertionError('must not import')), \
         patch.object(study_service, '_make_registry', side_effect=AssertionError('must not create storage')):
        assert cli.main(['plan', '--manifest', str(manifest)]) == 0
    assert json.loads(capsys.readouterr().out)['phase'] == 'before_dataset_audit'
    assert set(tmp_path.iterdir()) == before
    assert cli.main(['plan', '--manifest', str(manifest), '--config', '{"method_overrides":{"deg":"missing"}}']) == 1
    assert not json.loads(capsys.readouterr().out)['valid']


def test_cli_real_import_strict_resume_and_snapshot_report(tmp_path, monkeypatch, audit_only):
    source = input_file(tmp_path)
    run_dir = tmp_path / 'run'
    first = cli.run_study(manifest=study(), config={'method_profile': 'baseline'}, data=source, run_dir=run_dir)
    assert first['status'] == 'success', first
    original_report = (run_dir / 'report.md').read_text(encoding='utf-8')
    assert (run_dir / 'snapshots' / (first['run_id'] + '.json')).is_file()
    saved = json.loads((run_dir / 'run_config.json').read_text())
    assert saved['import_completed'] is True
    assert saved['source']['sha256_before_import'] == saved['source']['sha256_after_import']
    source.write_bytes(b'original source changed; imported registry copy is authoritative')
    elsewhere = tmp_path / 'elsewhere'
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    with patch.object(study_service, '_import_h5ad', side_effect=AssertionError('resume cannot re-import')):
        rebuilt = cli.report_study(run_dir=run_dir, output=tmp_path / 'rebuilt.md')
        assert Path(rebuilt['report']).read_text(encoding='utf-8') == original_report
        resumed = cli.resume_study(run_dir=run_dir)
    assert resumed['status'] == 'success', resumed
    kinds = [event.kind for event in read_run_events(resumed['event_log'])]
    assert 'task_resumed' in kinds and 'attempt_started' not in kinds
    assert len(list((run_dir / 'snapshots').glob('*.json'))) == 2
    assert json.loads((run_dir / 'run_config.json').read_text())['config'] == saved['config']


def test_run_refuses_existing_directory_without_modifying_it(tmp_path):
    folder = tmp_path / 'existing'
    folder.mkdir()
    marker = folder / 'keep.txt'
    marker.write_text('untouched')
    with pytest.raises(ValueError, match='already exists'):
        cli.run_study(manifest=study(), config={}, data=tmp_path / 'missing', run_dir=folder)
    assert marker.read_text() == 'untouched'
    assert set(folder.iterdir()) == {marker}


def test_snapshot_failure_makes_cli_delivery_fail(tmp_path, audit_only, monkeypatch):
    source = input_file(tmp_path)
    run_dir = tmp_path / 'run'
    with patch.object(study_service, '_persist_snapshot', side_effect=OSError('snapshot disk full')):
        with pytest.raises(OSError, match='snapshot disk full'):
            cli.run_study(manifest=study(), config={'method_profile': 'baseline'}, data=source, run_dir=run_dir)
    summary = json.loads((run_dir / 'summary.json').read_text())
    assert summary['status'] == 'failed'
    assert 'snapshot disk full' in summary['error']
    assert (run_dir / 'artifacts' / ArtifactRegistry.INDEX_FILENAME).is_file()


def test_report_rejects_changed_payload_without_overwriting_report(tmp_path, audit_only):
    run_dir = tmp_path / 'run'
    cli.run_study(manifest=study(), config={'method_profile': 'baseline'}, data=input_file(tmp_path), run_dir=run_dir)
    old = (run_dir / 'report.md').read_bytes()
    registry = ArtifactRegistry(str(run_dir / 'artifacts'))
    raw = registry.get_metadata('adata://cli_test/raw/v1')
    Path(raw.storage_path).write_bytes(b'corrupt')
    assert cli.main(['report', '--run-dir', str(run_dir)]) == 1
    assert (run_dir / 'report.md').read_bytes() == old


def test_resume_rejects_missing_registry_and_future_schema(tmp_path):
    path = tmp_path / 'run'
    path.mkdir()
    envelope = {'schema_version': 999, 'import_completed': True}
    config = path / 'run_config.json'
    config.write_text(json.dumps(envelope))
    assert cli.main(['resume', '--run-dir', str(path)]) == 1
    envelope['schema_version'] = 1
    config.write_text(json.dumps(envelope))
    assert cli.main(['resume', '--run-dir', str(path)]) == 1
    assert not (path / 'artifacts').exists()


def test_cleanup_corrupt_index_returns_failure(tmp_path, capsys):
    registry = ArtifactRegistry(str(tmp_path))
    registry.register('json://s/input/v1', {}, ArtifactType.JSON, 's', 'input', 'import')
    (tmp_path / ArtifactRegistry.INDEX_FILENAME).write_text('{corrupt')
    assert cli.main(['cleanup', str(tmp_path)]) == 1
    assert json.loads(capsys.readouterr().out)['blocked'] is True


def test_config_paths_preserve_cell_ids_and_checkpoint_cwd(tmp_path):
    config = {'capability_parameters': {
        'background_removal': {'run_cwd': 'worker', 'extra_args': ['--checkpoint', 'model.pt']},
        'clustering': {'marker_reference': {'network_path': ['GeneA', 'GeneB']}},
    }}
    result = cli._normalise_config_paths(config, base_dir=tmp_path)
    params = result['capability_parameters']
    assert params['background_removal']['run_cwd'] == str((tmp_path / 'worker').resolve())
    assert params['background_removal']['extra_args'][1] == str((tmp_path / 'worker/model.pt').resolve())
    assert params['clustering'] == config['capability_parameters']['clustering']

def test_cli_bundle_move_rebuilds_report_without_original_directory(tmp_path, audit_only, capsys):
    original = tmp_path / 'original'
    cli.run_study(manifest=study(), config={'method_profile': 'baseline'},
                  data=input_file(tmp_path), run_dir=original)
    bundle = tmp_path / 'study.zip'
    # Historical snapshots are attachments, not claims of current admission.
    (original / 'snapshots' / 'historical.json').write_text('{"historical": true}')
    assert cli.main(['export', '--run-dir', str(original), '--output', str(bundle)]) == 0
    capsys.readouterr()
    restored = tmp_path / 'restored'
    assert cli.main(['import', '--bundle', str(bundle), '--run-dir', str(restored)]) == 0
    capsys.readouterr()
    assert Path(cli.inspect_run(run_dir=restored)['event_log']).is_relative_to(restored)
    assert (restored / 'snapshots' / 'historical.json').read_text() == '{"historical": true}'
    original.rename(tmp_path / 'original_offline')
    assert cli.main(['report', '--run-dir', str(restored), '--output', str(restored / 'rebuilt.md')]) == 0
    rebuilt = (restored / 'rebuilt.md').read_text(encoding='utf-8')
    assert 'cli_test' in rebuilt
    assert 'Evidence' in rebuilt
    # Original metadata remains provenance, not permission to re-import it.
    saved = json.loads((restored / 'run_config.json').read_text())
    assert saved['import_completed'] is True
    with pytest.raises(FileExistsError):
        study_service.import_study(bundle=bundle, run_dir=restored)


def test_application_service_accepts_cancellation_context(tmp_path, audit_only):
    from eacbp.execution_context import ExecutionContext
    context = ExecutionContext(threads=2)
    context.cancel()
    result = study_service.run_study(manifest=study(), config={'method_profile': 'baseline'},
                                    data=input_file(tmp_path), run_dir=tmp_path / 'cancelled',
                                    execution_context=context)
    assert result['status'] == 'failed'
    snapshot = json.loads((tmp_path / 'cancelled' / 'snapshot.json').read_text())
    assert 'cancelled' in json.dumps(snapshot)


@pytest.mark.parametrize('invalid', ['incomplete', 'version', 'manifest'])
def test_application_export_rejects_invalid_envelope_without_publishing(tmp_path, audit_only, invalid):
    from eacbp.artifact.portability import ArtifactBundleError
    run_dir = tmp_path / 'run'
    study_service.run_study(manifest=study(), config={'method_profile': 'baseline'},
                            data=input_file(tmp_path), run_dir=run_dir)
    config_path = run_dir / 'run_config.json'
    envelope = json.loads(config_path.read_text())
    if invalid == 'incomplete':
        envelope['import_completed'] = False
    elif invalid == 'version':
        envelope['schema_version'] = 999
    else:
        envelope['manifest']['biological_design']['tissue'] = 'heart'
    config_path.write_text(json.dumps(envelope))
    output = tmp_path / 'invalid-study.zip'
    with pytest.raises(ArtifactBundleError):
        study_service.export_study(run_dir=run_dir, output=output)
    assert not output.exists()
