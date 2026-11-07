"""Evidence CLI publication preserves immutable inputs and scientific failure."""
from types import SimpleNamespace
from pathlib import Path
import json

import pytest

from eacbp._atomic_files import publish_new_text
from eacbp.cli import main


def test_atomic_publication_never_replaces_existing_result(tmp_path):
    target = tmp_path / 'receipt.json'
    target.write_text('original evidence', encoding='utf-8')
    with pytest.raises(FileExistsError):
        publish_new_text(target, 'replacement')
    assert target.read_text(encoding='utf-8') == 'original evidence'


def test_competing_publication_preserves_winner_and_removes_staging_file(tmp_path, monkeypatch):
    import eacbp._atomic_files as publication
    target = tmp_path / 'receipt.json'
    def competing_link(source, destination):
        Path(destination).write_text('other completed version', encoding='utf-8')
        raise FileExistsError(destination)
    monkeypatch.setattr(publication.os, 'link', competing_link)
    with pytest.raises(FileExistsError):
        publish_new_text(target, 'new version')
    assert target.read_text(encoding='utf-8') == 'other completed version'
    assert list(tmp_path.iterdir()) == [target]


def test_failed_publication_never_exposes_partial_receipt(tmp_path, monkeypatch):
    import eacbp._atomic_files as publication
    def failed_link(*args):
        raise OSError('injected publication failure')
    monkeypatch.setattr(publication.os, 'link', failed_link)
    with pytest.raises(OSError, match='injected'):
        publish_new_text(tmp_path / 'receipt.json', '{"verified": true}\n')
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize('case', ['existing', 'input', 'same_outputs'])
def test_cli_rejects_output_aliases_before_verification(tmp_path, monkeypatch, case):
    import eacbp.application.verification as verification
    def forbidden(*args, **kwargs):
        pytest.fail('Output collision must be detected before evidence verification')
    monkeypatch.setattr(verification, 'verify_acceptance_index', forbidden)
    source = tmp_path / 'input.json'
    source.write_text('original input', encoding='utf-8')
    output = tmp_path / 'output.json'
    args = ['verify-evidence', '--index', str(source), '--index-sha256', 'a' * 64]
    if case == 'existing':
        output.write_text('original output', encoding='utf-8')
        args += ['--output', str(output)]
    elif case == 'input':
        args += ['--output', str(source)]
    else:
        args += ['--output', str(output), '--report', str(output)]
    assert main(args) == 1
    assert source.read_text(encoding='utf-8') == 'original input'
    if case == 'existing':
        assert output.read_text(encoding='utf-8') == 'original output'
    elif case == 'same_outputs':
        assert not output.exists()


@pytest.mark.parametrize('verification_status,exit_code', [('passed', 0), ('failed', 1)])
def test_scientific_failure_does_not_replace_integrity_exit_status(tmp_path, monkeypatch, capsys, verification_status, exit_code):
    import eacbp.application.verification as verification
    payload = {'verification_status': verification_status, 'science_combined': {'status': 'failed'},
               'failure_conditions': ['absolute ARI below frozen minimum']}
    result = SimpleNamespace(verification_status=verification_status,
                             model_dump=lambda **kwargs: payload)
    calls = []
    def verified(*args, **kwargs):
        calls.append((args, kwargs))
        return result
    monkeypatch.setattr(verification, 'verify_acceptance_index', verified)
    monkeypatch.setattr(verification, 'format_verification_summary', lambda value: 'Scientific hypothesis: failed\n')
    source, archive = tmp_path / 'input.json', tmp_path / 'archive.json'
    output, report = tmp_path / 'new.json', tmp_path / 'new.md'
    code = main(['verify-evidence', '--index', str(source), '--index-sha256', 'a' * 64,
                 '--source-archive', str(archive), '--source-archive-sha256', 'b' * 64,
                 '--output', str(output), '--report', str(report)])
    assert code == exit_code
    assert json.loads(output.read_text(encoding='utf-8')) == payload
    assert json.loads(capsys.readouterr().out) == payload
    assert report.read_text(encoding='utf-8') == 'Scientific hypothesis: failed\n'
    assert calls[0][1]['archive_receipt_path'] == archive
    assert calls[0][1]['expected_archive_receipt_sha256'] == 'b' * 64
    assert not source.exists() and not archive.exists()


@pytest.mark.parametrize('role', ['evidence', 'source', 'archive_original', 'archive_copy'])
@pytest.mark.parametrize('argument', ['--output', '--report'])
def test_missing_declared_inputs_cannot_be_used_as_failure_output(tmp_path, monkeypatch, role, argument):
    import eacbp.application.verification as verification
    def forbidden(*args, **kwargs):
        pytest.fail('Declared missing-input collision must be rejected before verification')
    monkeypatch.setattr(verification, 'verify_acceptance_index', forbidden)
    missing = tmp_path / 'missing-evidence.json'
    source, archive = tmp_path / 'index.json', tmp_path / 'archive.json'
    index, mapping = {}, {}
    if role == 'evidence':
        index['evidence_files'] = [{'path': str(missing), 'sha256': 'a' * 64}]
    elif role == 'source':
        index['source_files'] = {str(missing): 'a' * 64}
    else:
        mapping['sources'] = [{
            'original_path': str(missing if role == 'archive_original' else tmp_path / 'original.py'),
            'archive_path': str(missing if role == 'archive_copy' else tmp_path / 'copy.py')}]
    source.write_text(json.dumps(index), encoding='utf-8')
    archive.write_text(json.dumps(mapping), encoding='utf-8')
    before = source.read_bytes(), archive.read_bytes()
    assert main(['verify-evidence', '--index', str(source), '--index-sha256', 'a' * 64,
                 '--source-archive', str(archive), '--source-archive-sha256', 'b' * 64,
                 argument, str(missing)]) == 1
    assert not missing.exists()
    assert (source.read_bytes(), archive.read_bytes()) == before
