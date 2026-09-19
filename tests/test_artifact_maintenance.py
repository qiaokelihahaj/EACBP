from datetime import datetime, timedelta, timezone
import subprocess
import sys

import pytest

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.artifact.maintenance import cleanup_artifacts, plan_cleanup, ArtifactMaintenanceError
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskResult, TaskStatus


def fixture_registry(path):
    registry = ArtifactRegistry(str(path))
    registry.register('json://s/seed/v1', {'seed': True}, ArtifactType.JSON, 's', 'seed', 'fixture')
    return registry


def abandoned(registry, name='nested/output'):
    transaction = TaskArtifactTransaction(registry)
    uri = f'json://s/{name}/v1'
    transaction.register(uri, {'value': 1}, ArtifactType.JSON, 's', 't', 'fixture')
    return transaction, uri


def future():
    return datetime.now(timezone.utc) + timedelta(days=31)


def test_preview_then_apply_only_abandoned_transaction(tmp_path):
    registry = fixture_registry(tmp_path)
    tx, _ = abandoned(registry)
    tx.close()
    preview = cleanup_artifacts(registry, now=future())
    assert preview.planned_count == 1, preview.as_dict()
    assert preview.planned_bytes > 0
    path = preview.plan.eligible[0].path
    assert path.exists()
    result = cleanup_artifacts(registry, apply=True, plan=preview.plan, now=future())
    assert result.deleted_count == 1, result.as_dict()
    assert not path.exists()
    assert registry.get('json://s/seed/v1')[1] == {'seed': True}


def test_active_and_committed_transactions_are_retained(tmp_path):
    registry = fixture_registry(tmp_path)
    tx, uri = abandoned(registry)
    assert plan_cleanup(registry, now=future()).planned_count == 0
    result = TaskResult(task_id='t', capability='fixture', method_used='fixture',
                        status=TaskStatus.SUCCESS, output_artifacts=[uri])
    registry.commit_task(tx, 'signature', result)
    tx.close()
    report = cleanup_artifacts(registry, apply=True, now=future())
    assert report.deleted_count == 0
    assert registry.get(uri)[1] == {'value': 1}


def test_retention_unknown_files_and_stale_preview(tmp_path):
    registry = fixture_registry(tmp_path)
    tx, _ = abandoned(registry)
    tx.close()
    assert plan_cleanup(registry).planned_count == 0
    preview = plan_cleanup(registry, now=future())
    assert preview.planned_count == 1
    registry.register('json://s/new/v1', {}, ArtifactType.JSON, 's', 'new', 'fixture')
    report = cleanup_artifacts(registry, apply=True, plan=preview, now=future())
    assert report.stale_plan and report.deleted_count == 0
    path = preview.eligible[0].path
    (path / 'keep-user-file.txt').write_text('keep')
    assert cleanup_artifacts(registry, apply=True, now=future()).deleted_count == 0
    assert (path / 'keep-user-file.txt').read_text() == 'keep'


def test_missing_or_corrupt_registry_fails_closed(tmp_path):
    missing = tmp_path / 'missing'
    with pytest.raises(ArtifactMaintenanceError):
        plan_cleanup(missing)
    assert not missing.exists()
    registry = fixture_registry(tmp_path)
    tx, _ = abandoned(registry)
    tx.close()
    (tmp_path / registry.INDEX_FILENAME).write_text('{broken')
    report = cleanup_artifacts(registry, apply=True, now=future())
    assert report.blocked and report.deleted_count == 0


def test_process_death_releases_transaction_lease(tmp_path):
    registry = fixture_registry(tmp_path)
    code = '''
import os, sys
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.schemas.artifact import ArtifactType
tx = TaskArtifactTransaction(ArtifactRegistry(sys.argv[1]))
tx.register('json://s/crashed/v1', {}, ArtifactType.JSON, 's', 'crashed', 'fixture')
os._exit(23)
'''
    result = subprocess.run([sys.executable, '-c', code, str(tmp_path)], capture_output=True, timeout=30)
    assert result.returncode == 23, result.stderr
    report = cleanup_artifacts(registry, apply=True, now=future())
    assert report.deleted_count == 1, report.as_dict()
