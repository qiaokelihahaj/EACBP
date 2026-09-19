import subprocess
import sys
from unittest.mock import patch
import pytest
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskResult, TaskStatus


def write_pair(tx):
    tx.register("json://s/a/v1", {"n": 1}, ArtifactType.JSON, "s", "t", "write")
    tx.register("json://s/b/v1", {"n": 2}, ArtifactType.JSON, "s", "t", "write", parent_uris=["json://s/a/v1"])
    return TaskResult(task_id="t", capability="test", method_used="test", status=TaskStatus.SUCCESS,
                      output_artifacts=["json://s/a/v1", "json://s/b/v1"])


@pytest.mark.parametrize("phase", ["first_output", "committed"])
def test_process_death_never_publishes_partial_task(tmp_path, phase):
    code = '''
import os, sys
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.transaction import TaskArtifactTransaction
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.task import TaskResult, TaskStatus
parent = ArtifactRegistry(sys.argv[1])
tx = TaskArtifactTransaction(parent)
tx.register("json://s/a/v1", {"n": 1}, ArtifactType.JSON, "s", "t", "write")
if sys.argv[2] == "first_output":
    os._exit(23)
tx.register("json://s/b/v1", {"n": 2}, ArtifactType.JSON, "s", "t", "write", parent_uris=["json://s/a/v1"])
result = TaskResult(task_id="t", capability="test", method_used="test", status=TaskStatus.SUCCESS,
                    output_artifacts=["json://s/a/v1", "json://s/b/v1"])
parent.commit_task(tx, "signature", result)
os._exit(23)
'''
    process = subprocess.run([sys.executable, "-c", code, str(tmp_path), phase], capture_output=True, timeout=40)
    assert process.returncode == 23, process.stderr
    parent = ArtifactRegistry(str(tmp_path))
    if phase == "first_output":
        assert parent.list_artifacts() == []
        assert parent.get_task_commit("signature") is None
        tx = TaskArtifactTransaction(parent)
        result = write_pair(tx)
        assert parent.list_artifacts() == []
        parent.commit_task(tx, "signature", result)
    assert len(parent.list_artifacts()) == 2
    assert parent.get("json://s/a/v1")[1] == {"n": 1}
    assert parent.get("json://s/b/v1")[1] == {"n": 2}
    assert set(parent.lineage.get_ancestors("json://s/b/v1")) == {"json://s/a/v1"}
    assert parent.get_task_commit("signature")["result"]["status"] == "success"


def test_failed_index_commit_preserves_existing_data_and_allows_retry(tmp_path):
    parent = ArtifactRegistry(str(tmp_path))
    parent.register("json://other/input/v1", {"stable": True}, ArtifactType.JSON, "other", "ingest", "write")
    tx = TaskArtifactTransaction(parent)
    result = write_pair(tx)
    with patch.object(parent, "_persist_index", side_effect=OSError("disk full")):
        with pytest.raises(OSError):
            parent.commit_task(tx, "signature", result)
    fresh = ArtifactRegistry(str(tmp_path))
    assert len(fresh.list_artifacts()) == 1
    assert fresh.get("json://other/input/v1")[1] == {"stable": True}
    assert fresh.get_task_commit("signature") is None
    fresh.commit_task(tx, "signature", result)
    assert len(fresh.list_artifacts()) == 3
