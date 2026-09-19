from concurrent.futures import ThreadPoolExecutor, TimeoutError
from threading import Event
from unittest.mock import patch
import subprocess
import sys
import pandas as pd
import pytest
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.orchestrator.checkpoint import StudyJournal
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.auditor.base import ValidationReport
from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import StudyManifest, BiologicalDesign
from eacbp.schemas.task import TaskContract, TaskResult, TaskStatus


def test_registry_parallel_processes_preserve_all_commits(tmp_path):
    code = '''
import sys
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.schemas.artifact import ArtifactType
reg = ArtifactRegistry(sys.argv[1])
sys.stdin.readline()
for i in range(2):
    reg.register(f"json://s/worker{sys.argv[2]}_{i}/v1", {"worker": sys.argv[2]}, ArtifactType.JSON, "s", "t", "write")
'''
    workers = [subprocess.Popen([sys.executable, "-c", code, str(tmp_path), str(i)],
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True) for i in range(4)]
    try:
        for worker in workers:
            worker.stdin.write("start\n")
            worker.stdin.flush()
        for worker in workers:
            out, err = worker.communicate(timeout=45)
            assert worker.returncode == 0, out + err
        reg = ArtifactRegistry(str(tmp_path))
        assert len(reg.list_artifacts()) == 8
        for i in range(4):
            for j in range(2):
                assert reg.get(f"json://s/worker{i}_{j}/v1")[1] == {"worker": str(i)}
    finally:
        for worker in workers:
            if worker.poll() is None:
                worker.kill()
            worker.communicate()


def test_study_lock_released_after_process_crash(tmp_path):
    code = '''
import os, sys
from eacbp.orchestrator.checkpoint import StudyJournal
journal = StudyJournal(sys.argv[1], "s")
journal.__enter__()
print("locked", flush=True)
sys.stdin.readline()
os._exit(23)
'''
    worker = subprocess.Popen([sys.executable, "-c", code, str(tmp_path)],
                              stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, text=True)
    try:
        assert worker.stdout.readline().strip() == "locked"
        with pytest.raises(RuntimeError, match="already running"):
            with StudyJournal(tmp_path, "s"):
                pass
        worker.communicate("crash\n", timeout=20)
        assert worker.returncode == 23
        with StudyJournal(tmp_path, "s") as journal:
            journal.entries = {"recovered": True}
            journal.save()
        with StudyJournal(tmp_path, "s") as journal:
            assert journal.entries == {"recovered": True}
    finally:
        if worker.poll() is None:
            worker.kill()
        worker.communicate()


def test_same_instance_reader_cannot_erase_inflight_registration(tmp_path):
    reg = ArtifactRegistry(str(tmp_path))
    committing, release, reading = Event(), Event(), Event()
    original = reg._persist_index
    def commit():
        committing.set()
        assert release.wait(5)
        original()
    def read():
        reading.set()
        return reg.list_artifacts()
    with patch.object(reg, "_persist_index", side_effect=commit), ThreadPoolExecutor(2) as pool:
        writer = pool.submit(reg.register, "json://s/a/v1", {"a": 1}, ArtifactType.JSON, "s", "t", "write")
        assert committing.wait(5)
        reader = pool.submit(read)
        assert reading.wait(5)
        try:
            with pytest.raises(TimeoutError):
                reader.result(timeout=.1)
        finally:
            release.set()
        writer.result(timeout=5)
        assert len(reader.result(timeout=5)) == 1
    assert ArtifactRegistry(str(tmp_path)).get("json://s/a/v1")[1] == {"a": 1}


@pytest.mark.parametrize("failure", ["audit", "commit", "computed"])
def test_computed_outputs_resume_after_failure_without_admitting_evidence(tmp_path, failure):
    reg = ArtifactRegistry(str(tmp_path))
    orch = ScientificOrchestrator(artifact_registry=reg)
    manifest = StudyManifest(study_id="s", biological_design=BiologicalDesign(species="test", tissue="test"))
    uri = "table://s/audit/v1"
    task = TaskContract(task_id="audit", capability="dataset_audit", expected_outputs=[uri])
    def execute(task, registry):
        registry.register(uri, pd.DataFrame({"n_cells": [3]}), ArtifactType.TABLE, "s", "audit", "audit")
        return TaskResult(task_id="audit", capability="dataset_audit", method_used="sc_audit_v1",
                          status=TaskStatus.SUCCESS, output_artifacts=[uri], metrics={"n_cells": 3})
    original_save = StudyJournal.save
    def save(journal):
        if failure == "commit" and journal.entries["audit"].get("phase") == "audited":
            raise OSError("injected journal write failure")
        if failure == "computed":
            raise OSError("injected first journal write failure")
        original_save(journal)
    audit = ValidationReport(auditor_name="test", target_task_id="audit")
    with patch.object(ComputationalDAGPlanner, "build_study_plan", return_value=[task]), \
         patch.object(orch.capability_registry, "execute_contract", side_effect=execute) as compute:
        with patch.object(StudyJournal, "save", save), \
             patch.object(orch.auditor, "audit_task", side_effect=RuntimeError("injected audit crash") if failure == "audit" else None, return_value=audit):
            first = orch.run_study(manifest)
        assert first["status"] == "failed"
        assert first["evidence_nodes_count"] == first["claims_count"] == 0
        assert "n_cells" not in orch.current_state
        with patch.object(orch.auditor, "audit_task", return_value=audit):
            resumed = orch.run_study(manifest, {"resume": True})
        assert resumed["status"] == "success", resumed["failures"]
        assert compute.call_count == 1
        assert resumed["evidence_nodes_count"] == 1
