"""A process tree remains owned after the command's original parent exits."""

import os
import subprocess
import sys
import threading
import time

import pytest

from eacbp.execution_context import (
    ExecutionCancelled, ExecutionContext, ExecutionDeadlineExceeded,
    activate_execution, run_external,
)
from eacbp.process_supervisor import ProcessSupervisor, _WindowsJob


def _parent_command(marker, parent_finished):
    child = "import pathlib,time;time.sleep(1.2);pathlib.Path(" + repr(str(marker)) + ").write_text('escaped')"
    parent = ("import pathlib,subprocess,sys;subprocess.Popen([sys.executable,'-c'," + repr(child)
              + "]);pathlib.Path(" + repr(str(parent_finished)) + ").write_text('spawned')")
    return [sys.executable, "-c", parent]


@pytest.mark.parametrize("interruption", ["timeout", "cancel"])
def test_exited_parent_does_not_let_child_survive_interruption(tmp_path, interruption):
    marker, parent_finished = tmp_path / "child.txt", tmp_path / "parent.txt"
    context = ExecutionContext()
    timer = threading.Timer(0.4, context.cancel) if interruption == "cancel" else None
    if timer is not None:
        timer.start()
    try:
        with activate_execution(context), pytest.raises(
            ExecutionCancelled if interruption == "cancel" else ExecutionDeadlineExceeded
        ):
            run_external(_parent_command(marker, parent_finished),
                         timeout=0.4 if interruption == "timeout" else None, capture_output=True)
    finally:
        if timer is not None:
            timer.cancel()
    assert parent_finished.exists(), "The parent must spawn its child before the interruption"
    time.sleep(1.3)
    assert not marker.exists(), "An orphaned worker must not write after cancellation/timeout"


def test_completed_command_cleans_up_background_workers(tmp_path):
    marker, parent_finished = tmp_path / "child.txt", tmp_path / "parent.txt"
    with activate_execution(ExecutionContext()):
        completed = run_external(_parent_command(marker, parent_finished),
                                 stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert completed.returncode == 0
    assert parent_finished.exists()
    time.sleep(1.3)
    assert not marker.exists()


@pytest.mark.skipif(os.name != "nt", reason="Windows must assign the job before resuming the process")
def test_failed_job_assignment_kills_suspended_parent(tmp_path, monkeypatch):
    marker = tmp_path / "must-not-start.txt"
    processes = []
    popen = subprocess.Popen

    def track_process(*args, **kwargs):
        process = popen(*args, **kwargs)
        processes.append(process)
        return process

    def reject_assignment(self, process):
        raise OSError("simulated job assignment failure")

    monkeypatch.setattr(subprocess, "Popen", track_process)
    monkeypatch.setattr(_WindowsJob, "assign_and_resume", reject_assignment)
    with pytest.raises(OSError, match="assignment failure"):
        ProcessSupervisor([sys.executable, "-c", "import pathlib;pathlib.Path(" + repr(str(marker)) + ").write_text('bad')"],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert len(processes) == 1 and processes[0].poll() is not None
    assert not marker.exists()


@pytest.mark.parametrize("kwargs", [
    {"deadline": float("nan")}, {"deadline": float("inf")}, {"deadline": True},
    {"deadline": "later"}, {"gpu_allowed": "false"}, {"gpu_allowed": 0},
    {"threads": 0}, {"threads": -1}, {"threads": 1.5}, {"threads": True},
])
def test_direct_execution_context_rejects_invalid_policy(kwargs):
    with pytest.raises(ValueError):
        ExecutionContext(**kwargs)


def test_expired_finite_context_remains_valid_to_construct():
    context = ExecutionContext(deadline=-1, gpu_allowed=False, threads=1)
    with pytest.raises(ExecutionDeadlineExceeded):
        context.check()
