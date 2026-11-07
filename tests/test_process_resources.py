"""Evidence semantics plus a real memory-allocating scientific grandchild."""
import json
import os
import subprocess
import sys
import time

import pytest

from eacbp.application import process_resources as resources
from eacbp.process_supervisor import ProcessSupervisor


class FakeProbe:
    name = "test_probe"
    creation_time_unit = "test_ticks"

    def __init__(self):
        self.parents = {100: 1, 101: 100, 102: 101}
        self.values = {100: (10, 4), 101: (11, 20), 102: (12, 70)}

    def processes(self):
        return self.parents.copy()

    def read(self, pid):
        value = self.values[pid]
        if isinstance(value, Exception):
            raise value
        return {"pid": pid, "creation_time": value[0], "rss_bytes": value[1]}


def _fake(monkeypatch, tmp_path):
    probe = FakeProbe()
    monkeypatch.setattr(resources, "_make_probe", lambda: probe)
    sampler = resources.ProcessTreeSampler(tmp_path / "attempt", root_pid=100, interval_seconds=30)
    return sampler, probe


def _rows(sampler):
    return [json.loads(line) for line in sampler.sample_path.read_text(encoding="utf-8").splitlines()]


def test_peak_is_maximum_of_sample_sums_and_paths_are_exclusive(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    with sampler:
        assert len(_rows(sampler)) == 1  # First probe is immediate and durable.
        with pytest.raises(FileExistsError):
            resources.ProcessTreeSampler(sampler.directory).start()
        probe.values[101], probe.values[102] = (11, 80), (12, 5)
        sampler._sample_safely()
    assert [row["sum_rss_bytes"] for row in _rows(sampler)] == [94, 89]
    assert sampler.summary["observed_peak_sum_rss_bytes"] == 94  # Not 4+80+70.
    assert sampler.summary["max_observed_processes"] == 3
    assert sampler.summary["coverage"] == "sampled_observations"
    assert sampler.summary["successful_sample_count"] == 2
    assert sampler.summary["root_creation_time"] == 10
    assert sampler.summary["elapsed_seconds"] >= 0
    assert sampler.summary["observed_max_sample_gap_seconds"] > 0
    assert json.loads(sampler.summary_path.read_text()) == sampler.summary
    assert not sampler.summary_path.with_suffix(".json.tmp").exists()


def test_known_grandchild_remains_after_parent_exits(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    with sampler:
        del probe.parents[101]
        sampler._sample_safely()
    assert {record["pid"] for record in _rows(sampler)[1]["processes"]} == {100, 102}
    assert _rows(sampler)[1]["sum_rss_bytes"] == 74


def test_root_pid_reuse_never_admits_replacement_tree(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    with sampler:
        probe.values[100] = (20, 9999)
        probe.parents[200], probe.values[200] = 100, (21, 9999)
        sampler._sample_safely()
        sampler._sample_safely()
    assert [row["sum_rss_bytes"] for row in _rows(sampler)] == [94, None, None]
    assert sampler.summary["observed_peak_sum_rss_bytes"] == 94
    assert sampler.summary["coverage"] == "partial"
    assert 200 not in {record["pid"] for record in sampler.summary["process_identities"]}
    assert any(error["phase"] == "root_pid_reused" for error in sampler.summary["errors"])


def test_descendant_pid_reuse_does_not_admit_its_new_children(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    with sampler:
        probe.values[101] = (20, 9999)
        probe.parents[200], probe.values[200] = 101, (21, 9999)
        sampler._sample_safely()
    # The original known grandchild remains independently identity-verified.
    assert _rows(sampler)[1]["sum_rss_bytes"] == 74
    assert {record["pid"] for record in _rows(sampler)[1]["processes"]} == {100, 102}
    assert 200 not in {record["pid"] for record in sampler.summary["process_identities"]}
    assert sampler.summary["errors"][-1]["pid"] == 101


def test_all_process_probes_failed_never_fabricate_zero(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    probe.values[100] = PermissionError(13, "Denied")
    with sampler:
        pass
    assert sampler.summary["observed_peak_sum_rss_bytes"] is None
    assert sampler.summary["successful_sample_count"] == 0
    assert sampler.summary["coverage"] == "unavailable"
    assert sampler.summary["root_creation_time"] is None
    assert _rows(sampler)[0]["sum_rss_bytes"] is None
    assert _rows(sampler)[0]["errors"][0]["pid"] == 100


def test_probe_errors_preserve_partial_rss_and_diagnostic_attribution(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    probe.values[102] = PermissionError(13, "Denied")
    with sampler:
        pass
    row = _rows(sampler)[0]
    assert row["sum_rss_bytes"] == 24
    assert row["status"] == "partial"
    assert row["errors"] == [{"phase": "process_probe", "pid": 102,
                              "error_type": "PermissionError", "os_error_code": 13}]
    assert sampler.summary["error_free_sample_count"] == 0
    assert sampler.summary["coverage"] == "partial"


def test_no_success_is_null_and_unsupported_platform_is_explicit(monkeypatch, tmp_path):
    monkeypatch.setattr(resources.sys, "platform", "unsupported-test")
    with resources.ProcessTreeSampler(tmp_path / "unsupported") as sampler:
        pass
    assert sampler.summary["backend"] == "unavailable"
    assert sampler.summary["coverage"] == "unavailable"
    assert sampler.summary["observed_peak_sum_rss_bytes"] is None
    assert sampler.summary["successful_sample_count"] == 0
    assert _rows(sampler)[0]["sum_rss_bytes"] is None
    assert sampler.summary["errors"][0]["error_type"] == "NotImplementedError"


def test_enumeration_failure_preserves_verified_root(monkeypatch, tmp_path):
    sampler, probe = _fake(monkeypatch, tmp_path)
    def fail():
        raise OSError(5, "Denied")
    monkeypatch.setattr(probe, "processes", fail)
    with sampler:
        pass
    assert sampler.summary["observed_peak_sum_rss_bytes"] == 4
    assert sampler.summary["coverage"] == "partial"
    assert sampler.summary["errors"][0]["phase"] == "tree_enumeration"


def test_exception_stops_thread_and_keeps_final_receipt(monkeypatch, tmp_path):
    sampler, _ = _fake(monkeypatch, tmp_path)
    with pytest.raises(RuntimeError, match="science failed"):
        with sampler:
            raise RuntimeError("science failed")
    assert not sampler._thread.is_alive()
    assert sampler.summary_path.exists()
    assert sampler.summary["elapsed_seconds"] >= 0
    assert len(_rows(sampler)) == 1
    assert sampler.stop() is sampler.summary


def test_formatting_failure_keeps_prior_log_and_original_exception(monkeypatch, tmp_path):
    sampler, _ = _fake(monkeypatch, tmp_path)
    with pytest.raises(RuntimeError, match="science failed"):
        with sampler:
            def fail(*args, **kwargs):
                raise TypeError("Cannot format")
            monkeypatch.setattr(resources.json, "dumps", fail)
            sampler._sample_safely()
            assert sampler._stop.is_set()
            raise RuntimeError("science failed")
    assert len(_rows(sampler)) == 1
    assert sampler.summary["errors"][-1]["phase"] == "sample_persistence_or_format"
    assert sampler.summary_path.exists()


def test_summary_formatter_failure_preserves_sample_log_and_science_exception(monkeypatch, tmp_path):
    sampler, _ = _fake(monkeypatch, tmp_path)
    with pytest.raises(RuntimeError, match="science failed"):
        with sampler:
            def fail(*args, **kwargs):
                raise TypeError("Cannot format summary")
            monkeypatch.setattr(resources.json, "dump", fail)
            raise RuntimeError("science failed")
    assert not sampler._thread.is_alive()
    assert len(_rows(sampler)) == 1
    assert not sampler.summary_path.exists()
    assert sampler.summary["summary_output_error"]["error_type"] == "TypeError"


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_interval_rejected(tmp_path, value):
    with pytest.raises(ValueError):
        resources.ProcessTreeSampler(tmp_path, interval_seconds=value)


def test_native_scientific_grandchild_working_set_is_recorded(tmp_path):
    if not (sys.platform == "win32" or sys.platform.startswith("linux")):
        with resources.ProcessTreeSampler(tmp_path / "native") as sampler:
            pass
        assert sampler.summary["coverage"] == "unavailable"
        assert sampler.summary["observed_peak_sum_rss_bytes"] is None
        return
    allocation = 64 * 1024 * 1024
    grandchild_code = (
        "import json,os,time; data=bytearray(" + str(allocation) + ");"
        "data[::4096]=b'x'*len(data[::4096]);"
        "print(json.dumps({'role':'grandchild','pid':os.getpid()}),flush=True);time.sleep(3)"
    )
    child_code = (
        "import subprocess,sys; child=subprocess.Popen([sys.executable,'-B','-c'," + repr(grandchild_code) +
        "],stdin=subprocess.DEVNULL);child.wait()"
    )
    parent_code = (
        "import json,os,subprocess,sys;print(json.dumps({'role':'root','pid':os.getpid()}),flush=True);"
        "input();child=subprocess.Popen([sys.executable,'-B','-c'," + repr(child_code) +
        "],stdin=subprocess.DEVNULL);child.wait()"
    )
    supervisor = ProcessSupervisor([sys.executable, "-B", "-c", parent_code],
                                   stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True)
    try:
        worker = supervisor.process
        root_pid = json.loads(worker.stdout.readline())["pid"]
        with resources.ProcessTreeSampler(tmp_path / "native", root_pid=root_pid,
                                         interval_seconds=.03) as sampler:
            worker.stdin.write("start\n")
            worker.stdin.flush()
            grandchild_pid = json.loads(worker.stdout.readline())["pid"]
            deadline = time.monotonic() + 2
            observed = []
            while time.monotonic() < deadline:
                observed = [record for row in _rows(sampler) for record in row["processes"]
                            if record["pid"] == grandchild_pid and record["rss_bytes"] >= allocation]
                if observed:
                    break
                time.sleep(.03)
            assert observed, "The allocating scientific grandchild must be measured, not only its launcher"
        rows = _rows(sampler)
        assert all(row["sum_rss_bytes"] == sum(record["rss_bytes"] for record in row["processes"])
                   for row in rows if row["sum_rss_bytes"] is not None)
        assert sampler.summary["root_pid"] == root_pid
        assert sampler.summary["max_observed_processes"] >= 3
        assert sampler.summary["observed_peak_sum_rss_bytes"] >= allocation
        assert sampler.summary["successful_sample_count"] >= 2
        assert all("cmdline" not in record for row in rows for record in row["processes"])
    finally:
        supervisor.close()
