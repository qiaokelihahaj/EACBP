"""Durable, sampled process-tree RSS evidence using only the standard library.

The sum includes shared pages once per process. It is neither a physical-memory
peak nor a memory limit. Sequential probes can miss short peaks and descendants
whose parent disappeared before discovery; already observed descendants persist.
"""
from __future__ import annotations

import ctypes
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
import threading
import time


class _WindowsProbe:
    name = "windows_toolhelp32_working_set"
    creation_time_unit = "100ns_since_1601_utc"

    def __init__(self):
        from ctypes import wintypes as w
        class Entry(ctypes.Structure):
            _fields_ = [("size", w.DWORD), ("usage", w.DWORD), ("pid", w.DWORD),
                        ("heap", ctypes.c_size_t), ("module", w.DWORD), ("threads", w.DWORD),
                        ("ppid", w.DWORD), ("priority", w.LONG), ("flags", w.DWORD),
                        ("exe", w.WCHAR * 260)]
        class Memory(ctypes.Structure):
            _fields_ = [("cb", w.DWORD), ("faults", w.DWORD)] + [
                (name, ctypes.c_size_t) for name in ("peak", "rss", "peak_paged", "paged",
                                                     "peak_nonpaged", "nonpaged", "commit", "peak_commit")]
        self.Entry, self.Memory = Entry, Memory
        self.kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        self.psapi = ctypes.WinDLL("psapi", use_last_error=True)
        definitions = {
            "CreateToolhelp32Snapshot": ([w.DWORD, w.DWORD], w.HANDLE),
            "Process32FirstW": ([w.HANDLE, ctypes.POINTER(Entry)], w.BOOL),
            "Process32NextW": ([w.HANDLE, ctypes.POINTER(Entry)], w.BOOL),
            "OpenProcess": ([w.DWORD, w.BOOL, w.DWORD], w.HANDLE),
            "CloseHandle": ([w.HANDLE], w.BOOL),
            "GetProcessTimes": ([w.HANDLE] + [ctypes.POINTER(w.FILETIME)] * 4, w.BOOL),
        }
        for name, (args, result) in definitions.items():
            function = getattr(self.kernel, name)
            function.argtypes, function.restype = args, result
        self.psapi.GetProcessMemoryInfo.argtypes = [w.HANDLE, ctypes.POINTER(Memory), w.DWORD]
        self.psapi.GetProcessMemoryInfo.restype = w.BOOL

    def processes(self):
        handle = self.kernel.CreateToolhelp32Snapshot(2, 0)
        if handle == ctypes.c_void_p(-1).value:
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            entry = self.Entry()
            entry.size = ctypes.sizeof(entry)
            if not self.kernel.Process32FirstW(handle, ctypes.byref(entry)):
                raise ctypes.WinError(ctypes.get_last_error())
            result = {}
            while True:
                result[int(entry.pid)] = int(entry.ppid)
                if not self.kernel.Process32NextW(handle, ctypes.byref(entry)):
                    if ctypes.get_last_error() != 18:  # ERROR_NO_MORE_FILES
                        raise ctypes.WinError(ctypes.get_last_error())
                    return result
        finally:
            self.kernel.CloseHandle(handle)

    def read(self, pid):
        from ctypes import wintypes as w
        handle = self.kernel.OpenProcess(0x1000, False, pid)  # QUERY_LIMITED_INFORMATION
        if not handle:
            raise ctypes.WinError(ctypes.get_last_error())
        try:
            times = [w.FILETIME() for _ in range(4)]
            if not self.kernel.GetProcessTimes(handle, *(ctypes.byref(value) for value in times)):
                raise ctypes.WinError(ctypes.get_last_error())
            memory = self.Memory()
            memory.cb = ctypes.sizeof(memory)
            if not self.psapi.GetProcessMemoryInfo(handle, ctypes.byref(memory), memory.cb):
                raise ctypes.WinError(ctypes.get_last_error())
            created = (times[0].dwHighDateTime << 32) | times[0].dwLowDateTime
            return {"pid": pid, "creation_time": created, "rss_bytes": int(memory.rss)}
        finally:
            self.kernel.CloseHandle(handle)


class _LinuxProbe:
    name = "linux_proc_stat_rss"
    creation_time_unit = "clock_ticks_since_boot"

    @staticmethod
    def _stat(pid):
        value = Path(f"/proc/{pid}/stat").read_text(encoding="utf-8")
        return value[value.rfind(")") + 2:].split()

    def processes(self):
        result = {}
        for path in Path("/proc").iterdir():
            if path.name.isdecimal():
                try:
                    result[int(path.name)] = int(self._stat(path.name)[1])
                except FileNotFoundError:  # Unrelated processes can exit during enumeration.
                    continue
        return result

    def read(self, pid):
        fields = self._stat(pid)
        rss = int(fields[21])
        if rss < 0:
            raise ValueError("Invalid negative RSS")
        return {"pid": pid, "creation_time": int(fields[19]),
                "rss_bytes": rss * os.sysconf("SC_PAGE_SIZE")}


def _make_probe():
    if sys.platform == "win32":
        return _WindowsProbe()
    if sys.platform.startswith("linux"):
        return _LinuxProbe()
    raise NotImplementedError("Process-tree RSS is unavailable on this platform")


def _error(phase, error, pid=None):
    return {"phase": phase, "pid": pid, "error_type": type(error).__name__,
            "os_error_code": getattr(error, "winerror", None) or getattr(error, "errno", None)}


class _IdentityMismatch(ValueError):
    pass


class ProcessTreeSampler:
    """An attempt's exclusive fsynced JSONL and atomic final resource summary."""

    def __init__(self, directory, root_pid=None, interval_seconds=.25):
        if not math.isfinite(interval_seconds) or interval_seconds <= 0:
            raise ValueError("interval_seconds must be finite and positive")
        self.directory = Path(directory)
        self.root_pid = os.getpid() if root_pid is None else int(root_pid)
        if self.root_pid <= 0:
            raise ValueError("root_pid must be positive")
        self.interval_seconds = float(interval_seconds)
        self.sample_path = self.directory / "process_resources.samples.jsonl"
        self.summary_path = self.directory / "process_resources.summary.json"
        self.summary = None
        self._stop = threading.Event()
        self._thread = self._stream = None
        self._known = {}
        self._root_creation = None
        self._root_reused = False
        self._errors = []
        self._error_count = 0
        self._samples = self._successful = self._complete = 0
        self._peak = None
        self._max_processes = 0
        self._last_elapsed = self._max_gap = None

    def start(self):
        if self._stream is not None or self.summary is not None:
            raise RuntimeError("A sampler cannot be restarted")
        self.directory.mkdir(parents=True, exist_ok=True)
        if self.summary_path.exists():
            raise FileExistsError(self.summary_path)
        self._stream = self.sample_path.open("x", encoding="utf-8")
        self._started = time.monotonic()
        self._started_utc = datetime.now(timezone.utc).isoformat()
        try:
            self._probe = _make_probe()
        except Exception as error:
            self._probe = None
            self._record_errors([_error("probe_initialization", error)])
        self._sample_safely()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def _run(self):
        while not self._stop.wait(self.interval_seconds):
            self._sample_safely()

    def _sample_safely(self):
        try:
            self._sample()
        except Exception as error:
            self._record_errors([_error("sample_persistence_or_format", error)])
            self._stop.set()  # Retain earlier durable rows if formatting or writing fails.

    def _record_errors(self, errors):
        self._error_count += len(errors)
        self._errors = (self._errors + errors)[-100:]

    def _sample(self):
        records, errors, parents = {}, [], {}
        def read(pid):
            if pid in records:
                return records[pid]
            try:
                record = self._probe.read(pid)
                known = self._known.get(pid)
                if known and known["creation_time"] != record["creation_time"]:
                    raise _IdentityMismatch("process PID reused")
                if record["rss_bytes"] < 0:
                    raise ValueError("Invalid negative RSS")
                records[pid] = record
                return record
            except Exception as error:
                errors.append(_error("process_probe", error, pid))
                return None
        if self._probe is not None and not self._root_reused:
            root = read(self.root_pid)
            if root is not None:
                if self._root_creation is None:
                    self._root_creation = root["creation_time"]
                elif root["creation_time"] != self._root_creation:
                    self._root_reused = True
            elif any(e["pid"] == self.root_pid and e["error_type"] == "_IdentityMismatch" for e in errors):
                self._root_reused = True
            if self._root_reused:
                records.clear()
                errors.append({"phase": "root_pid_reused", "pid": self.root_pid,
                               "error_type": "IdentityMismatch", "os_error_code": None})
            else:
                try:
                    parents = self._probe.processes()
                except Exception as error:
                    errors.append(_error("tree_enumeration", error))
                for pid in self._known:
                    if pid != self.root_pid and pid in parents:
                        read(pid)  # Keep verified descendants after their original parent exits.
                pending = list(records)
                while pending:
                    parent = pending.pop()
                    for pid, ppid in parents.items():
                        if ppid == parent and pid not in records:
                            child = read(pid)
                            if child and child["creation_time"] >= records[parent]["creation_time"]:
                                pending.append(pid)
                            elif child:
                                records.pop(pid)
                                errors.append({"phase": "parent_identity_order", "pid": pid,
                                               "error_type": "IdentityMismatch", "os_error_code": None})
        for pid, record in records.items():
            record["parent_pid"] = parents.get(pid)
            self._known[pid] = record.copy()
        total = sum(record["rss_bytes"] for record in records.values()) if records else None
        row_errors = errors if self._probe else self._errors[:1] if self._samples == 0 else []
        row = {"schema_version": 1, "sample_index": self._samples,
               "backend": self._probe.name if self._probe else "unavailable",
               "creation_time_unit": self._probe.creation_time_unit if self._probe else None,
               "root_creation_time": self._root_creation, "root_pid_reused": self._root_reused,
               "utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.monotonic() - self._started, "root_pid": self.root_pid,
               "processes": sorted(records.values(), key=lambda value: value["pid"]),
               "sum_rss_bytes": total, "errors": row_errors,
               "status": "partial" if total is not None and errors else "observed" if total is not None else "unavailable"}
        self._stream.write(json.dumps(row, allow_nan=False) + "\n")
        self._stream.flush()
        os.fsync(self._stream.fileno())
        if self._last_elapsed is not None:
            self._max_gap = max(self._max_gap or 0, row["elapsed_seconds"] - self._last_elapsed)
        self._last_elapsed = row["elapsed_seconds"]
        self._samples += 1
        self._record_errors(errors)
        if total is not None:
            self._successful += 1
            self._complete += int(not errors)
            self._peak = max(self._peak or 0, total)
            self._max_processes = max(self._max_processes, len(records))

    def stop(self):
        if self.summary is not None:
            return self.summary
        if self._stream is None:
            raise RuntimeError("The sampler has not started")
        self._stop.set()
        if self._thread is not None:
            self._thread.join()
        try:
            self._stream.close()
        except Exception as error:
            self._record_errors([_error("sample_close", error)])
        self.summary = {
            "schema_version": 1, "backend": self._probe.name if self._probe else "unavailable",
            "root_pid": self.root_pid, "root_creation_time": self._root_creation,
            "creation_time_unit": self._probe.creation_time_unit if self._probe else None,
            "started_utc": self._started_utc, "elapsed_seconds": time.monotonic() - self._started,
            "ended_utc": datetime.now(timezone.utc).isoformat(),
            "sampling_interval_seconds": self.interval_seconds, "sample_count": self._samples,
            "observed_max_sample_gap_seconds": self._max_gap,
            "successful_sample_count": self._successful, "error_free_sample_count": self._complete,
            "observed_peak_sum_rss_bytes": self._peak, "max_observed_processes": self._max_processes,
            "coverage": "unavailable" if not self._successful else "partial" if self._errors else "sampled_observations",
            "measurement_scope": "Root Python process and identity-verified observed descendants; launch parents excluded. "
                "Non-atomic sequential sampled SUM RSS includes shared pages repeatedly, misses short peaks and unobserved orphan descendants; "
                "not a physical-memory peak, complete job accounting, or enforced memory limit.",
            "process_identities": sorted(self._known.values(), key=lambda value: value["pid"]),
            "error_count": self._error_count, "errors": self._errors,
            "sample_path": str(self.sample_path),
        }
        temporary = self.summary_path.with_suffix(".json.tmp")
        try:
            with temporary.open("x", encoding="utf-8") as stream:
                json.dump(self.summary, stream, allow_nan=False, indent=2)
                stream.write("\n")
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, self.summary_path)
        except Exception as error:
            self.summary["summary_output_error"] = _error("summary_persistence_or_format", error)
        return self.summary

    __enter__ = start

    def __exit__(self, exc_type, exc, tb):
        self.stop()
        return False
