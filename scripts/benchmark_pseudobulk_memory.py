"""Measure sampled RSS for sparse pseudobulk preparation or audit input reconstruction.

Run this script in a fresh Python process. The generated sparse input has one
nonzero per cell, so its dense int64 equivalent is intentionally very large.
RSS is sampled during the selected operation; sampling can miss short-lived peaks and is
not a substitute for an operating-system high-water mark.
"""

from __future__ import annotations

import argparse
import ctypes
import json
import os
import threading
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy import sparse

from eacbp.capabilities.advanced_statistics import _prepare_pseudobulk
from eacbp.capabilities.sc_data import SCData


def _windows_rss_reader() -> Callable[[], int]:
    class ProcessMemoryCounters(ctypes.Structure):
        _fields_ = [
            ("cb", ctypes.c_ulong),
            ("PageFaultCount", ctypes.c_ulong),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    psapi = ctypes.WinDLL("Psapi.dll", use_last_error=True)
    get_process_memory_info = psapi.GetProcessMemoryInfo
    get_process_memory_info.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ProcessMemoryCounters),
        ctypes.c_ulong,
    ]
    get_process_memory_info.restype = ctypes.c_int
    kernel32 = ctypes.WinDLL("Kernel32.dll")
    get_current_process = kernel32.GetCurrentProcess
    get_current_process.restype = ctypes.c_void_p
    process = get_current_process()

    def read() -> int:
        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        if not get_process_memory_info(process, ctypes.byref(counters), counters.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return int(counters.WorkingSetSize)

    return read


def _rss_reader() -> Callable[[], int]:
    if os.name == "nt":
        return _windows_rss_reader()
    statm = Path("/proc/self/statm")
    if statm.is_file():
        page_size = int(os.sysconf("SC_PAGE_SIZE"))

        def read_proc_rss() -> int:
            resident_pages = int(statm.read_text(encoding="ascii").split()[1])
            return resident_pages * page_size

        return read_proc_rss
    try:
        import psutil
    except ImportError as exc:  # pragma: no cover - platform-specific fallback
        raise RuntimeError("RSS sampling needs /proc/self/statm or psutil on this platform") from exc
    process = psutil.Process()
    return lambda: int(process.memory_info().rss)


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=_positive_int, default=200_000)
    parser.add_argument("--genes", type=_positive_int, default=30_000)
    parser.add_argument("--sample-interval-ms", type=_positive_int, default=2)
    parser.add_argument("--operation", choices=("aggregate", "audit"), default="aggregate")
    args = parser.parse_args()
    if args.cells < 6:
        parser.error("--cells must be at least 6 to populate three paired donors")

    row = np.arange(args.cells, dtype=np.int64)
    col = row % args.genes
    counts = sparse.csr_matrix(
        (np.ones(args.cells, dtype=np.int64), (row, col)),
        shape=(args.cells, args.genes),
    )
    donor_labels = np.asarray(["d0", "d1", "d2"], dtype=object)
    obs = pd.DataFrame(
        {
            "donor_id": donor_labels[row % 3],
            "condition": np.where((row // 3) % 2 == 0, "A", "B"),
        }
    )
    data = SCData(
        X=sparse.csr_matrix((args.cells, args.genes), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"g{index}" for index in range(args.genes)]),
        layers={"counts": counts},
    )
    del row, col, counts, obs

    read_rss = _rss_reader()
    baseline_rss = read_rss()
    sampled_peak = [baseline_rss]
    stop = threading.Event()

    def sample_rss() -> None:
        while not stop.is_set():
            sampled_peak[0] = max(sampled_peak[0], read_rss())
            time.sleep(args.sample_interval_ms / 1000.0)

    sampler = threading.Thread(target=sample_rss, name="rss-sampler", daemon=True)
    sampler.start()
    params = {"condition_a": "A", "condition_b": "B", "donor_col": "donor_id", "paired": True}
    try:
        if args.operation == "audit":
            from eacbp.auditor.advanced_statistics import _reconstruct_input

            prepared = _reconstruct_input(data, params, {}, "pydeseq2_pseudobulk_v1")
        else:
            prepared = _prepare_pseudobulk(data, params)
    finally:
        stop.set()
        sampler.join()
    sampled_peak[0] = max(sampled_peak[0], read_rss())

    print(
        json.dumps(
            {
                "shape": [args.cells, args.genes],
                "operation": args.operation,
                "nnz": args.cells,
                "pseudobulk_shape": [len(prepared.metadata), args.genes],
                "dense_cell_matrix_bytes": args.cells * args.genes * np.dtype(np.int64).itemsize,
                "pseudobulk_dense_bytes": getattr(prepared, "pseudobulk_dense_bytes", None),
                "estimated_aggregation_working_bytes": getattr(prepared, "estimated_aggregation_working_bytes", None),
                "rss_baseline_bytes": baseline_rss,
                "sampled_peak_rss_bytes": sampled_peak[0],
                "observed_increment_bytes": sampled_peak[0] - baseline_rss,
                "sample_interval_ms": args.sample_interval_ms,
                "rss_source": "Windows WorkingSetSize" if os.name == "nt" else "process RSS sampler",
                "rss_measurement_limit": "sampled RSS may miss short-lived peaks; it is not an exact OS high-water mark",
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
