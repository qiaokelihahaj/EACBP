"""Run-wide resource policy and cooperative cancellation.

Python capabilities are checked at task boundaries. External commands are
supervised while running and their process trees are terminated on interruption.
This is deliberately not a claim of hard preemption for in-process Python code.
"""
from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
import math
import os
import subprocess
import threading
import time
from typing import Any, Mapping

from eacbp.process_supervisor import ProcessSupervisor


class ExecutionControlError(RuntimeError):
    failure_type = "resource_policy"


class ExecutionCancelled(ExecutionControlError):
    failure_type = "cancelled"


class ExecutionDeadlineExceeded(ExecutionControlError):
    failure_type = "timeout"


@dataclass
class ExecutionContext:
    deadline: float | None = None
    gpu_allowed: bool = True
    threads: int | None = None
    cancellation: threading.Event = field(default_factory=threading.Event)
    interruption: ExecutionControlError | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.deadline is not None and (
                isinstance(self.deadline, bool) or not isinstance(self.deadline, (int, float))
                or not math.isfinite(self.deadline)):
            raise ValueError("deadline must be a finite monotonic timestamp")
        if not isinstance(self.gpu_allowed, bool):
            raise ValueError("gpu_allowed must be a boolean")
        if self.threads is not None and (
                isinstance(self.threads, bool) or not isinstance(self.threads, int) or self.threads < 1):
            raise ValueError("threads must be a positive integer")

    @classmethod
    def from_constraints(cls, constraints, *, threads=None):
        hours = float(constraints.max_runtime_hours)
        if not math.isfinite(hours) or hours <= 0:
            raise ValueError("max_runtime_hours must be finite and positive")
        if threads is not None and (isinstance(threads, bool) or not isinstance(threads, int) or threads < 1):
            raise ValueError("threads must be a positive integer")
        return cls(time.monotonic() + hours * 3600, constraints.gpu_allowed, threads)

    def cancel(self):
        self.cancellation.set()

    def check(self):
        if self.cancellation.is_set():
            raise ExecutionCancelled("Study execution was cancelled")
        if self.deadline is not None and time.monotonic() >= self.deadline:
            raise ExecutionDeadlineExceeded("Study runtime limit exceeded")

    def validate_parameters(self, parameters):
        if self.gpu_allowed:
            return
        def visit(value):
            if isinstance(value, Mapping):
                for key, item in value.items():
                    name = str(key).casefold()
                    if name in {"cuda", "use_gpu", "gpu", "use_cuda"} and item:
                        raise ExecutionControlError("GPU use is forbidden by study constraints")
                    if name == "device" and str(item).casefold().startswith(("cuda", "gpu", "mps")):
                        raise ExecutionControlError("GPU device is forbidden by study constraints")
                    visit(item)
            elif isinstance(value, (list, tuple)):
                for item in value:
                    if isinstance(item, str) and item.split("=", 1)[0].casefold() in {"--cuda", "--gpu", "--use-gpu"}:
                        raise ExecutionControlError("GPU command flag is forbidden by study constraints")
                    visit(item)
        visit(parameters)


_active: ContextVar[ExecutionContext | None] = ContextVar("eacbp_execution_context", default=None)


@contextmanager
def activate_execution(context: ExecutionContext | None):
    token = _active.set(context)
    try:
        yield
    finally:
        _active.reset(token)


def current_execution_context():
    return _active.get()


def run_external(command, *, timeout=None, **kwargs):
    """subprocess.run-compatible subset used by the built-in scientific tools."""
    context = current_execution_context()
    if context is None:
        return subprocess.run(command, timeout=timeout, **kwargs)
    context.check()
    context.validate_parameters(command)
    if timeout is not None and (not math.isfinite(float(timeout)) or float(timeout) <= 0):
        raise ValueError("External command timeout must be finite and positive")
    limit = time.monotonic() + float(timeout) if timeout is not None else None
    capture = kwargs.pop("capture_output", False)
    check = kwargs.pop("check", False)
    if capture:
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    environment = dict(kwargs.pop("env", os.environ))
    if not context.gpu_allowed:
        environment.update(CUDA_VISIBLE_DEVICES="-1", HIP_VISIBLE_DEVICES="-1", ROCR_VISIBLE_DEVICES="-1")
    if context.threads is not None:
        environment.update({key: str(context.threads) for key in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})
    kwargs["env"] = environment
    supervisor = ProcessSupervisor(command, **kwargs)
    process = supervisor.process
    try:
        while True:
            context.check()
            if limit is not None and time.monotonic() >= limit:
                raise ExecutionDeadlineExceeded("External command timeout exceeded")
            try:
                stdout, stderr = process.communicate(timeout=0.2)
                break
            except subprocess.TimeoutExpired:
                continue
        context.check()
    except BaseException as exc:
        if isinstance(exc, ExecutionControlError):
            context.interruption = exc
        raise
    finally:
        supervisor.close()
    result = subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    if check and process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command, output=stdout, stderr=stderr)
    return result
