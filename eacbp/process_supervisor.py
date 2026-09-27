"""Own an external command's process tree for its entire lifetime.

Windows jobs retain descendants even after their original parent exits. Start
the parent suspended so it cannot create an unowned child before job assignment.
POSIX commands receive a private process group, which can outlive its leader.
"""

from __future__ import annotations

import ctypes
import os
import signal
import subprocess


class _WindowsJob:
    def __init__(self):
        from ctypes import wintypes

        class BasicLimitInformation(ctypes.Structure):
            _fields_ = [
                ("PerProcessUserTimeLimit", ctypes.c_int64),
                ("PerJobUserTimeLimit", ctypes.c_int64),
                ("LimitFlags", wintypes.DWORD),
                ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t),
                ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t),
                ("PriorityClass", wintypes.DWORD),
                ("SchedulingClass", wintypes.DWORD),
            ]

        class IOCounters(ctypes.Structure):
            _fields_ = [(name, ctypes.c_uint64) for name in (
                "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                "ReadTransferCount", "WriteTransferCount", "OtherTransferCount",
            )]

        class ExtendedLimitInformation(ctypes.Structure):
            _fields_ = [
                ("BasicLimitInformation", BasicLimitInformation),
                ("IoInfo", IOCounters),
                ("ProcessMemoryLimit", ctypes.c_size_t),
                ("JobMemoryLimit", ctypes.c_size_t),
                ("PeakProcessMemoryUsed", ctypes.c_size_t),
                ("PeakJobMemoryUsed", ctypes.c_size_t),
            ]

        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.CreateJobObjectW.argtypes = [wintypes.LPVOID, wintypes.LPCWSTR]
        kernel.CreateJobObjectW.restype = wintypes.HANDLE
        kernel.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int,
                                                   wintypes.LPVOID, wintypes.DWORD]
        kernel.SetInformationJobObject.restype = wintypes.BOOL
        kernel.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        kernel.AssignProcessToJobObject.restype = wintypes.BOOL
        kernel.CloseHandle.argtypes = [wintypes.HANDLE]
        kernel.CloseHandle.restype = wintypes.BOOL
        native = ctypes.WinDLL("ntdll")
        native.NtResumeProcess.argtypes = [wintypes.HANDLE]
        native.NtResumeProcess.restype = wintypes.LONG
        self.kernel = kernel
        self.native = native
        self.handle = kernel.CreateJobObjectW(None, None)
        if not self.handle:
            raise ctypes.WinError(ctypes.get_last_error())
        limits = ExtendedLimitInformation()
        limits.BasicLimitInformation.LimitFlags = 0x00002000  # JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE
        if not kernel.SetInformationJobObject(self.handle, 9, ctypes.byref(limits), ctypes.sizeof(limits)):
            error = ctypes.WinError(ctypes.get_last_error())
            self.close()
            raise error

    def assign_and_resume(self, process):
        handle = int(process._handle)  # Popen owns the Windows process handle.
        if not self.kernel.AssignProcessToJobObject(self.handle, handle):
            raise ctypes.WinError(ctypes.get_last_error())
        status = self.native.NtResumeProcess(handle)
        if status < 0:
            raise OSError(f"Unable to resume supervised process (NTSTATUS 0x{status & 0xffffffff:08x})")

    def close(self):
        if self.handle:
            handle, self.handle = self.handle, None
            if not self.kernel.CloseHandle(handle):
                raise ctypes.WinError(ctypes.get_last_error())


class ProcessSupervisor:
    """Process handle plus durable ownership of descendants, without new deps."""

    def __init__(self, command, **kwargs):
        self.process = None
        self.job = _WindowsJob() if os.name == "nt" else None
        try:
            if self.job is not None:
                kwargs["creationflags"] = (kwargs.get("creationflags", 0)
                    | 0x00000004  # CREATE_SUSPENDED
                    | subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW)
            else:
                kwargs["start_new_session"] = True
            self.process = subprocess.Popen(command, **kwargs)
            if self.job is not None:
                self.job.assign_and_resume(self.process)
        except BaseException:
            self.terminate()
            raise

    def terminate(self):
        """Stop the owned tree, including descendants of an exited parent."""
        try:
            if self.job is not None:
                self.job.close()
            elif self.process is not None:
                try:
                    os.killpg(self.process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        finally:
            if self.process is not None:
                # Assignment can fail before the suspended process joins its
                # job. Explicitly kill it as well so failure never leaks it.
                if self.process.poll() is None:
                    self.process.kill()
                self.process.communicate()

    def close(self):
        # A completed command must not leave background scientific workers
        # running. This also closes the Windows job handle on the success path.
        self.terminate()
