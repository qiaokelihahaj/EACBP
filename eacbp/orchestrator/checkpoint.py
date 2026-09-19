"""Atomic local task journals for explicit, integrity-checked resume."""
import hashlib
import json
import os
from pathlib import Path
from datetime import datetime
from enum import Enum
import numpy as np
from eacbp._atomic_json import atomic_write_json


def json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, datetime):
        return value.isoformat()
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Unsupported checkpoint value: {type(value).__name__}")


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=json_default).encode()).hexdigest()


def execution_environment():
    """Identify installed distributions and the actual package source for resume."""
    import platform
    from importlib.metadata import distributions
    package = Path(__file__).resolve().parents[1]
    sources = {str(path.relative_to(package)): hashlib.sha256(path.read_bytes()).hexdigest()
               for path in sorted(package.rglob("*.py"))}
    versions = {dist.metadata["Name"]: dist.version for dist in distributions() if dist.metadata["Name"]}
    return {"python": platform.python_version(), "platform": platform.platform(),
            "distributions": versions, "source_sha256": fingerprint(sources)}


class StudyJournal:
    def __init__(self, storage_dir, study_id):
        # The journal name cannot escape the storage root even for external study IDs.
        self.path = Path(storage_dir) / "_runs" / (fingerprint(study_id) + ".json")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock = self.path.with_suffix(".lock")
        self.entries = {}

    def __enter__(self):
        self._lock_handle = self.lock.open("a+b")
        if self._lock_handle.tell() == 0:
            self._lock_handle.write(b"0")
            self._lock_handle.flush()
        self._lock_handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(self._lock_handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(self._lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            self._lock_handle.close()
            raise RuntimeError("This study is already running in another executor.") from exc
        if self.path.exists():
            try:
                self.entries = json.loads(self.path.read_text(encoding="utf-8"))
            except Exception:
                self.__exit__()
                raise
        return self

    def save(self):
        atomic_write_json(
            self.path,
            self.entries,
            dump_kwargs={"default": json_default, "indent": 2, "allow_nan": False},
        )

    def __exit__(self, *args):
        # OS releases this lock even if the process terminates unexpectedly.
        # Keep the lock file to avoid an unlink/open race with another executor.
        try:
            self._lock_handle.seek(0)
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(self._lock_handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(self._lock_handle.fileno(), fcntl.LOCK_UN)
        finally:
            self._lock_handle.close()
