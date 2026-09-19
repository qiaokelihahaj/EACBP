"""Private task outputs; publication is one atomic parent-index update.

Uncommitted directories are intentionally retained for diagnosis. They are not
visible through the parent registry and do not reserve public artifact URIs.

Each transaction also owns a small, OS-backed lease.  The lease is deliberately
held for the lifetime of the transaction rather than being inferred from file
mtimes: an artifact maintainer can therefore distinguish an abandoned
transaction (the OS released its lease after a crash) from one that is still
being written.  The parent storage lock serializes creation with maintenance
enumeration, so a maintainer never observes a half-created transaction.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from pathlib import Path
from types import SimpleNamespace
from uuid import uuid4

from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.storage import ArtifactAlreadyExistsError


class TransactionLeaseError(RuntimeError):
    """Raised when a transaction lease cannot be acquired or released."""


def _lock_exclusive_nonblocking(handle) -> None:
    """Acquire one byte in *handle* without waiting for another owner."""

    handle.seek(0)
    if os.name == "nt":
        import msvcrt

        msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
    else:  # pragma: no cover - exercised on POSIX CI
        import fcntl

        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)


def _unlock(handle) -> None:
    handle.seek(0)
    if os.name == "nt":
        import msvcrt

        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    else:  # pragma: no cover - exercised on POSIX CI
        import fcntl

        fcntl.flock(handle.fileno(), fcntl.LOCK_UN
        )


class TaskArtifactTransaction(ArtifactRegistry):
    LEASE_FILENAME = ".transaction.lease"
    LEASE_VERSION = 1

    def __init__(self, parent):
        self.parent = parent
        self.transaction_id = uuid4().hex
        self.transaction_path = (
            parent.storage.base_dir / "_transactions" / self.transaction_id
        )
        self.lease_path = self.transaction_path / self.LEASE_FILENAME
        self._lease_handle = None
        self._closed = False

        # Parent maintenance takes this same lock before enumerating
        # ``_transactions``.  Keep it through both directory creation and
        # lease publication so there is never an unleased visible attempt.
        with parent.storage.lock():
            super().__init__(str(self.transaction_path))
            self._acquire_lease()

    def _acquire_lease(self) -> None:
        self.transaction_path.mkdir(parents=True, exist_ok=True)
        handle = None
        try:
            # ``a+b`` implies O_APPEND, which makes an in-place marker update
            # impossible after seeking back to byte zero.  The parent storage
            # lock serializes creation, so touch then open read/write is safe.
            self.lease_path.touch(exist_ok=True)
            handle = self.lease_path.open("r+b")
            if handle.seek(0, os.SEEK_END) == 0:
                handle.write(b"0")
                handle.flush()
            _lock_exclusive_nonblocking(handle)
        except (OSError, ValueError) as exc:
            try:
                handle.close()
            except (UnboundLocalError, AttributeError):
                pass
            raise TransactionLeaseError(
                f"Unable to acquire transaction lease at {self.lease_path}"
            ) from exc

        self._lease_handle = handle
        payload = {
            "lease_version": self.LEASE_VERSION,
            "transaction_id": self.transaction_id,
            "pid": os.getpid(),
            "created_at": datetime.now(timezone.utc).isoformat(),
        }
        try:
            # The lock is held while rewriting the marker.  The marker is
            # informational; the open handle is the liveness authority.
            encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
            handle.seek(0)
            handle.truncate()
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        except OSError as exc:
            self.close()
            raise TransactionLeaseError(
                f"Unable to publish transaction lease at {self.lease_path}"
            ) from exc

    @property
    def closed(self) -> bool:
        return self._closed

    def close(self) -> None:
        """Release the transaction lease; safe to call more than once."""

        if self._closed:
            return
        self._closed = True
        handle = self._lease_handle
        self._lease_handle = None
        if handle is None:
            return
        try:
            _unlock(handle)
        except (OSError, ValueError):
            # Closing an already-invalid handle still releases the OS lease.
            pass
        finally:
            try:
                handle.close()
            except OSError:
                pass

    def __enter__(self):
        return self

    def __exit__(self, *_exc_info):
        self.close()
        return False

    def __del__(self):  # pragma: no cover - interpreter shutdown fallback
        try:
            self.close()
        except Exception:
            pass

    def _inherit_provenance(self, parent_uris, summary_metrics):
        combined = dict(self.registry)
        for uri in parent_uris:
            if uri not in combined and self.parent.exists(uri):
                combined[uri] = self.parent.get_metadata(uri)
        return ArtifactRegistry._inherit_provenance(SimpleNamespace(registry=combined), parent_uris, summary_metrics)

    def register(self, uri_str, *args, **kwargs):
        if self._closed:
            raise TransactionLeaseError("Cannot write to a closed artifact transaction")
        if self.parent.exists(uri_str):
            raise ArtifactAlreadyExistsError(f"Artifact '{uri_str}' is already registered")
        return super().register(uri_str, *args, **kwargs)

    def get_metadata(self, uri_str):
        try:
            return super().get_metadata(uri_str)
        except KeyError:
            return self.parent.get_metadata(uri_str)

    def get(self, uri_str):
        if super().exists(uri_str):
            return super().get(uri_str)
        return self.parent.get(uri_str)

    def load_payload(self, uri_str):
        return self.get(uri_str)[1]

    def exists(self, uri_str, artifact_type=None):
        return super().exists(uri_str, artifact_type) or self.parent.exists(uri_str, artifact_type)
