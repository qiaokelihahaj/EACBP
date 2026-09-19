"""Small shared primitive for durable JSON publication.

Callers keep ownership of their locks, serialization defaults, and domain
specific error handling.  This helper only stages one JSON document, flushes
and fsyncs it, replaces the destination, and removes an abandoned temporary
file on every failure.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping


def atomic_write_json(
    path: str | Path,
    payload: Any,
    *,
    dump_kwargs: Mapping[str, Any] | None = None,
    create_parent: bool = False,
    newline: str | None = None,
    trailing_newline: bool = False,
) -> Path:
    """Atomically publish *payload* as JSON and return its destination.

    ``dump_kwargs`` is passed directly to :func:`json.dump`; this preserves
    each existing caller's serializer, ``allow_nan`` and key-order policy.
    The temporary file is created beside the destination with ``mkstemp`` so
    replacement remains on one filesystem and concurrent writers do not share
    a predictable staging path.
    """

    target = Path(path)
    if not target.name:
        raise ValueError("JSON path must name a file")
    if create_parent:
        target.parent.mkdir(parents=True, exist_ok=True)

    temporary: Path | None = None
    fd: int | None = None
    try:
        fd, name = tempfile.mkstemp(
            prefix=f".{target.name}.", suffix=".tmp", dir=str(target.parent)
        )
        temporary = Path(name)
        handle = os.fdopen(fd, "w", encoding="utf-8", newline=newline)
        fd = None
        with handle:
            json.dump(payload, handle, **dict(dump_kwargs or {}))
            if trailing_newline:
                handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
        temporary = None
        return target
    finally:
        if fd is not None:
            try:
                os.close(fd)
            except OSError:
                pass
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass


__all__ = ["atomic_write_json"]
