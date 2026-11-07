"""Publish complete UTF-8 outputs without replacing an existing version."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile


def publish_new_text(path: str | Path, text: str) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        raise FileExistsError(target)
    descriptor, name = tempfile.mkstemp(prefix='.' + target.name + '.', suffix='.tmp', dir=target.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8', newline='') as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, target)
        return target
    finally:
        temporary.unlink(missing_ok=True)
