"""Pin explicit external resource contents into task/resume signatures.

The planner stores the digest under the parameter name that the corresponding
adapter observes.  This is deliberately stricter than merely normalising a
path: an external executable or CellBender checkpoint is part of the method
contract, so changing its bytes must invalidate a saved task.
"""

from __future__ import annotations

import hashlib
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_file(value: Any, *, key: str) -> Path:
    """Resolve an explicitly supplied local file and fail loudly if invalid."""

    if value is None or not str(value).strip():
        raise ValueError(f"External resource {key} requires a non-empty local path")
    path = Path(str(value)).expanduser().resolve(strict=True)
    if not path.is_file():
        raise ValueError(f"External resource {key} must be a file: {path}")
    return path


def _resolve_executable(value: Any) -> Path:
    """Resolve a local executable path or an executable already on ``PATH``.

    CellBender's adapter accepts either spelling.  Resolving command names at
    planning time makes the selected executable content part of the resume
    signature and keeps the contract value compatible with the adapter's
    resolved provenance record.
    """

    if value is None or not str(value).strip():
        raise ValueError("External resource executable requires an explicit local executable")
    raw = str(value).strip()
    candidate = Path(raw).expanduser()
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError:
        located = shutil.which(raw)
        if located is None:
            raise FileNotFoundError(f"External executable does not exist or is not on PATH: {raw}")
        resolved = Path(located).resolve(strict=True)
    if not resolved.is_file():
        raise ValueError(f"External executable must be a file: {resolved}")
    return resolved


def _as_argument_list(value: Any, *, key: str) -> List[str]:
    if not isinstance(value, (list, tuple)) or any(not isinstance(item, str) for item in value):
        raise ValueError(f"CellBender {key} must be a list of strings")
    return list(value)


def _checkpoint_paths(parameters: Mapping[str, Any]) -> List[Tuple[str, Path]]:
    """Find local ``--checkpoint`` files without rewriting the CLI arguments.

    Relative paths are interpreted exactly as the CellBender subprocess does:
    relative to ``run_cwd``.  A logical key is used for each checkpoint so the
    adapter can expose the same key in its observed hash map.  The singular
    ``checkpoint`` key keeps the common one-checkpoint contract compact.
    """

    selected_key = None
    if "extra_args" in parameters:
        selected_key = "extra_args"
        raw_args = parameters["extra_args"]
    elif "cli_args" in parameters:
        selected_key = "cli_args"
        raw_args = parameters["cli_args"]
    else:
        return []

    args = _as_argument_list(raw_args, key=selected_key)
    run_cwd_value = parameters.get("run_cwd", Path.cwd())
    if run_cwd_value is None or not str(run_cwd_value).strip():
        raise ValueError("CellBender run_cwd must be a non-empty local directory")
    run_cwd = Path(str(run_cwd_value)).expanduser().resolve(strict=True)
    if not run_cwd.is_dir():
        raise ValueError(f"CellBender run_cwd must be a directory: {run_cwd}")

    paths: List[Tuple[str, Path]] = []
    index = 0
    while index < len(args):
        argument = args[index]
        checkpoint_value = None
        if argument == "--checkpoint":
            if index + 1 >= len(args):
                raise ValueError(f"CellBender {selected_key} has --checkpoint without a path")
            checkpoint_value = args[index + 1]
            index += 2
        elif argument.startswith("--checkpoint="):
            checkpoint_value = argument.split("=", 1)[1]
            index += 1
        else:
            index += 1
        if checkpoint_value is None:
            continue
        if not checkpoint_value.strip():
            raise ValueError(f"CellBender {selected_key} has an empty --checkpoint path")
        raw_path = Path(checkpoint_value).expanduser()
        if not raw_path.is_absolute():
            raw_path = run_cwd / raw_path
        checkpoint = raw_path.resolve(strict=True)
        if not checkpoint.is_file():
            raise ValueError(f"CellBender checkpoint must be a file: {checkpoint}")
        key = "checkpoint" if not paths else f"checkpoint[{len(paths)}]"
        paths.append((key, checkpoint))
    return paths


def pin_resource_files(parameters):
    """Return a copy with content hashes for explicit local resources.

    Every non-destination ``*_path`` is required to be an existing file.  The
    CellBender-specific ``executable`` and local ``--checkpoint`` arguments are
    handled as resources as well.  ``output_path`` is resolved but never
    hashed because it is produced by the task itself; its path remains in the
    task contract and therefore still participates in the resume fingerprint.

    No network access, package installation, or fallback resource discovery is
    performed here.  A command name already present on ``PATH`` is resolved to
    its concrete executable file so that the adapter can verify the same file.
    """

    if not isinstance(parameters, Mapping):
        raise ValueError("Capability parameters must be a mapping")
    result = dict(parameters)
    hashes: Dict[str, str] = {}

    existing_hashes = parameters.get("external_resource_sha256")
    if existing_hashes is not None:
        if not isinstance(existing_hashes, Mapping):
            raise ValueError("external_resource_sha256 must be a mapping of resource name to SHA-256")
        hashes.update({str(key): str(value) for key, value in existing_hashes.items()})

    for key, value in parameters.items():
        if not key.endswith("_path") or value is None:
            continue
        if key == "output_path":
            # A destination is not an input resource. Its existence/overwrite
            # policy belongs to the adapter; hashing it breaks first execution
            # and makes resume depend on a file produced by the task itself.
            if not str(value).strip():
                raise ValueError("External resource output_path requires a non-empty destination path")
            result[key] = str(Path(value).expanduser().resolve())
            continue
        path = _resolve_file(value, key=key)
        result[key] = str(path)
        hashes[key] = _sha256_file(path)

    if "executable" in parameters:
        executable = _resolve_executable(parameters["executable"])
        result["executable"] = str(executable)
        hashes["executable"] = _sha256_file(executable)

    for key, checkpoint in _checkpoint_paths(parameters):
        hashes[key] = _sha256_file(checkpoint)

    if hashes:
        result["external_resource_sha256"] = hashes
    return result
