"""Conservative cleanup for abandoned private artifact transactions.

The public artifact tree is immutable and is never a cleanup target.  This
module only considers direct children of ``<storage>/_transactions`` that have
the transaction shape produced by :class:`TaskArtifactTransaction`.  A
candidate is deleted only when its lease is no longer held, it is older than
the requested retention period, and no metadata, computation receipt, audit
receipt, or study journal points at it.

Planning and application both run under the parent storage lock.  An
application supplied with a preview plan compares the complete state snapshot
again before deleting anything; a stale plan is therefore a no-op.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import threading
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

from eacbp.artifact.lineage import LineageGraph
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.storage import ArtifactStorageBackend
from eacbp.artifact.uri import ArtifactURI
from eacbp.schemas.artifact import ArtifactMetadata


TRANSACTION_ID_RE = re.compile(r"^[0-9a-f]{32}$")
LEASE_FILENAME = ".transaction.lease"
INDEX_FILENAME = ArtifactRegistry.INDEX_FILENAME
LOCK_FILENAME = ".registry.lock"
DEFAULT_RETENTION_DAYS = 30


class ArtifactMaintenanceError(RuntimeError):
    """Base error for a maintenance operation that cannot be made safe."""


class RegistryUnavailable(ArtifactMaintenanceError):
    """Raised internally when the parent registry cannot be trusted."""


@dataclass(frozen=True)
class CleanupCandidate:
    """One transaction directory observed during planning."""

    transaction_id: str
    path: Path
    size_bytes: int
    age_days: float
    eligible: bool
    reason: str
    snapshot: str = ""
    references: tuple[str, ...] = ()

    @property
    def relative_path(self) -> str:
        return str(self.path)

    def as_dict(self) -> dict[str, Any]:
        return {
            "transaction_id": self.transaction_id,
            "path": str(self.path),
            "size_bytes": self.size_bytes,
            "age_days": round(self.age_days, 6),
            "eligible": self.eligible,
            "reason": self.reason,
            "references": list(self.references),
        }


@dataclass(frozen=True)
class CleanupPlan:
    """A preview that can be supplied to :func:`cleanup_artifacts`."""

    base_dir: Path
    older_than_days: int
    generated_at: datetime
    candidates: tuple[CleanupCandidate, ...] = ()
    blocked_reasons: tuple[str, ...] = ()
    snapshot_token: str = ""
    registry_available: bool = True

    @property
    def eligible(self) -> tuple[CleanupCandidate, ...]:
        return tuple(candidate for candidate in self.candidates if candidate.eligible)

    @property
    def planned_count(self) -> int:
        return len(self.eligible)

    @property
    def planned_bytes(self) -> int:
        return sum(candidate.size_bytes for candidate in self.eligible)

    @property
    def blocked(self) -> bool:
        return bool(self.blocked_reasons)

    def as_dict(self) -> dict[str, Any]:
        return {
            "base_dir": str(self.base_dir),
            "older_than_days": self.older_than_days,
            "generated_at": self.generated_at.isoformat(),
            "registry_available": self.registry_available,
            "planned_count": self.planned_count,
            "planned_bytes": self.planned_bytes,
            "blocked_reasons": list(self.blocked_reasons),
            "candidates": [candidate.as_dict() for candidate in self.candidates],
        }


@dataclass
class CleanupReport:
    """Result of a preview or an application attempt."""

    mode: str
    plan: CleanupPlan
    planned_count: int = 0
    planned_bytes: int = 0
    deleted_count: int = 0
    deleted_bytes: int = 0
    skipped: list[dict[str, Any]] = field(default_factory=list)
    reasons: dict[str, int] = field(default_factory=dict)
    stale_plan: bool = False
    errors: list[str] = field(default_factory=list)

    @property
    def blocked(self) -> bool:
        return bool(self.plan.blocked_reasons) or self.stale_plan or bool(self.errors)

    @property
    def bytes_reclaimed(self) -> int:
        return self.deleted_bytes

    def as_dict(self) -> dict[str, Any]:
        return {
            "mode": self.mode,
            "planned_count": self.planned_count,
            "planned_bytes": self.planned_bytes,
            "deleted_count": self.deleted_count,
            "deleted_bytes": self.deleted_bytes,
            "bytes_reclaimed": self.bytes_reclaimed,
            "stale_plan": self.stale_plan,
            "blocked_reasons": list(self.plan.blocked_reasons),
            "reasons": dict(self.reasons),
            "skipped": list(self.skipped),
            "errors": list(self.errors),
            "candidates": [candidate.as_dict() for candidate in self.plan.candidates],
        }

    def __getitem__(self, key: str) -> Any:
        return self.as_dict()[key]


@dataclass
class _TransactionInfo:
    transaction_id: str
    path: Path
    metadata: dict[str, ArtifactMetadata]
    snapshot: str
    size_bytes: int
    age_days: float
    valid: bool
    reason: str = ""
    references: set[str] = field(default_factory=set)


def _utc_now(value: Optional[datetime] = None) -> datetime:
    current = value or datetime.now(timezone.utc)
    if current.tzinfo is None:
        return current.replace(tzinfo=timezone.utc)
    return current.astimezone(timezone.utc)


def _normalise_base(storage_dir: str | Path | ArtifactRegistry) -> Path:
    if isinstance(storage_dir, ArtifactRegistry):
        return storage_dir.storage.base_dir.resolve()
    return Path(storage_dir).expanduser().resolve()


def _is_reparse_or_symlink(path: Path) -> bool:
    """Return true for symlinks and Windows junction/reparse points."""

    try:
        info = os.lstat(path)
    except OSError:
        return False
    if stat.S_ISLNK(info.st_mode):
        return True
    attributes = getattr(info, "st_file_attributes", 0)
    reparse = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0x400)
    return bool(attributes & reparse)


def _inside(path: Path, root: Path, *, reject_reparse: bool = True) -> Optional[Path]:
    """Resolve *path* below *root*, rejecting path tricks and junctions."""

    root = root.resolve()
    raw = path if path.is_absolute() else root / path
    try:
        relative = raw.relative_to(root)
    except ValueError:
        return None
    if reject_reparse:
        current = root
        for part in relative.parts:
            current = current / part
            if _is_reparse_or_symlink(current):
                return None
    try:
        resolved = raw.resolve(strict=False)
        resolved.relative_to(root)
    except (OSError, ValueError):
        return None
    return resolved


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tree_snapshot(path: Path) -> tuple[str, int, str]:
    """Return ``(token, bytes, error)`` for one safe, link-free tree."""

    entries: list[tuple[str, str, int, int, str]] = []
    total = 0
    try:
        if _is_reparse_or_symlink(path) or not path.is_dir():
            return "", 0, "path is not a real directory"
        stack = [(path, "")]
        while stack:
            current, relative_prefix = stack.pop()
            with os.scandir(current) as iterator:
                children = sorted(iterator, key=lambda item: item.name)
            for entry in children:
                child = Path(entry.path)
                relative = f"{relative_prefix}/{entry.name}" if relative_prefix else entry.name
                if _is_reparse_or_symlink(child):
                    return "", 0, f"symlink or junction present at {relative}"
                info = entry.stat(follow_symlinks=False)
                if entry.is_dir(follow_symlinks=False):
                    entries.append((relative, "d", 0, info.st_mtime_ns, ""))
                    stack.append((child, relative))
                elif entry.is_file(follow_symlinks=False):
                    size = int(info.st_size)
                    total += size
                    # Metadata and lease contents are cheap to hash.  Payload
                    # bytes need only stable stat data because an active writer
                    # must hold the lease; hashing them would make previews of
                    # large scientific outputs unnecessarily expensive.
                    content_hash = ""
                    if entry.name in {INDEX_FILENAME, LEASE_FILENAME}:
                        content_hash = _hash_file(child)
                    entries.append((relative, "f", size, info.st_mtime_ns, content_hash))
                else:
                    return "", 0, f"unsupported filesystem entry at {relative}"
    except OSError as exc:
        return "", 0, f"unable to inspect tree: {exc}"
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest(), total, ""


def _try_lease(path: Path) -> bool:
    """Return true when a transaction lease is currently held elsewhere."""

    try:
        handle = path.open("a+b")
        if handle.seek(0, os.SEEK_END) == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        if os.name == "nt":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:  # pragma: no cover - exercised on POSIX CI
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        handle.close()
        return False
    except OSError:
        try:
            handle.close()
        except (UnboundLocalError, AttributeError, OSError):
            pass
        return True


def _load_parent_registry(storage: ArtifactStorageBackend) -> ArtifactRegistry:
    index = storage.base_dir / INDEX_FILENAME
    if _is_reparse_or_symlink(index) or not index.is_file():
        raise RegistryUnavailable(f"parent registry index is missing: {index}")
    registry = ArtifactRegistry.__new__(ArtifactRegistry)
    registry._mutex = threading.RLock()
    registry.storage = storage
    registry._index_path = index
    registry.lineage = LineageGraph()
    registry.registry = {}
    registry.task_commits = {}
    registry.audit_records = {}
    try:
        registry._load_persisted()
    except Exception as exc:
        raise RegistryUnavailable(f"parent registry is corrupt: {index}: {exc}") from exc
    return registry


def _validate_transaction_marker(path: Path, transaction_id: str) -> tuple[bool, str]:
    marker = path / LEASE_FILENAME
    if _is_reparse_or_symlink(marker):
        return False, "transaction lease marker is a symlink or junction"
    if not marker.exists():
        return False, "transaction has no lease marker"
    try:
        payload = json.loads(marker.read_text(encoding="utf-8"))
    except Exception as exc:
        return False, f"transaction lease marker is corrupt: {exc}"
    if not isinstance(payload, Mapping):
        return False, "transaction lease marker is not an object"
    if payload.get("lease_version") != 1 or payload.get("transaction_id") != transaction_id:
        return False, "transaction lease marker identity is invalid"
    return True, ""


def _parse_transaction_index(path: Path) -> tuple[dict[str, ArtifactMetadata], str]:
    index = path / INDEX_FILENAME
    if _is_reparse_or_symlink(index):
        return {}, "transaction registry index is a symlink or junction"
    if not index.exists():
        return {}, ""
    try:
        payload = json.loads(index.read_text(encoding="utf-8"))
        if not isinstance(payload, Mapping):
            return {}, "transaction registry root is not an object"
        records = payload.get("metadata", [])
        if not isinstance(records, list):
            return {}, "transaction registry metadata is not a list"
        metadata: dict[str, ArtifactMetadata] = {}
        for raw in records:
            item = ArtifactMetadata.model_validate(raw)
            canonical = ArtifactURI.parse(item.uri).to_string()
            if canonical != item.uri:
                return {}, f"transaction metadata URI is not canonical: {item.uri}"
            target = _inside(Path(item.storage_path), path)
            if target is None:
                return {}, f"transaction metadata path escapes transaction: {item.storage_path}"
            metadata[item.uri] = item
        return metadata, ""
    except Exception as exc:
        return {}, f"transaction registry index is corrupt: {exc}"


def _allowed_transaction_entry(relative: str, metadata: Mapping[str, ArtifactMetadata], root: Path) -> bool:
    if relative in {INDEX_FILENAME, LOCK_FILENAME, LEASE_FILENAME}:
        return True
    payload_paths = {Path(item.storage_path).relative_to(root).as_posix() for item in metadata.values()}
    if relative in payload_paths or any(value.startswith(relative + "/") for value in payload_paths):
        return True
    payload_names = {Path(value).name for value in payload_paths}
    name = Path(relative).name
    # A serializer may have been interrupted after creating a temporary file.
    # These are safe to remove with the abandoned transaction itself, but only
    # when their prefix identifies a registered payload or the registry index.
    if name.startswith(f".{INDEX_FILENAME}.") and name.endswith(".tmp"):
        return True
    return any(name.startswith(f".{payload_name}.") and name.endswith(".tmp") for payload_name in payload_names)


def _inspect_transaction(path: Path, now: datetime, base: Path) -> _TransactionInfo:
    transaction_id = path.name
    if not TRANSACTION_ID_RE.fullmatch(transaction_id):
        return _TransactionInfo(transaction_id, path, {}, "", 0, 0.0, False, "unknown transaction directory")
    if _is_reparse_or_symlink(path) or not path.is_dir():
        return _TransactionInfo(transaction_id, path, {}, "", 0, 0.0, False, "transaction path is a symlink or junction")
    marker_ok, marker_reason = _validate_transaction_marker(path, transaction_id)
    if not marker_ok:
        return _TransactionInfo(transaction_id, path, {}, "", 0, 0.0, False, marker_reason)
    metadata, index_reason = _parse_transaction_index(path)
    if index_reason:
        return _TransactionInfo(transaction_id, path, {}, "", 0, 0.0, False, index_reason)
    snapshot, size, tree_reason = _tree_snapshot(path)
    if tree_reason:
        return _TransactionInfo(transaction_id, path, metadata, snapshot, size, 0.0, False, tree_reason)
    try:
        newest = max([path.stat().st_mtime_ns, *[entry[3] for entry in _snapshot_entries(path)]]) / 1e9
        age_days = max(0.0, (now.timestamp() - newest) / 86400.0)
    except OSError as exc:
        return _TransactionInfo(transaction_id, path, metadata, snapshot, size, 0.0, False, f"unable to stat transaction: {exc}")
    # Ensure a transaction cannot smuggle an arbitrary file into the deletion
    # scope.  Recognized serializer leftovers are allowed above.
    try:
        for relative, *_ in _snapshot_entries(path):
            if not _allowed_transaction_entry(relative, metadata, path):
                return _TransactionInfo(transaction_id, path, metadata, snapshot, size, age_days, False, f"unknown transaction entry: {relative}")
    except OSError as exc:
        return _TransactionInfo(transaction_id, path, metadata, snapshot, size, age_days, False, f"unable to inspect transaction entries: {exc}")
    return _TransactionInfo(transaction_id, path, metadata, snapshot, size, age_days, True)


def _snapshot_entries(path: Path) -> list[tuple[str, str, int, int, str]]:
    """Return snapshot records and raise on links/special files."""
    entries: list[tuple[str, str, int, int, str]] = []
    stack = [(path, "")]
    while stack:
        current, prefix = stack.pop()
        with os.scandir(current) as iterator:
            children = sorted(iterator, key=lambda item: item.name)
        for entry in children:
            child = Path(entry.path)
            relative = f"{prefix}/{entry.name}" if prefix else entry.name
            if _is_reparse_or_symlink(child):
                raise OSError(f"symlink or junction present at {relative}")
            info = entry.stat(follow_symlinks=False)
            if entry.is_dir(follow_symlinks=False):
                entries.append((relative, "d", 0, info.st_mtime_ns, ""))
                stack.append((child, relative))
            elif entry.is_file(follow_symlinks=False):
                digest = _hash_file(child) if entry.name in {INDEX_FILENAME, LEASE_FILENAME} else ""
                entries.append((relative, "f", int(info.st_size), info.st_mtime_ns, digest))
            else:
                raise OSError(f"unsupported filesystem entry at {relative}")
    return entries


def _candidate_refs(
    registry: ArtifactRegistry,
    infos: Mapping[str, _TransactionInfo],
    base: Path,
    *,
    journal_values: Sequence[tuple[Path, Any]] = (),
) -> tuple[dict[str, set[str]], list[str]]:
    """Map durable references to transaction IDs and report unknown references."""

    refs: dict[str, set[str]] = {transaction_id: set() for transaction_id in infos}
    reasons: list[str] = []
    uri_map: dict[str, set[str]] = {}
    for transaction_id, info in infos.items():
        for uri in info.metadata:
            uri_map.setdefault(uri, set()).add(transaction_id)

    def add_path(raw_path: Any, source: str) -> None:
        if not isinstance(raw_path, (str, os.PathLike)):
            reasons.append(f"{source} contains a non-path reference")
            return
        raw = Path(raw_path)
        target = _inside(raw, base)
        if target is None:
            reasons.append(f"{source} path escapes artifact storage: {raw_path}")
            return
        for transaction_id, info in infos.items():
            if target == info.path or info.path in target.parents:
                refs[transaction_id].add(source)

    def add_uri(raw_uri: Any, source: str, unknown_is_error: bool = True) -> None:
        if not isinstance(raw_uri, str):
            reasons.append(f"{source} contains a non-string artifact URI")
            return
        try:
            uri = ArtifactURI.parse(raw_uri).to_string()
        except Exception:
            if unknown_is_error:
                reasons.append(f"{source} contains an invalid artifact URI: {raw_uri!r}")
            return
        if uri in registry.registry:
            add_path(registry.registry[uri].storage_path, source)
            return
        if uri in uri_map:
            for transaction_id in uri_map[uri]:
                refs[transaction_id].add(source)
            return
        if unknown_is_error:
            reasons.append(f"{source} references unknown artifact URI: {uri}")

    # The parent metadata is authoritative for every published payload.
    for uri, metadata in registry.registry.items():
        add_path(metadata.storage_path, f"metadata:{uri}")

    # Computation receipts and audit receipts carry URI references even if a
    # metadata record was accidentally pruned.  Such ambiguity blocks all
    # deletion rather than guessing that the payload is disposable.
    for signature, receipt in registry.task_commits.items():
        source = f"receipt:{signature}"
        if not isinstance(receipt, Mapping):
            reasons.append(f"{source} is not an object")
            continue
        hashes = receipt.get("output_hashes", {})
        if not isinstance(hashes, Mapping):
            reasons.append(f"{source}.output_hashes is not an object")
        else:
            for raw_uri in hashes:
                add_uri(raw_uri, f"{source}.output_hashes")
        result = receipt.get("result", {})
        if isinstance(result, Mapping):
            for key in ("input_artifacts", "output_artifacts"):
                values = result.get(key, [])
                if not isinstance(values, list):
                    reasons.append(f"{source}.result.{key} is not a list")
                    continue
                for raw_uri in values:
                    add_uri(raw_uri, f"{source}.result.{key}")
        elif result is not None:
            reasons.append(f"{source}.result is not an object")

    for signature, record in registry.audit_records.items():
        source = f"audit:{signature}"
        for field_name in ("input_artifacts", "output_artifacts"):
            values = getattr(record, field_name, [])
            if not isinstance(values, list):
                reasons.append(f"{source}.{field_name} is not a list")
                continue
            for raw_uri in values:
                add_uri(raw_uri, f"{source}.{field_name}")
        hashes = getattr(record, "artifact_hashes", {})
        if isinstance(hashes, Mapping):
            for raw_uri in hashes:
                add_uri(raw_uri, f"{source}.artifact_hashes")
        else:
            reasons.append(f"{source}.artifact_hashes is not an object")

    # Journal entries are intentionally read only from the well-known _runs
    # directory.  We inspect artifact-looking fields and direct storage paths;
    # log strings are not treated as durable references.
    def walk(value: Any, source: str, key: str = "") -> None:
        if isinstance(value, Mapping):
            for child_key, child in value.items():
                child_name = str(child_key)
                walk(child, f"{source}.{child_name}", child_name.lower())
            return
        if isinstance(value, list):
            for position, child in enumerate(value):
                walk(child, f"{source}[{position}]", key)
            return
        if not key:
            return
        if "storage_path" in key or key in {"path", "payload_path"}:
            add_path(value, source)
        elif "artifact" in key or key.endswith("_uri") or key == "uri":
            if isinstance(value, list):
                for item in value:
                    add_uri(item, source)
            else:
                add_uri(value, source)

    for journal_path, value in journal_values:
        walk(value, f"journal:{journal_path.name}")
    return refs, reasons


def _journal_snapshot(base: Path) -> tuple[list[tuple[Path, Any]], list[tuple[str, int, int, str]], list[str]]:
    runs = base / "_runs"
    if not runs.exists():
        return [], [], []
    if _is_reparse_or_symlink(runs) or not runs.is_dir():
        return [], [], ["study journal directory is a symlink or not a directory"]
    values: list[tuple[Path, Any]] = []
    snapshot: list[tuple[str, int, int, str]] = []
    reasons: list[str] = []
    try:
        for path in sorted(runs.iterdir(), key=lambda item: item.name):
            if _is_reparse_or_symlink(path):
                reasons.append(f"study journal entry is a symlink or junction: {path.name}")
                continue
            if not path.is_file() or path.suffix.lower() != ".json":
                continue
            try:
                raw = path.read_bytes()
                values.append((path, json.loads(raw.decode("utf-8"))))
                info = path.stat()
                snapshot.append((path.name, int(info.st_size), info.st_mtime_ns, hashlib.sha256(raw).hexdigest()))
            except Exception as exc:
                reasons.append(f"study journal is corrupt: {path}: {exc}")
    except OSError as exc:
        reasons.append(f"unable to inspect study journals: {exc}")
    return values, snapshot, reasons


def _build_plan_locked(storage: ArtifactStorageBackend, older_than_days: int, now: datetime) -> CleanupPlan:
    base = storage.base_dir
    generated_at = now
    if not isinstance(older_than_days, int) or older_than_days < 0:
        raise ValueError("older_than_days must be a non-negative integer")
    index = base / INDEX_FILENAME
    registry_hash = ""
    try:
        registry_hash = _hash_file(index)
    except OSError:
        pass

    blocked_reasons: list[str] = []
    try:
        registry = _load_parent_registry(storage)
        registry_available = True
    except RegistryUnavailable as exc:
        registry = None
        registry_available = False
        blocked_reasons.append(str(exc))

    transactions_root = base / "_transactions"
    infos: dict[str, _TransactionInfo] = {}
    if transactions_root.exists():
        if _is_reparse_or_symlink(transactions_root) or not transactions_root.is_dir():
            blocked_reasons.append("transaction root is a symlink, junction, or not a directory")
        else:
            try:
                for path in sorted(transactions_root.iterdir(), key=lambda item: item.name):
                    if not path.is_dir() and not _is_reparse_or_symlink(path):
                        continue
                    info = _inspect_transaction(path, now, base)
                    infos[info.transaction_id] = info
            except OSError as exc:
                blocked_reasons.append(f"unable to inspect transaction root: {exc}")

    journal_values, journal_snapshot, journal_reasons = _journal_snapshot(base)
    blocked_reasons.extend(journal_reasons)
    references: dict[str, set[str]] = {transaction_id: set() for transaction_id in infos}
    if registry is not None:
        refs, ref_reasons = _candidate_refs(registry, infos, base, journal_values=journal_values)
        references.update(refs)
        blocked_reasons.extend(ref_reasons)

    candidates: list[CleanupCandidate] = []
    for transaction_id, info in infos.items():
        reasons: list[str] = []
        if not info.valid:
            reasons.append(info.reason)
        if not registry_available:
            reasons.append("parent registry unavailable; fail closed")
        if references.get(transaction_id):
            reasons.append("live reference: " + ", ".join(sorted(references[transaction_id])[:4]))
        if info.valid and _try_lease(info.path / LEASE_FILENAME):
            reasons.append("transaction lease is active")
        if info.valid and info.age_days < older_than_days:
            reasons.append(f"younger than retention ({older_than_days} days)")
        # A global ambiguity blocks application even for candidates without a
        # direct reference; preserving every candidate is the fail-closed rule.
        if blocked_reasons:
            reasons.append("maintenance is blocked by registry/reference validation")
        eligible = info.valid and not reasons and not blocked_reasons
        reason = "; ".join(reasons) if reasons else "eligible abandoned transaction"
        candidates.append(
            CleanupCandidate(
                transaction_id=transaction_id,
                path=info.path,
                size_bytes=info.size_bytes,
                age_days=info.age_days,
                eligible=eligible,
                reason=reason,
                snapshot=info.snapshot,
                references=tuple(sorted(references.get(transaction_id, set()))),
            )
        )

    state = {
        "registry": registry_hash,
        "journals": journal_snapshot,
        "transactions": [(item.transaction_id, item.snapshot, item.size_bytes, item.reason) for item in candidates],
        "blocked": blocked_reasons,
        "older_than_days": older_than_days,
    }
    token = hashlib.sha256(json.dumps(state, sort_keys=True, default=str, separators=(",", ":")).encode("utf-8")).hexdigest()
    return CleanupPlan(
        base_dir=base,
        older_than_days=older_than_days,
        generated_at=generated_at,
        candidates=tuple(candidates),
        blocked_reasons=tuple(dict.fromkeys(blocked_reasons)),
        snapshot_token=token,
        registry_available=registry_available,
    )


def plan_cleanup(
    storage_dir: str | Path | ArtifactRegistry,
    older_than_days: int = DEFAULT_RETENTION_DAYS,
    *,
    now: Optional[datetime] = None,
) -> CleanupPlan:
    """Build a conservative cleanup plan without deleting anything."""

    base = _normalise_base(storage_dir)
    if not base.is_dir() or not (base / INDEX_FILENAME).is_file():
        raise ArtifactMaintenanceError("Cleanup requires an existing artifact root and registry index")
    storage = ArtifactStorageBackend(str(base))
    with storage.lock():
        return _build_plan_locked(storage, older_than_days, _utc_now(now))


def _report_from_plan(plan: CleanupPlan, *, mode: str = "preview") -> CleanupReport:
    reasons: dict[str, int] = {}
    skipped: list[dict[str, Any]] = []
    for candidate in plan.candidates:
        if not candidate.eligible:
            reasons[candidate.reason] = reasons.get(candidate.reason, 0) + 1
            skipped.append(candidate.as_dict())
    for reason in plan.blocked_reasons:
        reasons[reason] = reasons.get(reason, 0) + 1
    return CleanupReport(
        mode=mode,
        plan=plan,
        planned_count=plan.planned_count,
        planned_bytes=plan.planned_bytes,
        skipped=skipped,
        reasons=reasons,
    )


def _remove_tree(path: Path) -> None:
    """Delete a previously validated tree without following links."""

    if _is_reparse_or_symlink(path) or not path.is_dir():
        raise ArtifactMaintenanceError(f"refusing to remove unsafe transaction path: {path}")
    # Validate immediately before mutating; the lease check and parent lock
    # prevent legitimate transaction writers from changing this tree.
    _, _, error = _tree_snapshot(path)
    if error:
        raise ArtifactMaintenanceError(f"refusing to remove unsafe transaction tree: {error}")
    for entry in sorted(path.iterdir(), key=lambda item: item.name, reverse=True):
        if _is_reparse_or_symlink(entry):
            raise ArtifactMaintenanceError(f"refusing to remove link entry: {entry}")
        if entry.is_dir():
            _remove_tree(entry)
        else:
            entry.unlink()
    path.rmdir()


def cleanup_artifacts(
    storage_dir: str | Path | ArtifactRegistry,
    *,
    apply: bool = False,
    older_than_days: int = DEFAULT_RETENTION_DAYS,
    plan: Optional[CleanupPlan] = None,
    now: Optional[datetime] = None,
) -> CleanupReport:
    """Preview or apply cleanup of old, abandoned transaction directories.

    ``apply=False`` is the default.  When *plan* is supplied with
    ``apply=True``, the plan is revalidated under the storage lock and becomes
    a no-op if any registry, journal, or transaction state changed.
    """

    base = _normalise_base(storage_dir)
    if plan is not None and plan.base_dir.resolve() != base:
        raise ValueError("cleanup plan belongs to a different storage directory")
    if not apply:
        selected = plan if plan is not None else plan_cleanup(base, older_than_days, now=now)
        return _report_from_plan(selected, mode="preview")

    if not base.is_dir() or not (base / INDEX_FILENAME).is_file():
        raise ArtifactMaintenanceError("Cleanup requires an existing artifact root and registry index")
    storage = ArtifactStorageBackend(str(base))
    with storage.lock():
        current = _build_plan_locked(storage, older_than_days if plan is None else plan.older_than_days, _utc_now(now))
        if plan is not None and current.snapshot_token != plan.snapshot_token:
            report = _report_from_plan(current, mode="apply")
            report.stale_plan = True
            report.planned_count = plan.planned_count
            report.planned_bytes = plan.planned_bytes
            report.reasons["stale_plan"] = report.reasons.get("stale_plan", 0) + 1
            report.errors.append("cleanup preview is stale; no directories were removed")
            return report

        report = _report_from_plan(current, mode="apply")
        for candidate in current.eligible:
            try:
                safe_path = _inside(candidate.path, base)
                if safe_path is None or safe_path.parent != base / "_transactions" or not TRANSACTION_ID_RE.fullmatch(safe_path.name):
                    raise ArtifactMaintenanceError("Cleanup target escaped the transaction root")
                if _try_lease(candidate.path / LEASE_FILENAME):
                    report.skipped.append({**candidate.as_dict(), "reason": "transaction lease is active"})
                    report.reasons["transaction lease is active"] = report.reasons.get("transaction lease is active", 0) + 1
                    continue
                latest_snapshot, latest_size, snapshot_error = _tree_snapshot(candidate.path)
                if snapshot_error or latest_snapshot != candidate.snapshot or latest_size != candidate.size_bytes:
                    report.skipped.append({**candidate.as_dict(), "reason": "transaction changed after revalidation"})
                    report.reasons["transaction changed after revalidation"] = report.reasons.get("transaction changed after revalidation", 0) + 1
                    continue
                _remove_tree(candidate.path)
                report.deleted_count += 1
                report.deleted_bytes += candidate.size_bytes
            except Exception as exc:
                report.skipped.append({**candidate.as_dict(), "reason": f"delete failed: {exc}"})
                report.reasons["delete failed"] = report.reasons.get("delete failed", 0) + 1
                report.errors.append(f"{candidate.path}: {exc}")
        return report


def apply_cleanup(
    storage_dir: str | Path | ArtifactRegistry,
    plan: CleanupPlan,
    *,
    now: Optional[datetime] = None,
) -> CleanupReport:
    """Apply a previously generated plan after revalidation."""

    return cleanup_artifacts(storage_dir, apply=True, plan=plan, now=now)


class ArtifactMaintenance:
    """Small object wrapper for integrations that prefer an instance API."""

    def __init__(self, storage_dir: str | Path | ArtifactRegistry):
        self.storage_dir = _normalise_base(storage_dir)

    def plan(self, older_than_days: int = DEFAULT_RETENTION_DAYS, *, now: Optional[datetime] = None) -> CleanupPlan:
        return plan_cleanup(self.storage_dir, older_than_days, now=now)

    def cleanup(
        self,
        *,
        apply: bool = False,
        older_than_days: int = DEFAULT_RETENTION_DAYS,
        plan: Optional[CleanupPlan] = None,
        now: Optional[datetime] = None,
    ) -> CleanupReport:
        return cleanup_artifacts(
            self.storage_dir,
            apply=apply,
            older_than_days=older_than_days,
            plan=plan,
            now=now,
        )


# Compatibility spellings make the maintenance boundary easy to discover for
# callers while keeping one implementation and one safety policy.
ArtifactCleanup = ArtifactMaintenance
build_cleanup_plan = plan_cleanup
preview_cleanup = plan_cleanup


__all__ = [
    "ArtifactCleanup",
    "ArtifactMaintenance",
    "ArtifactMaintenanceError",
    "CleanupCandidate",
    "CleanupPlan",
    "CleanupReport",
    "DEFAULT_RETENTION_DAYS",
    "apply_cleanup",
    "build_cleanup_plan",
    "cleanup_artifacts",
    "plan_cleanup",
    "preview_cleanup",
]
