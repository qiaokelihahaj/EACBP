"""Versioned, integrity-checked research bundles for offline relocation.

A bundle preserves computation/audit receipts and snapshots. It does not grant
new admission or promise that a changed host can resume the computation.
"""
from __future__ import annotations

from contextlib import ExitStack, contextmanager
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import zipfile

from eacbp._atomic_json import atomic_write_json
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.artifact.storage import _sha256_file
from eacbp.artifact.transaction import _lock_exclusive_nonblocking, _unlock
from eacbp.evidence.snapshot import (
    load_study_snapshot, write_study_snapshot, _verify_registry_scope,
)

BUNDLE_FORMAT = "eacbp.research_bundle"
BUNDLE_VERSION = 1
MANIFEST_NAME = "bundle-manifest.json"


class ArtifactBundleError(ValueError):
    """A bundle cannot be published/imported without losing integrity."""


def _key(value: str) -> str:
    if (not isinstance(value, str) or not value or value == "."
            or any(character in value for character in '\\:<>|?*')
            or any(ord(character) < 32 for character in value)):
        raise ArtifactBundleError(f"invalid bundle object key: {value!r}")
    path = PurePosixPath(value)
    if path.is_absolute() or value != path.as_posix() or any(part in (".", "..") for part in path.parts):
        raise ArtifactBundleError(f"unsafe bundle object key: {value!r}")
    for part in path.parts:
        if part.endswith((".", " ")) or part.split(".")[0].upper() in {
            "CON", "PRN", "AUX", "NUL", *(f"COM{i}" for i in range(1, 10)),
            *(f"LPT{i}" for i in range(1, 10)),
        }:
            raise ArtifactBundleError(f"nonportable bundle object key: {value!r}")
    return value


def _regular_file(path: Path) -> None:
    # Check parents as well: resolving first would conceal directory symlinks.
    for item in (path, *path.parents):
        info = item.lstat()
        if item.is_symlink() or getattr(info, "st_file_attributes", 0) & 0x400:
            raise ArtifactBundleError(f"symlink/junction is not a bundle source: {item}")
    if not path.is_file():
        raise ArtifactBundleError(f"bundle source is not a regular file: {path}")


@contextmanager
def _study_lock(path):
    _regular_file(path)
    with path.open("r+b") as handle:
        try:
            _lock_exclusive_nonblocking(handle)
        except OSError as exc:
            raise ArtifactBundleError("Cannot export a study while it is running") from exc
        try:
            yield
        finally:
            _unlock(handle)


def _validate_registry(registry):
    metadata = registry.verify_many([item.uri for item in registry.list_artifacts()])
    for signature, receipt in registry.task_commits.items():
        if not isinstance(receipt, dict) or receipt.get("signature") != signature:
            raise ArtifactBundleError("invalid computation receipt")
        hashes = receipt.get("output_hashes")
        if not isinstance(hashes, dict):
            raise ArtifactBundleError("invalid computation output hashes")
        for uri, digest in hashes.items():
            if uri not in metadata or metadata[uri].sha256_hash != digest:
                raise ArtifactBundleError(f"computation receipt does not match {uri}")
    for record in registry.audit_records.values():
        receipt = registry.task_commits.get(record.signature)
        if receipt is None or receipt.get("output_hashes") != record.artifact_hashes:
            raise ArtifactBundleError("audit and computation receipts disagree")
        if record.status.value == "passed":
            # Reuse the existing gate, including context/report and sibling
            # deserialization checks. No new receipt is produced.
            registry._require_audited_locked(record.signature, record.contract)


def _validate_tree(root, snapshots):
    registry = ArtifactRegistry(str(root / "artifacts"))
    _validate_registry(registry)
    for key in snapshots:
        snapshot = load_study_snapshot(root / _key(key))
        _verify_registry_scope(snapshot, registry)
    return registry


def _validate_study_envelope(root, snapshots, files):
    """Validate the application study contract before publishing a run tree."""
    if "snapshot.json" not in snapshots or "run_config.json" not in files:
        raise ArtifactBundleError("Study bundle requires snapshot.json and run_config.json")
    from eacbp.schemas.study import StudyManifest

    try:
        envelope = json.loads((root / "run_config.json").read_text(encoding="utf-8"))
        if (not isinstance(envelope, dict) or type(envelope.get("schema_version")) is not int
                or envelope["schema_version"] != 1 or envelope.get("import_completed") is not True
                or not isinstance(envelope.get("manifest"), dict)
                or not isinstance(envelope.get("config"), dict)):
            raise ArtifactBundleError("Invalid or incomplete study run_config envelope")
        manifest = StudyManifest.model_validate(envelope["manifest"])
        snapshot = load_study_snapshot(root / "snapshot.json")
        if (envelope.get("study_id") != snapshot.manifest.study_id
                or manifest.model_dump(mode="json") != snapshot.manifest.model_dump(mode="json")):
            raise ArtifactBundleError("Study run_config manifest does not match snapshot.json")
    except (OSError, TypeError, ValueError) as exc:
        if isinstance(exc, ArtifactBundleError):
            raise
        raise ArtifactBundleError(f"Invalid study run_config envelope: {exc}") from exc


def export_artifact_bundle(storage_dir, bundle_path, *, snapshot_paths=(), run_files=None,
                           require_study=False):
    """Export published artifacts plus specified snapshots into a new ZIP.

    ``snapshot_paths`` accepts paths (stored under snapshots/<basename>) or a
    mapping of relative bundle names to paths (e.g. snapshot.json). ``run_files``
    optionally includes report/config files under explicit relative names.
    Existing destinations are never replaced. Private/uncommitted payloads and
    lease files are excluded. Study journals/event files are preserved.
    ``require_study`` validates the same application envelope required on study
    import before a completed ZIP is published.
    """
    source = Path(storage_dir).expanduser().absolute()
    _regular_file(source / ArtifactRegistry.INDEX_FILENAME)
    registry = ArtifactRegistry(str(source))
    target = Path(bundle_path).expanduser().absolute()
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(snapshot_paths, "items"):
        snapshot_map = dict(snapshot_paths)
    else:
        paths = list(snapshot_paths)
        snapshot_map = {f"snapshots/{Path(path).name}": path for path in paths}
        if len(snapshot_map) != len(paths):
            raise ArtifactBundleError("snapshot filenames must be unique")
    extras = dict(run_files or {})
    names = list(snapshot_map) + list(extras)
    if len({key.casefold() for key in names}) != len(names):
        raise ArtifactBundleError("duplicate bundle paths")
    for key in names:
        _key(key)
        if key.casefold() == MANIFEST_NAME or key.casefold() == "artifacts" or key.casefold().startswith("artifacts/"):
            raise ArtifactBundleError("reserved bundle path")
    with tempfile.TemporaryDirectory(prefix=".eacbp-export-", dir=target.parent) as temporary:
        workspace = Path(temporary)
        tree = workspace / "tree"
        artifact_root = tree / "artifacts"
        artifact_root.mkdir(parents=True)
        with ExitStack() as stack:
            for lock in sorted((source / "_runs").glob("*.lock")):
                stack.enter_context(_study_lock(lock))
            with registry.storage.lock():
                registry._load_persisted()
                payload = registry._index_payload()
                for item in registry.registry.values():
                    original = Path(item.storage_path)
                    _regular_file(original)
                    registry.storage.verify_payload(item.uri, item.type, item.sha256_hash,
                                                    item.storage_path, expected_size=item.size_bytes)
                    key = _key(registry.storage.relative_path(original))
                    copied = artifact_root / key
                    copied.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(original, copied)
                atomic_write_json(artifact_root / registry.INDEX_FILENAME, payload)
                for original in sorted((source / "_runs").rglob("*")):
                    if original.is_file() and original.suffix in {".json", ".jsonl"}:
                        _regular_file(original)
                        copied = artifact_root / original.relative_to(source)
                        copied.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(original, copied)
            for key, original in snapshot_map.items():
                _regular_file(Path(original).absolute())
                snapshot = load_study_snapshot(original)
                _verify_registry_scope(snapshot, registry)
                write_study_snapshot(
                    tree / key, manifest=snapshot.manifest, config=snapshot.config,
                    summary=snapshot.summary, evidence_graph=snapshot.evidence_graph,
                    artifact_registry=registry, task_history=snapshot.task_results,
                    audit_reports=snapshot.audit_reports,
                )
            for key, original in extras.items():
                original = Path(original).absolute()
                _regular_file(original)
                copied = tree / key
                copied.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(original, copied)
        _validate_tree(tree, snapshot_map)
        files = {}
        for path in sorted(tree.rglob("*")):
            if path.is_file() and path.name != ".registry.lock":
                digest, size = _sha256_file(path)
                files[path.relative_to(tree).as_posix()] = {"sha256": digest, "size_bytes": size}
        manifest = {"format": BUNDLE_FORMAT, "version": BUNDLE_VERSION,
                    "files": files, "snapshots": sorted(snapshot_map)}
        if require_study:
            _validate_study_envelope(tree, snapshot_map, files)
        archive_path = workspace / "bundle.zip"
        with zipfile.ZipFile(archive_path, "w", compression=zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
            archive.writestr(MANIFEST_NAME, json.dumps(manifest, sort_keys=True))
            for key in files:
                archive.write(tree / key, key)
        # Atomic create-if-absent, matching artifact publication semantics.
        os.link(archive_path, target)
    return target


def import_artifact_bundle(bundle_path, destination, *, require_study=False):
    """Validate/extract into a new directory; failures publish nothing.

    ``require_study`` additionally requires the latest snapshot and a complete,
    matching application run-config envelope before the target is published.
    """
    target = Path(destination).expanduser().absolute()
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".eacbp-import-", dir=target.parent) as temporary:
        tree = Path(temporary) / "run"
        tree.mkdir()
        with zipfile.ZipFile(bundle_path) as archive:
            infos = archive.infolist()
            names = [item.filename for item in infos]
            if len({name.casefold() for name in names}) != len(names):
                raise ArtifactBundleError("duplicate archive members")
            for item in infos:
                _key(item.filename)
                file_type = stat.S_IFMT(item.external_attr >> 16)
                if item.is_dir() or file_type not in (0, stat.S_IFREG):
                    raise ArtifactBundleError("archive entries must be regular files")
            if MANIFEST_NAME not in names:
                raise ArtifactBundleError("bundle manifest is missing")
            manifest = json.loads(archive.read(MANIFEST_NAME))
            if (not isinstance(manifest, dict) or manifest.get("format") != BUNDLE_FORMAT or type(manifest.get("version")) is not int
                    or manifest["version"] != BUNDLE_VERSION):
                raise ArtifactBundleError("unsupported research bundle version")
            files, snapshots = manifest.get("files"), manifest.get("snapshots")
            if not isinstance(files, dict) or not isinstance(snapshots, list):
                raise ArtifactBundleError("invalid bundle manifest")
            if (MANIFEST_NAME in files or any(not isinstance(key, str) for key in snapshots)
                    or len(set(snapshots)) != len(snapshots)):
                raise ArtifactBundleError("invalid bundle file/snapshot manifest")
            if set(names) != set(files) | {MANIFEST_NAME}:
                raise ArtifactBundleError("archive does not match manifest file set")
            for key, expected in files.items():
                path = tree / _key(key)
                if (not isinstance(expected, dict) or type(expected.get("size_bytes")) is not int
                        or expected["size_bytes"] < 0
                        or archive.getinfo(key).file_size != expected["size_bytes"]):
                    raise ArtifactBundleError(f"invalid bundle file size: {key}")
                path.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(key) as source, path.open("xb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                digest, size = _sha256_file(path)
                if not isinstance(expected, dict) or expected != {"sha256": digest, "size_bytes": size}:
                    raise ArtifactBundleError(f"bundle content hash mismatch: {key}")
            if any(key not in files for key in snapshots):
                raise ArtifactBundleError("snapshot missing from bundle")
        index_path = tree / "artifacts" / ArtifactRegistry.INDEX_FILENAME
        index = json.loads(index_path.read_text(encoding="utf-8"))
        if index.get("index_version") != ArtifactRegistry.INDEX_VERSION:
            raise ArtifactBundleError("bundle requires a portable registry index")
        for item in index.get("metadata", []):
            key = _key(item["storage_path"])
            if f"artifacts/{key}" not in files:
                raise ArtifactBundleError("registered artifact missing from bundle")
        _validate_tree(tree, snapshots)
        if require_study:
            _validate_study_envelope(tree, snapshots, files)
        # Windows rename rejects an existing destination. POSIX can replace an
        # empty directory, so reserve that destination exclusively first.
        if os.name == "nt":
            os.rename(tree, target)
        else:
            target.mkdir()
            try:
                os.replace(tree, target)
            except BaseException:
                target.rmdir()
                raise
    return target


__all__ = ["ArtifactBundleError", "export_artifact_bundle", "import_artifact_bundle"]
