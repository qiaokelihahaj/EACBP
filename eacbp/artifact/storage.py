"""Durable, immutable storage for EACBP artifacts.

The storage layer deliberately separates three concerns that used to be mixed
together: serialization, atomic publication, and integrity verification.  A
payload is always serialized to a temporary file in the destination directory,
hashed, and published with an exclusive hard-link.  Consequently a competing
writer can never replace an existing artifact and a failed serializer cannot
leave a partially written target behind.
"""

from __future__ import annotations

import hashlib
import json
import os
import base64
import tempfile
import threading
from contextlib import contextmanager
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, Iterator, Mapping, Optional, Tuple

import numpy as np
import pandas as pd

from eacbp.artifact.uri import ArtifactURI
from eacbp.schemas.artifact import ArtifactType


class ArtifactStorageError(Exception):
    """Base class for storage and serialization failures."""


class ArtifactAlreadyExistsError(ArtifactStorageError):
    """Raised when attempting to overwrite an immutable artifact version."""


class ArtifactIntegrityError(ValueError, ArtifactStorageError):
    """Raised when a persisted artifact no longer matches its registered hash."""


class ArtifactSerializationError(ArtifactStorageError):
    """Raised when a payload cannot be represented by its declared format."""


class ArtifactDependencyError(ArtifactSerializationError):
    """Raised when an optional serializer dependency is unavailable."""


def _sha256_file(path: Path) -> Tuple[str, int]:
    """Return a streaming SHA-256 digest and byte count for *path*."""

    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
            size += len(chunk)
    return f"sha256:{digest.hexdigest()}", size


def _json_default(value: Any) -> Any:
    """Convert common scientific Python scalars without stringifying structure.

    The old serializer used ``default=str``.  That made nested arrays, dates,
    and pandas values indistinguishable strings on a round trip.  This helper
    performs only explicit, loss-aware conversions and raises for unknown
    values so callers see a real serialization failure.
    """

    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (pd.Timestamp, datetime, date)):
        return value.isoformat()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (set, frozenset)):
        return list(value)
    if isinstance(value, pd.DataFrame):
        return {"__eacbp_dataframe__": value.to_dict(orient="split")}
    if isinstance(value, pd.Series):
        return {"__eacbp_series__": value.to_dict()}
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


class PayloadSerializer:
    """Serialize and deserialize the supported artifact payload types."""

    _ANNDATA_TYPES = (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA)

    @staticmethod
    def _get_anndata_module():
        try:
            import anndata as ad
        except ImportError as exc:  # pragma: no cover - depends on optional bio extra
            raise ArtifactDependencyError(
                "AnnData artifacts require the optional 'anndata' dependency; "
                "use the explicit .npz fallback by passing SCData or a dictionary payload"
            ) from exc
        return ad

    @classmethod
    def preferred_format(cls, payload: Any, artifact_type: ArtifactType) -> str:
        """Return ``h5ad`` or the explicit structured ``npz`` fallback.

        A fallback is selected before writing, so an NPZ is never hidden inside
        a file whose suffix says ``.h5ad``.  If AnnData is available, payloads
        that can be converted to AnnData use h5ad; plain dictionaries without
        the required matrix/table keys use the structured fallback.
        """

        if artifact_type not in cls._ANNDATA_TYPES:
            return "native"
        try:
            ad = cls._get_anndata_module()
        except ArtifactDependencyError:
            if isinstance(payload, dict) or payload.__class__.__name__ == "SCData":
                return "npz"
            raise

        if isinstance(payload, ad.AnnData):
            return "h5ad"
        if payload.__class__.__name__ == "SCData":
            return "h5ad"
        if hasattr(payload, "to_anndata"):
            return "h5ad"
        if isinstance(payload, dict) and {"X", "obs", "var"}.issubset(payload):
            return "h5ad"
        # A dictionary such as {"X": ...} is still useful in lightweight
        # environments; it is stored with the explicit structured fallback.
        return "npz"

    @staticmethod
    def _as_anndata(payload: Any) -> Any:
        ad = PayloadSerializer._get_anndata_module()
        if isinstance(payload, ad.AnnData):
            return payload
        if payload.__class__.__name__ == "SCData":
            converted = payload.to_anndata()
            if not isinstance(converted, ad.AnnData):
                raise ArtifactSerializationError(
                    "SCData.to_anndata() did not return an anndata.AnnData instance"
                )
            return converted
        if hasattr(payload, "to_anndata"):
            converted = payload.to_anndata()
            if not isinstance(converted, ad.AnnData):
                raise ArtifactSerializationError(
                    "to_anndata() did not return an anndata.AnnData instance"
                )
            return converted
        if isinstance(payload, dict) and {"X", "obs", "var"}.issubset(payload):
            from eacbp.capabilities.sc_data import SCData

            return SCData.from_dict(payload).to_anndata()
        raise ArtifactSerializationError(
            f"Cannot convert {type(payload).__name__} to AnnData; "
            "provide an AnnData/SCData payload or a complete dictionary"
        )

    @staticmethod
    def _fallback_payload(payload: Any) -> Any:
        # SCData is intentionally represented by its dictionary contract so the
        # lightweight fallback has the same public shape as older callers.
        if payload.__class__.__name__ == "SCData":
            return payload.to_dict()
        return payload

    @classmethod
    def _serialize_structured_npz(cls, payload: Any, target_path: Path) -> None:
        """Write a restricted typed tree plus numeric arrays to an NPZ.

        The fallback is intentionally independent of Python object deserialization
        (``allow_pickle=False`` on read).  Dictionaries, lists, pandas frames,
        sparse matrices, nested ``uns`` values, and primitive scientific scalars
        are represented by explicit tags; unsupported objects fail loudly.
        """

        arrays: Dict[str, np.ndarray] = {}
        counter = [0]

        def array_ref(value: np.ndarray) -> str:
            if value.dtype.hasobject:
                raise ArtifactSerializationError(
                    "Object-dtype arrays require an explicit element representation"
                )
            key = f"array_{counter[0]}"
            counter[0] += 1
            arrays[key] = np.array(value, copy=True)
            return key

        def encode(value: Any) -> Dict[str, Any]:
            if value is None:
                return {"t": "none"}
            if isinstance(value, (bool, int, str)):
                return {"t": "scalar", "v": value}
            if isinstance(value, float):
                return {"t": "scalar", "v": value}
            if isinstance(value, complex):
                return {"t": "complex", "real": value.real, "imag": value.imag}
            if isinstance(value, np.generic):
                if np.issubdtype(value.dtype, np.complexfloating):
                    item = complex(value.item())
                    return {"t": "complex", "real": item.real, "imag": item.imag}
                if np.issubdtype(value.dtype, np.datetime64):
                    return {"t": "datetime64", "v": str(value)}
                return {"t": "scalar", "v": value.item(), "dtype": str(value.dtype)}
            if isinstance(value, np.ndarray):
                if value.dtype.hasobject:
                    flat = [encode(item) for item in value.reshape(-1).tolist()]
                    return {
                        "t": "ndarray_object",
                        "shape": list(value.shape),
                        "items": flat,
                    }
                return {
                    "t": "ndarray",
                    "key": array_ref(value),
                    "shape": list(value.shape),
                    "dtype": str(value.dtype),
                }
            if isinstance(value, pd.Categorical):
                return {
                    "t": "categorical",
                    "categories": encode(value.categories.tolist()),
                    "ordered": bool(value.ordered),
                    "codes": encode(np.asarray(value.codes, dtype=np.int64)),
                }
            if isinstance(value, pd.DataFrame):
                columns = [encode(column) for column in value.columns.tolist()]
                encoded_columns = []
                for column in value.columns:
                    series = value[column]
                    if isinstance(series.dtype, pd.CategoricalDtype):
                        encoded_columns.append({
                            "categorical": True,
                            "categories": encode(series.cat.categories.tolist()),
                            "ordered": bool(series.cat.ordered),
                            "values": encode(np.asarray(series.cat.codes, dtype=np.int64)),
                        })
                    else:
                        encoded_columns.append({
                            "categorical": False,
                            "dtype": str(series.dtype),
                            "values": encode(series.tolist()),
                        })
                return {
                    "t": "dataframe",
                    "columns": columns,
                    "index": encode(value.index.tolist()),
                    "column_values": encoded_columns,
                }
            if isinstance(value, pd.Series):
                return {
                    "t": "series",
                    "name": encode(value.name),
                    "index": encode(value.index.tolist()),
                    "values": encode(value.tolist()),
                    "dtype": str(value.dtype),
                }
            if isinstance(value, pd.Timestamp):
                return {"t": "timestamp", "v": value.isoformat()}
            if isinstance(value, datetime):
                return {"t": "datetime", "v": value.isoformat()}
            if isinstance(value, date):
                return {"t": "date", "v": value.isoformat()}
            if isinstance(value, bytes):
                return {
                    "t": "bytes",
                    "v": base64.b64encode(value).decode("ascii"),
                }
            if hasattr(value, "tocsr") and hasattr(value, "shape"):
                from scipy import sparse

                matrix = value
                if sparse.isspmatrix_csr(matrix):
                    return {
                        "t": "sparse",
                        "format": "csr",
                        "shape": list(matrix.shape),
                        "data": encode(np.asarray(matrix.data)),
                        "indices": encode(np.asarray(matrix.indices)),
                        "indptr": encode(np.asarray(matrix.indptr)),
                    }
                if sparse.isspmatrix_csc(matrix):
                    return {
                        "t": "sparse",
                        "format": "csc",
                        "shape": list(matrix.shape),
                        "data": encode(np.asarray(matrix.data)),
                        "indices": encode(np.asarray(matrix.indices)),
                        "indptr": encode(np.asarray(matrix.indptr)),
                    }
                converted = matrix.tocsr()
                return {
                    "t": "sparse",
                    "format": "csr",
                    "shape": list(converted.shape),
                    "data": encode(np.asarray(converted.data)),
                    "indices": encode(np.asarray(converted.indices)),
                    "indptr": encode(np.asarray(converted.indptr)),
                }
            if isinstance(value, Mapping):
                return {
                    "t": "dict",
                    "items": [[encode(key), encode(item)] for key, item in value.items()],
                }
            if isinstance(value, tuple):
                return {"t": "tuple", "items": [encode(item) for item in value]}
            if isinstance(value, list):
                return {"t": "list", "items": [encode(item) for item in value]}
            raise ArtifactSerializationError(
                f"Unsupported value in structured AnnData fallback: {type(value).__name__}"
            )

        manifest = {
            "version": 2,
            "payload": encode(cls._fallback_payload(payload)),
        }
        manifest_bytes = json.dumps(
            manifest,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=True,
        ).encode("utf-8")
        arrays["__eacbp_manifest__"] = np.frombuffer(manifest_bytes, dtype=np.uint8)
        with target_path.open("wb") as handle:
            np.savez_compressed(handle, **arrays)

    @staticmethod
    def _deserialize_structured_npz(target_path: Path) -> Any:
        try:
            with np.load(target_path, allow_pickle=False) as data:
                if "__eacbp_manifest__" not in data.files:
                    raise ArtifactSerializationError(
                        f"Unsupported structured artifact at {target_path}: "
                        "missing EACBP manifest"
                    )
                manifest_bytes = np.asarray(
                    data["__eacbp_manifest__"], dtype=np.uint8
                ).tobytes()
                manifest = json.loads(manifest_bytes.decode("utf-8"))

                def decode(node: Dict[str, Any]) -> Any:
                    if not isinstance(node, dict) or not isinstance(node.get("t"), str):
                        raise ArtifactSerializationError("Invalid structured payload node")
                    tag = node["t"]
                    if tag == "none":
                        return None
                    if tag == "scalar":
                        return node.get("v")
                    if tag == "complex":
                        return complex(node["real"], node["imag"])
                    if tag == "datetime64":
                        return np.datetime64(node["v"])
                    if tag == "timestamp":
                        return pd.Timestamp(node["v"])
                    if tag == "datetime":
                        return datetime.fromisoformat(node["v"])
                    if tag == "date":
                        return date.fromisoformat(node["v"])
                    if tag == "bytes":
                        return base64.b64decode(node["v"].encode("ascii"))
                    if tag == "ndarray":
                        key = node.get("key")
                        if not isinstance(key, str) or key not in data.files:
                            raise ArtifactSerializationError("Missing structured array payload")
                        value = np.array(data[key], copy=True)
                        expected_shape = tuple(node.get("shape", []))
                        if value.shape != expected_shape:
                            value = value.reshape(expected_shape)
                        return value
                    if tag == "ndarray_object":
                        values = [decode(item) for item in node.get("items", [])]
                        return np.asarray(values, dtype=object).reshape(tuple(node["shape"]))
                    if tag == "categorical":
                        categories = decode(node["categories"])
                        codes = np.asarray(decode(node["codes"]), dtype=np.int64)
                        return pd.Categorical.from_codes(
                            codes, categories=categories, ordered=bool(node.get("ordered", False))
                        )
                    if tag == "dataframe":
                        columns = [decode(item) for item in node.get("columns", [])]
                        index = decode(node["index"])
                        values = node.get("column_values", [])
                        frame = pd.DataFrame(index=index)
                        for position, column in enumerate(columns):
                            column_info = values[position]
                            if column_info.get("categorical"):
                                frame[column] = decode({
                                    "t": "categorical",
                                    "categories": column_info["categories"],
                                    "ordered": column_info["ordered"],
                                    "codes": column_info["values"],
                                })
                            else:
                                series_values = decode(column_info["values"])
                                frame[column] = series_values
                                dtype = column_info.get("dtype")
                                if dtype and dtype != "object":
                                    try:
                                        frame[column] = frame[column].astype(dtype)
                                    except (TypeError, ValueError):
                                        pass
                        frame.index = index
                        frame.columns = columns
                        return frame
                    if tag == "series":
                        series = pd.Series(
                            decode(node["values"]),
                            index=decode(node["index"]),
                            name=decode(node["name"]),
                        )
                        dtype = node.get("dtype")
                        if dtype and dtype != "object":
                            try:
                                series = series.astype(dtype)
                            except (TypeError, ValueError):
                                pass
                        return series
                    if tag == "sparse":
                        from scipy import sparse

                        shape = tuple(node["shape"])
                        data_array = decode(node["data"])
                        indices = decode(node["indices"])
                        indptr = decode(node["indptr"])
                        if node["format"] == "csc":
                            return sparse.csc_matrix((data_array, indices, indptr), shape=shape)
                        return sparse.csr_matrix((data_array, indices, indptr), shape=shape)
                    if tag == "dict":
                        return {decode(key): decode(value) for key, value in node.get("items", [])}
                    if tag == "tuple":
                        return tuple(decode(item) for item in node.get("items", []))
                    if tag == "list":
                        return [decode(item) for item in node.get("items", [])]
                    raise ArtifactSerializationError(f"Unknown structured payload tag: {tag!r}")

                if manifest.get("version") != 2:
                    raise ArtifactSerializationError("Unsupported structured payload version")
                return decode(manifest["payload"])
        except ArtifactSerializationError:
            raise
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ArtifactSerializationError(
                f"Unable to read structured artifact at {target_path}: {exc}"
            ) from exc

    @classmethod
    def serialize(
        cls,
        payload: Any,
        artifact_type: ArtifactType,
        target_path: Path,
        storage_format: Optional[str] = None,
    ) -> Tuple[str, int]:
        """Serialize *payload* to an already-created temporary path.

        The caller publishes the temporary path atomically after this method
        returns.  This method never catches a generic exception, allowing
        callers to distinguish unsupported data from disk failures.
        """

        target_path.parent.mkdir(parents=True, exist_ok=True)

        if artifact_type == ArtifactType.TABLE:
            if isinstance(payload, pd.DataFrame):
                payload.to_csv(target_path, index=True)
            elif isinstance(payload, (list, dict)):
                pd.DataFrame(payload).to_csv(target_path, index=True)
            else:
                raise ArtifactSerializationError(
                    f"Cannot serialize Table payload of type {type(payload).__name__}"
                )
        elif artifact_type in (ArtifactType.JSON, ArtifactType.GENE_LIST, ArtifactType.FASTQ, ArtifactType.GRAPH):
            with target_path.open("w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2, ensure_ascii=False, default=_json_default)
        elif artifact_type == ArtifactType.REPORT:
            if not isinstance(payload, str):
                raise ArtifactSerializationError("Report payload must be a string")
            target_path.write_text(payload, encoding="utf-8")
        elif artifact_type == ArtifactType.FIGURE:
            if isinstance(payload, bytes):
                target_path.write_bytes(payload)
            elif hasattr(payload, "savefig"):
                payload.savefig(target_path, bbox_inches="tight", dpi=150)
            else:
                raise ArtifactSerializationError(
                    f"Cannot serialize Figure payload of type {type(payload).__name__}"
                )
        elif artifact_type in cls._ANNDATA_TYPES:
            storage_format = storage_format or (
                "h5ad" if target_path.suffix.lower() == ".h5ad" else "npz"
            )
            if storage_format == "h5ad":
                adata = cls._as_anndata(payload)
                adata.write_h5ad(target_path)
            elif storage_format == "npz":
                cls._serialize_structured_npz(payload, target_path)
            else:
                raise ArtifactSerializationError(
                    f"Unknown AnnData storage format: {storage_format!r}"
                )
        else:
            raise ArtifactSerializationError(
                f"Unsupported artifact type: {artifact_type!r}"
            )

        return _sha256_file(target_path)

    @classmethod
    def deserialize(cls, target_path: Path, artifact_type: ArtifactType) -> Any:
        if not target_path.exists():
            raise FileNotFoundError(f"Artifact payload not found at {target_path}")

        if artifact_type == ArtifactType.TABLE:
            return pd.read_csv(target_path, index_col=0)
        if artifact_type in (ArtifactType.JSON, ArtifactType.GENE_LIST, ArtifactType.FASTQ, ArtifactType.GRAPH):
            with target_path.open("r", encoding="utf-8") as handle:
                return json.load(handle)
        if artifact_type == ArtifactType.REPORT:
            return target_path.read_text(encoding="utf-8")
        if artifact_type == ArtifactType.FIGURE:
            return target_path.read_bytes()
        if artifact_type in cls._ANNDATA_TYPES:
            suffix = target_path.suffix.lower()
            if suffix == ".h5ad":
                ad = cls._get_anndata_module()
                from eacbp.capabilities.sc_data import SCData

                # Read failures are intentionally propagated.  A corrupt h5ad
                # must not silently become an unrelated NPZ payload.
                return SCData.from_anndata(ad.read_h5ad(target_path))
            if suffix == ".npz":
                return cls._deserialize_structured_npz(target_path)
            raise ArtifactSerializationError(
                f"Unsupported AnnData artifact suffix {target_path.suffix!r}"
            )
        raise ArtifactSerializationError(f"Unsupported artifact type: {artifact_type!r}")


class ArtifactStorageBackend:
    """Manage immutable artifact files under one validated base directory."""

    _thread_locks: Dict[str, threading.RLock] = {}
    _thread_locks_guard = threading.Lock()

    def __init__(self, base_dir: str = ".artifacts"):
        self.base_dir = Path(base_dir).expanduser().resolve()
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self._lock_path = self.base_dir / ".registry.lock"

    def _validated_path(self, path: Path) -> Path:
        """Resolve *path* and reject paths outside ``base_dir``."""

        resolved_base = self.base_dir.resolve()
        try:
            resolved = path.expanduser().resolve(strict=False)
            resolved.relative_to(resolved_base)
        except (OSError, ValueError) as exc:
            raise ValueError(
                f"Artifact storage path escapes configured base directory: {path}"
            ) from exc
        return resolved

    def _get_path_for_uri(
        self,
        uri: ArtifactURI,
        artifact_type: ArtifactType,
        storage_format: Optional[str] = None,
    ) -> Path:
        ext_map = {
            ArtifactType.FASTQ: ".json",
            ArtifactType.ANNDATA: ".h5ad",
            ArtifactType.SPATIAL_DATA: ".h5ad",
            ArtifactType.TABLE: ".csv",
            ArtifactType.FIGURE: ".png",
            ArtifactType.GRAPH: ".json",
            ArtifactType.REPORT: ".md",
            ArtifactType.JSON: ".json",
            ArtifactType.GENE_LIST: ".json",
        }
        if artifact_type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA) and storage_format == "npz":
            ext = ".npz"
        else:
            ext = ext_map.get(artifact_type, ".dat")
        candidate = self.base_dir / uri.study_id / uri.name / f"{uri.version}{ext}"
        return self._validated_path(candidate)

    @contextmanager
    def lock(self) -> Iterator[None]:
        """Acquire a process/thread and cross-process lock for registry writes."""

        key = str(self.base_dir).casefold()
        with self._thread_locks_guard:
            thread_lock = self._thread_locks.setdefault(key, threading.RLock())
        with thread_lock:
            self._lock_path.parent.mkdir(parents=True, exist_ok=True)
            with self._lock_path.open("a+b") as handle:
                handle.seek(0, os.SEEK_END)
                if handle.tell() == 0:
                    handle.write(b"0")
                    handle.flush()
                handle.seek(0)
                if os.name == "nt":
                    import msvcrt

                    msvcrt.locking(handle.fileno(), msvcrt.LK_LOCK, 1)
                else:  # pragma: no cover - exercised on POSIX CI
                    import fcntl

                    fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
                try:
                    yield
                finally:
                    if os.name == "nt":
                        import msvcrt

                        handle.seek(0)
                        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
                    else:  # pragma: no cover - exercised on POSIX CI
                        import fcntl

                        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    @staticmethod
    def _lexists(path: Path) -> bool:
        return os.path.lexists(str(path))

    def _find_existing_path(self, uri: ArtifactURI, artifact_type: ArtifactType) -> Path:
        primary = self._get_path_for_uri(uri, artifact_type)
        if self._lexists(primary):
            return self._validated_path(primary)
        if artifact_type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
            explicit_fallback = self._get_path_for_uri(uri, artifact_type, "npz")
            if self._lexists(explicit_fallback):
                return self._validated_path(explicit_fallback)
        return primary

    def save(
        self,
        uri_str: str,
        payload: Any,
        artifact_type: ArtifactType,
    ) -> Tuple[Path, str, int]:
        """Serialize and publish one immutable artifact atomically."""

        uri = ArtifactURI.parse(uri_str)
        storage_format = PayloadSerializer.preferred_format(payload, artifact_type)
        target_path = self._get_path_for_uri(uri, artifact_type, storage_format)
        target_path.parent.mkdir(parents=True, exist_ok=True)

        if self._lexists(target_path):
            raise ArtifactAlreadyExistsError(
                f"Artifact '{uri_str}' already exists at '{target_path}'. "
                "Immutable artifact versions cannot be overwritten."
            )

        temp_path: Optional[Path] = None
        try:
            fd, temp_name = tempfile.mkstemp(
                prefix=f".{target_path.name}.",
                suffix=".tmp",
                dir=str(target_path.parent),
            )
            os.close(fd)
            temp_path = Path(temp_name)
            sha256_hash, size_bytes = PayloadSerializer.serialize(
                payload,
                artifact_type,
                temp_path,
                storage_format=storage_format,
            )

            # A hard-link is an atomic create-if-absent operation on the same
            # filesystem.  Unlike os.replace(), it cannot clobber a competing
            # writer that won the race after our initial existence check.
            if self._lexists(target_path):
                raise ArtifactAlreadyExistsError(
                    f"Artifact '{uri_str}' was created by another writer"
                )
            try:
                os.link(temp_path, target_path)
            except FileExistsError as exc:
                raise ArtifactAlreadyExistsError(
                    f"Artifact '{uri_str}' was created by another writer"
                ) from exc
            except OSError as exc:
                raise ArtifactStorageError(
                    f"Atomic publication of artifact '{uri_str}' failed: {exc}"
                ) from exc
            return target_path, sha256_hash, size_bytes
        finally:
            if temp_path is not None:
                try:
                    temp_path.unlink(missing_ok=True)
                except OSError:
                    pass

    def remove(self, storage_path: str | Path) -> None:
        """Remove a just-published path after a failed metadata commit."""

        path = self._validated_path(Path(storage_path))
        if path == self.base_dir or path == self._lock_path:
            raise ValueError("Refusing to remove the storage root or registry lock")
        path.unlink(missing_ok=True)

    def load(
        self,
        uri_str: str,
        artifact_type: ArtifactType,
        expected_sha256: Optional[str] = None,
        storage_path: Optional[str] = None,
    ) -> Any:
        uri = ArtifactURI.parse(uri_str)
        if storage_path is None:
            target_path = self._find_existing_path(uri, artifact_type)
        else:
            raw_path = Path(storage_path)
            if not raw_path.is_absolute():
                raw_path = self.base_dir / raw_path
            target_path = self._validated_path(raw_path)
        if not target_path.exists():
            raise FileNotFoundError(f"Artifact payload not found at {target_path}")

        if expected_sha256:
            actual_hash, _ = _sha256_file(target_path)
            if actual_hash != expected_sha256:
                raise ArtifactIntegrityError(
                    f"SHA-256 mismatch for '{uri_str}': expected {expected_sha256}, "
                    f"observed {actual_hash}"
                )
        return PayloadSerializer.deserialize(target_path, artifact_type)

    def exists(self, uri_str: str, artifact_type: ArtifactType) -> bool:
        uri = ArtifactURI.parse(uri_str)
        path = self._find_existing_path(uri, artifact_type)
        return self._lexists(path)
