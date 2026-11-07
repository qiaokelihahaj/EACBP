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
from pathlib import Path, PurePosixPath
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
    TABLE_FORMAT = "table_npz_v1"
    TABLE_CONTRACT_KEY = "__eacbp_table_contract__"

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

        if artifact_type == ArtifactType.TABLE:
            return cls.TABLE_FORMAT
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

    @staticmethod
    def _table_require(condition: bool, message: str) -> None:
        if not condition:
            raise ArtifactSerializationError("Invalid Table v1 contract: " + message)

    @classmethod
    def _serialize_table_npz(cls, payload: Any, target_path: Path) -> None:
        """Persist table values and identity without CSV type/NA inference.

        This is a separate versioned contract from the legacy AnnData fallback.
        Native numeric arrays retain their binary dtype. Nullable columns carry
        an explicit mask; object cells use primitive tags that distinguish an
        empty string, NA-looking text, None, NaN and pandas.NA. Unsupported
        types fail before publication rather than being stringified.
        """
        if isinstance(payload, pd.DataFrame):
            frame = payload
        elif isinstance(payload, (list, dict)):
            frame = pd.DataFrame(payload)
        else:
            raise ArtifactSerializationError(f"Cannot serialize Table payload of type {type(payload).__name__}")
        if frame.attrs:
            raise ArtifactSerializationError("Table v1 does not support DataFrame.attrs; use artifact metadata")
        arrays: Dict[str, np.ndarray] = {}

        def array(value: np.ndarray) -> Dict[str, Any]:
            value = np.asarray(value)
            cls._table_require(not value.dtype.hasobject, "binary arrays cannot contain objects")
            key = f"array_{len(arrays)}"
            arrays[key] = np.array(value, copy=True)
            return {"key": key, "dtype": value.dtype.str, "shape": list(value.shape)}

        def scalar(value: Any) -> Dict[str, Any]:
            if value is None:
                return {"t": "none"}
            if value is pd.NA:
                return {"t": "pd_na"}
            if isinstance(value, np.generic):
                cls._table_require(value.dtype.kind in "biufU", "unsupported numpy object-cell scalar")
                return {"t": "numpy_scalar", "array": array(np.asarray(value))}
            if type(value) is bool:
                return {"t": "bool", "v": value}
            if type(value) is int:
                return {"t": "int", "v": str(value)}
            if type(value) is float:
                return {"t": "float", "v": value.hex()}
            if type(value) is str:
                return {"t": "str", "v": value}
            raise ArtifactSerializationError(f"Unsupported Table v1 scalar: {type(value).__name__}")

        def vector(values: Any) -> Dict[str, Any]:
            dtype = values.dtype
            if isinstance(dtype, pd.CategoricalDtype):
                categorical = values.array if isinstance(values, (pd.Series, pd.Index)) else values
                return {"kind": "categorical", "dtype": "category", "ordered": bool(dtype.ordered),
                        "categories": axis(categorical.categories), "codes": array(np.asarray(categorical.codes))}
            if isinstance(dtype, pd.StringDtype):
                cls._table_require(dtype.storage == "python", "only python-backed StringDtype is supported")
                missing = "pd.NA" if dtype.na_value is pd.NA else "NaN"
                mask = np.asarray(pd.isna(values), dtype=bool)
                text = np.asarray(["" if flag else value for value, flag in zip(values, mask)], dtype=str)
                return {"kind": "string", "dtype": "string", "storage": "python", "missing": missing,
                        "values": array(text), "mask": array(mask)}
            nullable = {"Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64", "Float32", "Float64", "boolean"}
            if isinstance(dtype, pd.api.extensions.ExtensionDtype):
                name = str(dtype)
                cls._table_require(name in nullable, "unsupported extension dtype: " + name)
                native = np.dtype("bool" if name == "boolean" else name.lower())
                mask = np.asarray(pd.isna(values), dtype=bool)
                return {"kind": "nullable", "dtype": name, "missing": "pd.NA", "mask": array(mask),
                        "values": array(values.to_numpy(dtype=native, na_value=0))}
            dtype = np.dtype(dtype)
            if dtype == np.dtype("object"):
                return {"kind": "object", "dtype": "object", "values": [scalar(value) for value in values]}
            cls._table_require(dtype.kind in "biuf" and dtype.itemsize in (1, 2, 4, 8), "unsupported native dtype: " + str(dtype))
            return {"kind": "native", "dtype": dtype.str, "values": array(np.asarray(values))}

        def axis(value: pd.Index) -> Dict[str, Any]:
            if isinstance(value, pd.RangeIndex):
                return {"kind": "range", "start": value.start, "stop": value.stop, "step": value.step, "name": scalar(value.name)}
            cls._table_require(not isinstance(value, pd.MultiIndex), "MultiIndex is not supported")
            return {"kind": "categorical_index" if isinstance(value, pd.CategoricalIndex) else "index",
                    "name": scalar(value.name), "values": vector(value)}

        manifest = {"format": "eacbp.table.npz", "version": 1, "decoder": "eacbp.table.v1", "shape": list(frame.shape),
                    "index": axis(frame.index), "columns": axis(frame.columns),
                    "column_values": [vector(frame.iloc[:, position]) for position in range(frame.shape[1])],
                    "allows_duplicate_labels": bool(frame.flags.allows_duplicate_labels)}
        encoded = json.dumps(manifest, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")
        arrays[cls.TABLE_CONTRACT_KEY] = np.frombuffer(encoded, dtype=np.uint8)
        with target_path.open("wb") as handle:
            np.savez_compressed(handle, **arrays)

    @classmethod
    def _deserialize_table_npz(cls, target_path: Path) -> pd.DataFrame:
        """Decode Table v1 strictly; declared type/shape never triggers coercion."""
        try:
            with np.load(target_path, allow_pickle=False) as stored:
                cls._table_require(cls.TABLE_CONTRACT_KEY in stored.files, "missing table contract")
                header = stored[cls.TABLE_CONTRACT_KEY]
                cls._table_require(header.dtype == np.dtype("uint8") and header.ndim == 1, "contract must be UTF-8 bytes")
                manifest = json.loads(header.tobytes().decode("utf-8"))
                cls._table_require(isinstance(manifest, dict) and manifest.get("format") == "eacbp.table.npz"
                                   and type(manifest.get("version")) is int and manifest["version"] == 1
                                   and manifest.get("decoder") == "eacbp.table.v1", "unsupported format/version/decoder")
                shape = manifest.get("shape")
                cls._table_require(isinstance(shape, list) and len(shape) == 2 and all(type(n) is int and n >= 0 for n in shape), "invalid table shape")
                used = {cls.TABLE_CONTRACT_KEY}

                def array(info: Dict[str, Any], expected_shape: Tuple[int, ...]) -> np.ndarray:
                    cls._table_require(isinstance(info, dict) and set(info) == {"key", "dtype", "shape"}, "invalid array declaration")
                    key = info["key"]
                    cls._table_require(isinstance(key, str) and key in stored.files and key not in used, "missing/reused array key")
                    value = np.array(stored[key], copy=True)
                    cls._table_require(not value.dtype.hasobject and value.dtype.str == info["dtype"]
                                       and list(value.shape) == info["shape"] and value.shape == expected_shape, "array dtype/shape does not match contract")
                    used.add(key)
                    return value

                def scalar(node: Dict[str, Any]) -> Any:
                    cls._table_require(isinstance(node, dict), "invalid scalar tag")
                    tag = node.get("t")
                    if tag in ("none", "pd_na"):
                        cls._table_require(set(node) == {"t"}, "invalid missing scalar")
                        return None if tag == "none" else pd.NA
                    if tag == "numpy_scalar":
                        cls._table_require(set(node) == {"t", "array"}, "invalid numpy scalar")
                        value = array(node["array"], ())
                        cls._table_require(value.dtype.kind in "biufU", "unsupported numpy scalar dtype")
                        return value[()]
                    cls._table_require(set(node) == {"t", "v"}, "invalid primitive scalar")
                    value = node["v"]
                    if tag == "str":
                        cls._table_require(type(value) is str, "text scalar must remain text")
                        return value
                    if tag == "bool":
                        cls._table_require(type(value) is bool, "boolean scalar must remain bool")
                        return value
                    if tag == "int":
                        cls._table_require(type(value) is str and str(int(value)) == value, "invalid integer scalar")
                        return int(value)
                    if tag == "float":
                        cls._table_require(type(value) is str and float.fromhex(value).hex() == value, "invalid floating scalar")
                        return float.fromhex(value)
                    raise ArtifactSerializationError("Invalid Table v1 contract: unknown scalar tag")

                def vector(node: Dict[str, Any], length: int) -> Any:
                    cls._table_require(isinstance(node, dict), "invalid column declaration")
                    kind = node.get("kind")
                    if kind == "object":
                        cls._table_require(set(node) == {"kind", "dtype", "values"} and node["dtype"] == "object"
                                           and isinstance(node["values"], list) and len(node["values"]) == length, "invalid object dtype/length")
                        value = np.empty(length, dtype=object)
                        value[:] = [scalar(item) for item in node["values"]]
                        return value
                    if kind == "native":
                        cls._table_require(set(node) == {"kind", "dtype", "values"}, "invalid native declaration")
                        value = array(node["values"], (length,))
                        cls._table_require(value.dtype.str == node["dtype"] and value.dtype.kind in "biuf"
                                           and value.dtype.itemsize in (1, 2, 4, 8), "native dtype mismatch")
                        return value
                    if kind in ("nullable", "string"):
                        expected_keys = {"kind", "dtype", "missing", "values", "mask"} | ({"storage"} if kind == "string" else set())
                        cls._table_require(set(node) == expected_keys, "invalid nullable declaration")
                        mask = array(node["mask"], (length,))
                        cls._table_require(mask.dtype == np.dtype("bool"), "missing mask must be boolean")
                        value = array(node["values"], (length,))
                        if kind == "string":
                            cls._table_require(node["dtype"] == "string" and node["storage"] == "python" and node["missing"] in ("pd.NA", "NaN")
                                               and value.dtype.kind == "U", "invalid string dtype/decoder")
                            dtype = pd.StringDtype(storage="python", na_value=pd.NA if node["missing"] == "pd.NA" else np.nan)
                            result = pd.array(value, dtype=dtype)
                            result[mask] = dtype.na_value
                            return result
                        nullable = {"Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64", "Float32", "Float64", "boolean"}
                        name = node["dtype"]
                        cls._table_require(name in nullable and node["missing"] == "pd.NA"
                                           and value.dtype == np.dtype("bool" if name == "boolean" else name.lower()), "nullable dtype mismatch")
                        constructor = pd.arrays.BooleanArray if name == "boolean" else pd.arrays.FloatingArray if name.startswith("Float") else pd.arrays.IntegerArray
                        return constructor(value, mask)
                    if kind == "categorical":
                        cls._table_require(set(node) == {"kind", "dtype", "ordered", "categories", "codes"}
                                           and node["dtype"] == "category" and type(node["ordered"]) is bool, "invalid categorical declaration")
                        categories = axis(node["categories"])
                        codes = array(node["codes"], (length,))
                        cls._table_require(codes.dtype.kind == "i" and np.all(codes >= -1) and np.all(codes < len(categories)), "invalid categorical codes")
                        return pd.Categorical.from_codes(codes, categories=categories, ordered=node["ordered"])
                    raise ArtifactSerializationError("Invalid Table v1 contract: unsupported column kind")

                def axis(node: Dict[str, Any], length: Optional[int] = None) -> pd.Index:
                    cls._table_require(isinstance(node, dict), "invalid axis declaration")
                    kind = node.get("kind")
                    if kind == "range":
                        cls._table_require(set(node) == {"kind", "start", "stop", "step", "name"}
                                           and all(type(node[key]) is int for key in ("start", "stop", "step")) and node["step"] != 0, "invalid range axis")
                        result = pd.RangeIndex(node["start"], node["stop"], node["step"], name=scalar(node["name"]))
                    else:
                        cls._table_require(kind in ("index", "categorical_index") and set(node) == {"kind", "name", "values"}, "unsupported axis kind")
                        values_info = node["values"]
                        if length is None:
                            if values_info.get("kind") == "object":
                                length = len(values_info["values"])
                            else:
                                ref = values_info.get("codes") if values_info.get("kind") == "categorical" else values_info.get("values")
                                length = ref["shape"][0]
                        values = vector(values_info, length)
                        result = pd.CategoricalIndex(values, name=scalar(node["name"])) if kind == "categorical_index" else pd.Index(values, dtype=values.dtype, name=scalar(node["name"]))
                    cls._table_require(length is None or len(result) == length, "axis length differs from table shape")
                    return result

                index, columns = axis(manifest["index"], shape[0]), axis(manifest["columns"], shape[1])
                column_values = manifest.get("column_values")
                cls._table_require(isinstance(column_values, list) and len(column_values) == shape[1], "column count differs from table shape")
                decoded = [vector(node, shape[0]) for node in column_values]
                frame = pd.DataFrame({position: pd.Series(value, dtype=value.dtype) for position, value in enumerate(decoded)}, index=pd.RangeIndex(shape[0]))
                frame.index, frame.columns = index, columns
                cls._table_require(type(manifest.get("allows_duplicate_labels")) is bool, "invalid duplicate-label policy")
                frame.flags.allows_duplicate_labels = manifest["allows_duplicate_labels"]
                cls._table_require(set(stored.files) == used, "unbound array payload")
                return frame
        except ArtifactSerializationError:
            raise
        except (OSError, ValueError, TypeError, KeyError, IndexError, OverflowError) as exc:
            raise ArtifactSerializationError(f"Unable to read Table v1 artifact at {target_path}: {exc}") from exc

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
                    "index_name": encode(value.index.name),
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
                        if "index_name" in node:
                            frame.index.name = decode(node["index_name"])
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
            storage_format = storage_format or cls.TABLE_FORMAT
            if storage_format != cls.TABLE_FORMAT:
                raise ArtifactSerializationError(f"Unknown Table storage format: {storage_format!r}")
            cls._serialize_table_npz(payload, target_path)
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
            if target_path.name.lower().endswith(".table.npz"):
                return cls._deserialize_table_npz(target_path)
            if target_path.suffix.lower() != ".csv":
                raise ArtifactSerializationError(f"Unsupported Table artifact suffix {target_path.suffix!r}")
            # Legacy CSV's inference is deliberately unchanged. New tables
            # carry a binary versioned contract; old files are never migrated.
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
                return SCData.from_anndata(
                    ad.read_h5ad(target_path), add_missing_annotations=False
                )
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

    def resolve_path(self, storage_path: str | Path) -> Path:
        """Resolve a stored object key (or a legacy absolute path) safely."""
        path = Path(storage_path)
        return self._validated_path(path if path.is_absolute() else self.base_dir / path)

    @staticmethod
    def validate_object_key(value: str) -> str:
        """Reject platform-dependent or traversing portable object keys."""
        path = PurePosixPath(value)
        if (not value or value == "." or "\\" in value or ":" in value
                or path.is_absolute() or path.as_posix() != value
                or any(part in (".", "..") for part in path.parts)):
            raise ValueError(f"Invalid artifact object key: {value!r}")
        return value

    def relative_path(self, storage_path: str | Path) -> str:
        """Return a portable object key beneath this storage root."""
        return self.resolve_path(storage_path).relative_to(self.base_dir).as_posix()

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
        elif artifact_type == ArtifactType.TABLE and storage_format == PayloadSerializer.TABLE_FORMAT:
            ext = ".table.npz"
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
        if artifact_type == ArtifactType.TABLE:
            typed = self._get_path_for_uri(uri, artifact_type, PayloadSerializer.TABLE_FORMAT)
            if self._lexists(typed):
                return self._validated_path(typed)
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
        if artifact_type == ArtifactType.TABLE and self._lexists(self._get_path_for_uri(uri, artifact_type)):
            raise ArtifactAlreadyExistsError(f"Legacy Table artifact '{uri_str}' already exists; its version cannot be replaced by a new format")

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
        target_path = self._payload_path(uri_str, artifact_type, storage_path)
        self.verify_payload(uri_str, artifact_type, expected_sha256, storage_path)
        return PayloadSerializer.deserialize(target_path, artifact_type)

    def _payload_path(self, uri_str, artifact_type, storage_path=None) -> Path:
        uri = ArtifactURI.parse(uri_str)
        return (self._find_existing_path(uri, artifact_type) if storage_path is None
                else self.resolve_path(storage_path))

    def verify_payload(
        self, uri_str: str, artifact_type: ArtifactType,
        expected_sha256: Optional[str] = None, storage_path: Optional[str] = None,
        *, expected_size: Optional[int] = None,
    ) -> Tuple[str, int]:
        """Stream a payload's digest without constructing its scientific object."""
        target_path = self._payload_path(uri_str, artifact_type, storage_path)
        if not target_path.exists():
            raise FileNotFoundError(f"Artifact payload not found at {target_path}")
        actual_hash, actual_size = _sha256_file(target_path)
        if expected_sha256 and actual_hash != expected_sha256:
            raise ArtifactIntegrityError(
                f"SHA-256 mismatch for '{uri_str}': expected {expected_sha256}, "
                f"observed {actual_hash}"
            )
        if expected_size is not None and actual_size != expected_size:
            raise ArtifactIntegrityError(f"Payload size mismatch for '{uri_str}'")
        return actual_hash, actual_size

    def exists(self, uri_str: str, artifact_type: ArtifactType) -> bool:
        uri = ArtifactURI.parse(uri_str)
        path = self._find_existing_path(uri, artifact_type)
        return self._lexists(path)
