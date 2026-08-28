"""
Artifact Store with strict content-addressed hashing and immutability guarantees.
"""

import os
import io
import json
import hashlib
from typing import Any, Tuple, Optional, Dict
from pathlib import Path
import pandas as pd
import numpy as np

from eacbp.schemas.artifact import ArtifactType, ArtifactMetadata
from eacbp.artifact.uri import ArtifactURI


class ArtifactAlreadyExistsError(Exception):
    """Raised when attempting to overwrite an immutable artifact version."""
    pass


class PayloadSerializer:
    """Handles serialization and deserialization for various artifact payload types."""

    @staticmethod
    def serialize(payload: Any, artifact_type: ArtifactType, target_path: Path) -> Tuple[str, int]:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        hasher = hashlib.sha256()

        if artifact_type == ArtifactType.TABLE:
            if isinstance(payload, pd.DataFrame):
                payload.to_csv(target_path, index=True)
            elif isinstance(payload, list) or isinstance(payload, dict):
                pd.DataFrame(payload).to_csv(target_path, index=True)
            else:
                raise TypeError(f"Cannot serialize Table payload of type {type(payload)}")

        elif artifact_type in (ArtifactType.JSON, ArtifactType.GENE_LIST):
            with open(target_path, "w", encoding="utf-8") as f:
                json.dump(payload, f, indent=2, default=str)

        elif artifact_type == ArtifactType.REPORT:
            with open(target_path, "w", encoding="utf-8") as f:
                f.write(str(payload))

        elif artifact_type == ArtifactType.FIGURE:
            if isinstance(payload, bytes):
                with open(target_path, "wb") as f:
                    f.write(payload)
            elif hasattr(payload, "savefig"):  # Matplotlib figure
                payload.savefig(target_path, bbox_inches="tight", dpi=150)
            else:
                with open(target_path, "wb") as f:
                    f.write(bytes(payload))

        elif artifact_type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
            # Check if AnnData is available and object is AnnData
            try:
                import anndata as ad
                if isinstance(payload, ad.AnnData):
                    payload.write_h5ad(target_path)
                    with open(target_path, "rb") as f:
                        data = f.read()
                        hasher.update(data)
                    return f"sha256:{hasher.hexdigest()}", len(data)
            except ImportError:
                pass

            # Fallback lightweight AnnData dictionary / custom serializer
            if isinstance(payload, dict):
                npz_dict = {}
                for k, v in payload.items():
                    if k == "obsm" and isinstance(v, dict):
                        for obsm_k, obsm_v in v.items():
                            npz_dict[f"__obsm__{obsm_k}"] = np.asarray(obsm_v, dtype=np.float32)
                    elif isinstance(v, (np.ndarray, list)):
                        npz_dict[k] = np.asarray(v)
                    elif isinstance(v, pd.DataFrame):
                        # Convert dataframe to records
                        npz_dict[f"__df__{k}"] = np.array(json.dumps(v.to_dict(orient="records"), default=str))
                    elif isinstance(v, dict):
                        npz_dict[f"__dict__{k}"] = np.array(json.dumps(v, default=str))
                    else:
                        npz_dict[f"__val__{k}"] = np.array(str(v))
                
                with open(target_path, "wb") as f:
                    np.savez_compressed(f, **npz_dict)
            else:
                with open(target_path, "w", encoding="utf-8") as f:
                    json.dump(payload, f, default=str)

        else:
            with open(target_path, "w", encoding="utf-8") as f:
                f.write(str(payload))

        # Calculate hash and size
        with open(target_path, "rb") as f:
            data = f.read()
            hasher.update(data)
            size = len(data)

        return f"sha256:{hasher.hexdigest()}", size

    @staticmethod
    def deserialize(target_path: Path, artifact_type: ArtifactType) -> Any:
        if not target_path.exists():
            raise FileNotFoundError(f"Artifact payload not found at {target_path}")

        if artifact_type == ArtifactType.TABLE:
            return pd.read_csv(target_path, index_col=0)

        elif artifact_type in (ArtifactType.JSON, ArtifactType.GENE_LIST):
            with open(target_path, "r", encoding="utf-8") as f:
                return json.load(f)

        elif artifact_type == ArtifactType.REPORT:
            with open(target_path, "r", encoding="utf-8") as f:
                return f.read()

        elif artifact_type == ArtifactType.FIGURE:
            with open(target_path, "rb") as f:
                return f.read()

        elif artifact_type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA):
            try:
                import anndata as ad
                if target_path.suffix == ".h5ad":
                    return ad.read_h5ad(target_path)
            except (ImportError, Exception):
                pass

            try:
                data = np.load(target_path, allow_pickle=True)
                res = {"obsm": {}, "uns": {}}
                for k in data.files:
                    if k.startswith("__obsm__"):
                        real_k = k[8:]
                        res["obsm"][real_k] = data[k]
                    elif k.startswith("__df__"):
                        real_k = k[6:]
                        raw_str = str(data[k])
                        records = json.loads(raw_str)
                        res[real_k] = pd.DataFrame(records)
                    elif k.startswith("__dict__"):
                        real_k = k[8:]
                        raw_str = str(data[k])
                        res[real_k] = json.loads(raw_str)
                    elif k.startswith("__val__"):
                        real_k = k[7:]
                        res[real_k] = str(data[k])
                    else:
                        res[k] = data[k]
                return res
            except Exception:
                pass

        with open(target_path, "r", encoding="utf-8", errors="ignore") as f:
            return f.read()


class ArtifactStorageBackend:
    """Manages physical files and strict versioning on local disk."""

    def __init__(self, base_dir: str = ".artifacts"):
        self.base_dir = Path(base_dir).resolve()
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def _get_path_for_uri(self, uri: ArtifactURI, artifact_type: ArtifactType) -> Path:
        ext_map = {
            ArtifactType.ANNDATA: ".h5ad",
            ArtifactType.SPATIAL_DATA: ".h5ad",
            ArtifactType.TABLE: ".csv",
            ArtifactType.FIGURE: ".png",
            ArtifactType.GRAPH: ".json",
            ArtifactType.REPORT: ".md",
            ArtifactType.JSON: ".json",
            ArtifactType.GENE_LIST: ".json",
        }
        ext = ext_map.get(artifact_type, ".dat")
        return self.base_dir / uri.study_id / uri.name / f"{uri.version}{ext}"

    def save(
        self,
        uri_str: str,
        payload: Any,
        artifact_type: ArtifactType
    ) -> Tuple[Path, str, int]:
        uri = ArtifactURI.parse(uri_str)
        target_path = self._get_path_for_uri(uri, artifact_type)

        if target_path.exists():
            raise ArtifactAlreadyExistsError(
                f"Artifact '{uri_str}' already exists at '{target_path}'. "
                f"EACBP Invariant 2 prohibits in-place overwrites. Please use a new version (e.g. {uri.next_version()}) or branch."
            )

        sha256_hash, size_bytes = PayloadSerializer.serialize(payload, artifact_type, target_path)
        return target_path, sha256_hash, size_bytes

    def load(self, uri_str: str, artifact_type: ArtifactType) -> Any:
        uri = ArtifactURI.parse(uri_str)
        target_path = self._get_path_for_uri(uri, artifact_type)
        if not target_path.exists():
            # Check fallback npz/json if h5ad was saved as npz
            if artifact_type in (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA) and target_path.with_suffix(".npz").exists():
                target_path = target_path.with_suffix(".npz")
        return PayloadSerializer.deserialize(target_path, artifact_type)

    def exists(self, uri_str: str, artifact_type: ArtifactType) -> bool:
        uri = ArtifactURI.parse(uri_str)
        target_path = self._get_path_for_uri(uri, artifact_type)
        return target_path.exists() or target_path.with_suffix(".npz").exists()
