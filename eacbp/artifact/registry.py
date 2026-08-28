"""
Central Artifact Registry unifying physical storage, SHA-256 metadata verification, and lineage DAG.
"""

import sys
from typing import Any, Tuple, Optional, List, Dict
from eacbp.schemas.artifact import ArtifactMetadata, ArtifactType
from eacbp.artifact.uri import ArtifactURI
from eacbp.artifact.storage import ArtifactStorageBackend
from eacbp.artifact.lineage import LineageGraph


class ArtifactRegistry:
    """Singleton/Instance artifact manager ensuring immutability, lineage, and retrieval."""

    def __init__(self, storage_dir: str = ".artifacts"):
        self.storage = ArtifactStorageBackend(base_dir=storage_dir)
        self.lineage = LineageGraph()
        self.registry: Dict[str, ArtifactMetadata] = {}

    def register(
        self,
        uri_str: str,
        payload: Any,
        artifact_type: ArtifactType,
        study_id: str,
        created_by_task: str,
        operation: str,
        parent_uris: Optional[List[str]] = None,
        parameters: Optional[Dict[str, Any]] = None,
        software_versions: Optional[Dict[str, str]] = None,
        random_seed: int = 42,
        summary_metrics: Optional[Dict[str, Any]] = None,
        container_image: Optional[str] = None,
    ) -> ArtifactMetadata:
        parent_uris = parent_uris or []
        parameters = parameters or {}
        software_versions = software_versions or {
            "python": sys.version.split()[0],
            "eacbp": "0.1.0",
        }
        summary_metrics = summary_metrics or {}

        # Save to disk with immutability enforcement & sha256 calculation
        storage_path, sha256_hash, size_bytes = self.storage.save(
            uri_str=uri_str,
            payload=payload,
            artifact_type=artifact_type
        )

        metadata = ArtifactMetadata(
            artifact_id=uri_str,
            uri=uri_str,
            type=artifact_type,
            study_id=study_id,
            parent_uris=parent_uris,
            created_by_task=created_by_task,
            operation=operation,
            parameters=parameters,
            software_versions=software_versions,
            container_image=container_image,
            random_seed=random_seed,
            sha256_hash=sha256_hash,
            storage_path=str(storage_path),
            size_bytes=size_bytes,
            summary_metrics=summary_metrics,
        )

        self.registry[uri_str] = metadata
        self.lineage.add_artifact(metadata)
        return metadata

    def get_metadata(self, uri_str: str) -> ArtifactMetadata:
        if uri_str not in self.registry:
            raise KeyError(f"Artifact URI '{uri_str}' is not registered in the system.")
        return self.registry[uri_str]

    def load_payload(self, uri_str: str) -> Any:
        meta = self.get_metadata(uri_str)
        return self.storage.load(uri_str, meta.type)

    def get(self, uri_str: str) -> Tuple[ArtifactMetadata, Any]:
        meta = self.get_metadata(uri_str)
        payload = self.storage.load(uri_str, meta.type)
        return meta, payload

    def exists(self, uri_str: str, artifact_type: Optional[ArtifactType] = None) -> bool:
        if uri_str in self.registry:
            return True
        if artifact_type:
            return self.storage.exists(uri_str, artifact_type)
        return False

    def list_artifacts(self, study_id: Optional[str] = None, artifact_type: Optional[ArtifactType] = None) -> List[ArtifactMetadata]:
        results = list(self.registry.values())
        if study_id:
            results = [a for a in results if a.study_id == study_id]
        if artifact_type:
            results = [a for a in results if a.type == artifact_type]
        return results
