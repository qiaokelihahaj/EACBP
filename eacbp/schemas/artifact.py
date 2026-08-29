"""
Artifact schema defining metadata, immutability guarantees, hashes, and lineage tracking.
"""

from typing import List, Dict, Any, Optional
from enum import Enum
from datetime import datetime, timezone
from pydantic import BaseModel, Field


class ArtifactType(str, Enum):
    FASTQ = "FASTQ"
    ANNDATA = "AnnData"
    SPATIAL_DATA = "SpatialData"
    TABLE = "Table"
    FIGURE = "Figure"
    GRAPH = "Graph"
    REPORT = "Report"
    JSON = "JSON"
    GENE_LIST = "GeneList"


class ArtifactRef(BaseModel):
    uri: str = Field(..., description="Uniform Resource Identifier, e.g., adata://AD_mouse_001/microglia_subset/v4")
    artifact_id: str = Field(..., description="Unique identifier for the artifact")
    type: ArtifactType = Field(ArtifactType.ANNDATA)


class ArtifactMetadata(BaseModel):
    artifact_id: str = Field(..., description="Unique artifact ID or canonical URI")
    uri: str = Field(..., description="Canonical URI, e.g. adata://study/name/vN")
    type: ArtifactType = Field(...)
    study_id: str = Field(...)
    
    # Lineage / Provenance
    parent_uris: List[str] = Field(default_factory=list, description="Direct parent artifact URIs")
    created_by_task: str = Field(..., description="Task ID that produced this artifact")
    operation: str = Field(..., description="High-level operation performed, e.g., subset_cells, normalize")
    parameters: Dict[str, Any] = Field(default_factory=dict, description="Parameters used to generate artifact")
    
    # Environment & Reproducibility
    software_versions: Dict[str, str] = Field(default_factory=dict, description="Versions of python/scanpy/packages")
    container_image: Optional[str] = Field(None, description="Container tag if executed in container")
    random_seed: int = Field(42, description="Random seed configured")
    sha256_hash: str = Field(..., description="SHA-256 content checksum")
    
    # File Storage
    storage_path: str = Field(..., description="Local or remote absolute path to the persisted artifact payload")
    size_bytes: int = Field(0, description="Size of artifact file in bytes")
    created_at: datetime = Field(default_factory=lambda: datetime.now(timezone.utc))
    
    # Domain-specific quick summaries (e.g. n_obs, n_vars, clusters, deg_count)
    summary_metrics: Dict[str, Any] = Field(default_factory=dict)


class LineageNode(BaseModel):
    uri: str
    artifact_id: str
    type: ArtifactType
    parents: List[str] = Field(default_factory=list)
    children: List[str] = Field(default_factory=list)
    task_id: str
    operation: str
    created_at: datetime
