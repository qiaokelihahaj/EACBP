"""
Base Capability definition and implementation interface.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional
from enum import Enum
from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.registry import ArtifactRegistry


class ImplementationType(str, Enum):
    PYTHON_TOOL = "python_tool"
    R_TOOL = "r_tool"
    AGENT_ADAPTER = "agent_adapter"
    CONTAINER = "container"


class BaseCapability(ABC):
    """Abstract interface for all computational capabilities (tools, pipelines, agent adapters)."""

    def __init__(
        self,
        capability_name: str,
        implementation_id: str,
        implementation_type: ImplementationType = ImplementationType.PYTHON_TOOL,
        accepts_modalities: Optional[List[str]] = None,
        accepts_types: Optional[List[ArtifactType]] = None,
        requires_keys: Optional[List[str]] = None,
        suitable_for: Optional[List[str]] = None,
        output_types: Optional[List[ArtifactType]] = None,
    ):
        self.capability_name = capability_name
        self.implementation_id = implementation_id
        self.implementation_type = implementation_type
        self.accepts_modalities = accepts_modalities or ["scRNA"]
        self.accepts_types = accepts_types or [ArtifactType.ANNDATA]
        self.requires_keys = requires_keys or []
        self.suitable_for = suitable_for or []
        self.output_types = output_types or [ArtifactType.ANNDATA]

    @abstractmethod
    def execute(self, contract: TaskContract, registry: ArtifactRegistry) -> TaskResult:
        """Executes the contracted computation, registering new artifacts in registry and returning TaskResult."""
        pass
