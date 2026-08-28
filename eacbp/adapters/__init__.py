"""
EACBP External Agent Adapter Plane.
Standardized interfaces, contract guardrails, and versioned artifact generation for external biological AI agents.
"""

from eacbp.adapters.base import BaseAgentAdapter
from eacbp.adapters.spacell import SpaCellAgentAdapter
from eacbp.adapters.chatcell import ChatCellAdapter
from eacbp.adapters.geneagent import GeneAgentAdapter
from eacbp.capabilities.registry import CapabilityRegistry


def register_all_adapters(registry: CapabilityRegistry) -> None:
    """Registers all external agent adapters into the given CapabilityRegistry."""
    registry.register(SpaCellAgentAdapter())
    registry.register(ChatCellAdapter())
    registry.register(GeneAgentAdapter())


__all__ = [
    "BaseAgentAdapter",
    "SpaCellAgentAdapter",
    "ChatCellAdapter",
    "GeneAgentAdapter",
    "register_all_adapters",
]
