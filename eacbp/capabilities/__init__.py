"""
EACBP Compute Plane Capabilities.
"""

from eacbp.capabilities.base import BaseCapability, ImplementationType
from eacbp.capabilities.side_effect import SideEffectValidator
from eacbp.capabilities.registry import CapabilityRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.qc import DatasetAuditCapability, QCCapability
from eacbp.capabilities.normalization import NormalizationCapability
from eacbp.capabilities.integration import IntegrationCapability
from eacbp.capabilities.clustering import ClusteringCapability
from eacbp.capabilities.subset import SubsetCapability
from eacbp.capabilities.deg import DifferentialAbundanceCapability, DifferentialExpressionCapability
from eacbp.capabilities.trajectory import TrajectoryCapability


from eacbp.capabilities.quantification import FASTQQuantificationCapability


def create_default_capability_registry() -> CapabilityRegistry:
    """Instantiates and registers all standard built-in single-cell capabilities."""
    reg = CapabilityRegistry()
    reg.register(FASTQQuantificationCapability(implementation_id="kb_python_v1"))
    reg.register(FASTQQuantificationCapability(implementation_id="starsolo_v1"))
    reg.register(FASTQQuantificationCapability(implementation_id="sc_quant_v1"))
    reg.register(DatasetAuditCapability())
    reg.register(QCCapability())
    reg.register(NormalizationCapability())
    reg.register(IntegrationCapability(implementation_id="harmony"))
    reg.register(IntegrationCapability(implementation_id="no_correction"))
    reg.register(ClusteringCapability())
    reg.register(SubsetCapability())
    reg.register(DifferentialAbundanceCapability())
    reg.register(DifferentialExpressionCapability())
    reg.register(TrajectoryCapability(implementation_id="paga_dpt"))
    reg.register(TrajectoryCapability(implementation_id="cellrank"))
    return reg


__all__ = [
    "BaseCapability",
    "ImplementationType",
    "SideEffectValidator",
    "CapabilityRegistry",
    "SCData",
    "DatasetAuditCapability",
    "QCCapability",
    "NormalizationCapability",
    "IntegrationCapability",
    "ClusteringCapability",
    "SubsetCapability",
    "DifferentialAbundanceCapability",
    "DifferentialExpressionCapability",
    "TrajectoryCapability",
    "FASTQQuantificationCapability",
    "create_default_capability_registry",
]
