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
from eacbp.capabilities.real_methods import HarmonyIntegrationCapability, LeidenClusteringCapability, DPTTrajectoryCapability, CellRankFateCapability
from eacbp.capabilities.spatial import (
    SpatialDomainCapability,
    SpatialDEGCapability,
    CellCellCommunicationCapability,
)
from eacbp.capabilities.perturbation import (
    GeneticPerturbationCapability,
    CompoundPerturbationCapability,
)
from eacbp.adapters import (
    SpaCellAgentAdapter,
    ChatCellAdapter,
    GeneAgentAdapter,
)
from eacbp.knowledge.capability import KnowledgeRetrievalCapability


from eacbp.capabilities.quantification import FASTQQuantificationCapability


def create_default_capability_registry() -> CapabilityRegistry:
    """Build the complete built-in capability registry.

    This is the one default assembly entry point.  The orchestrator receives
    the resulting registry as-is; an explicitly supplied registry is never
    completed or overwritten implicitly.  Keeping optional planes here also
    makes capability discovery consistent for direct Python callers and for
    ``ScientificOrchestrator``.
    """
    reg = CapabilityRegistry()
    reg.register(FASTQQuantificationCapability(implementation_id="kb_python_v1"))
    reg.register(FASTQQuantificationCapability(implementation_id="starsolo_v1"))
    reg.register(FASTQQuantificationCapability(implementation_id="sc_quant_v1"))
    reg.register(DatasetAuditCapability())
    reg.register(QCCapability())
    reg.register(NormalizationCapability())
    reg.register(IntegrationCapability(implementation_id="batch_mean_centering_v1"))
    reg.register(IntegrationCapability(implementation_id="no_correction_v1"))
    reg.register(ClusteringCapability())
    reg.register(SubsetCapability())
    reg.register(DifferentialAbundanceCapability())
    reg.register(DifferentialExpressionCapability())
    reg.register(TrajectoryCapability())
    reg.register(HarmonyIntegrationCapability())
    reg.register(LeidenClusteringCapability())
    reg.register(DPTTrajectoryCapability())
    reg.register(CellRankFateCapability())
    reg.register(SpatialDomainCapability())
    reg.register(SpatialDEGCapability(implementation_id="spatial_deg_morans_i_v1"))
    reg.register(SpatialDEGCapability(implementation_id="spatial_moran_deg_v1"))
    reg.register(CellCellCommunicationCapability())
    reg.register(GeneticPerturbationCapability(implementation_id="in_silico_crispr_ko_v1"))
    reg.register(GeneticPerturbationCapability(implementation_id="in_silico_overexpression_v1"))
    reg.register(CompoundPerturbationCapability())
    reg.register(KnowledgeRetrievalCapability(implementation_id="knowledge_engine_discovery_v1"))
    reg.register(KnowledgeRetrievalCapability(implementation_id="knowledge_engine_prior_v1"))
    reg.register(KnowledgeRetrievalCapability(implementation_id="knowledge_engine_v1"))
    reg.register(SpaCellAgentAdapter())
    reg.register(ChatCellAdapter())
    reg.register(GeneAgentAdapter())
    # Optional imports inside the implementations are lazy. Discovery works
    # without installing heavy dependencies; execution reports their absence.
    from eacbp.capabilities.advanced_statistics import (
        PyDESeq2PseudobulkCapability, DecouplerFunctionalAnalysisCapability,
        PyDESeq2LeaveOneDonorOutCapability,
    )
    from eacbp.capabilities.advanced_qc import (
        ScrubletDoubletCapability, CellTypistAnnotationCapability,
        CellBenderBackgroundRemovalCapability,
    )
    from eacbp.capabilities.advanced_communication import LianaCommunicationCapability
    for capability_class in (
        PyDESeq2PseudobulkCapability, DecouplerFunctionalAnalysisCapability,
        PyDESeq2LeaveOneDonorOutCapability, ScrubletDoubletCapability,
        CellTypistAnnotationCapability, CellBenderBackgroundRemovalCapability,
        LianaCommunicationCapability,
    ):
        reg.register(capability_class())
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
    "HarmonyIntegrationCapability",
    "LeidenClusteringCapability",
    "DPTTrajectoryCapability",
    "CellRankFateCapability",
    "SpatialDomainCapability",
    "SpatialDEGCapability",
    "CellCellCommunicationCapability",
    "GeneticPerturbationCapability",
    "CompoundPerturbationCapability",
    "KnowledgeRetrievalCapability",
    "SpaCellAgentAdapter",
    "ChatCellAdapter",
    "GeneAgentAdapter",
    "create_default_capability_registry",
]
