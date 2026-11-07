"""Gradual migration of legacy built-ins to enforced parameter/type contracts."""

from eacbp.capabilities.base import CapabilityDescriptor
from eacbp.capabilities.parameters import (
    NormalizationParameters, PyDESeq2Parameters, QCParameters, WelchDEGParameters,
    LeidenParameters, NoCorrectionParameters, HarmonyParameters,
)
from eacbp.schemas.artifact import ArtifactType


_PARAMETERS = {
    ("qc", "sc_qc_v1"): QCParameters,
    ("normalization", "library_size_log1p_v1"): NormalizationParameters,
    ("deg", "donor_pseudobulk_welch_v1"): WelchDEGParameters,
    ("deg", "pydeseq2_pseudobulk_v1"): PyDESeq2Parameters,
    ("clustering", "scanpy_leiden_umap_v1"): LeidenParameters,
    ("integration", "no_correction_v1"): NoCorrectionParameters,
    ("integration", "harmonypy_v1"): HarmonyParameters,
}


_RUNTIME_PACKAGES = {
    ("deg", "pydeseq2_pseudobulk_v1"): ("pydeseq2",),
    ("integration", "harmonypy_v1"): ("anndata", "scanpy", "harmonypy"),
    ("clustering", "scanpy_leiden_umap_v1"): ("anndata", "scanpy", "igraph", "leidenalg"),
    ("trajectory_inference", "scanpy_dpt_v1"): ("anndata", "scanpy"),
    ("fate_mapping", "cellrank_fate_v1"): ("anndata", "scanpy", "cellrank"),
}


def descriptor_for_builtin(capability):
    key = (capability.capability_name, capability.implementation_id)
    model = _PARAMETERS.get(key)
    packages = _RUNTIME_PACKAGES.get(key, ())
    if model is None and not packages:
        return None
    descriptor = CapabilityDescriptor.from_capability(capability)
    descriptor.runtime_packages = packages
    if model is None:
        return descriptor
    descriptor.parameter_model = model
    if capability.implementation_id == "scanpy_leiden_umap_v1":
        from eacbp.auditor.clustering import ClusteringIntegrityValidator
        descriptor.validator = ClusteringIntegrityValidator()
        descriptor.required_audit_ids = ("leiden_execution_integrity",)
    # Both registry types deserialize to SCData. These matrix methods operate
    # on counts/observations and also accept spatial coordinates carried in
    # obsm; their historical constructor metadata listed only AnnData.
    descriptor.input_types = (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA)
    descriptor.validate_types = True
    descriptor.scope = "per_target" if capability.capability_name == "deg" else "shared"
    return descriptor
