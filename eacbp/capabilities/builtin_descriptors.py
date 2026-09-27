"""Gradual migration of legacy built-ins to enforced parameter/type contracts."""

from eacbp.capabilities.base import CapabilityDescriptor
from eacbp.capabilities.parameters import (
    NormalizationParameters, PyDESeq2Parameters, QCParameters, WelchDEGParameters,
)
from eacbp.schemas.artifact import ArtifactType


_PARAMETERS = {
    ("qc", "sc_qc_v1"): QCParameters,
    ("normalization", "library_size_log1p_v1"): NormalizationParameters,
    ("deg", "donor_pseudobulk_welch_v1"): WelchDEGParameters,
    ("deg", "pydeseq2_pseudobulk_v1"): PyDESeq2Parameters,
}


def descriptor_for_builtin(capability):
    model = _PARAMETERS.get((capability.capability_name, capability.implementation_id))
    if model is None:
        return None
    descriptor = CapabilityDescriptor.from_capability(capability)
    descriptor.parameter_model = model
    # Both registry types deserialize to SCData. These matrix methods operate
    # on counts/observations and also accept spatial coordinates carried in
    # obsm; their historical constructor metadata listed only AnnData.
    descriptor.input_types = (ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA)
    descriptor.validate_types = True
    descriptor.scope = "per_target" if capability.capability_name == "deg" else "shared"
    return descriptor
