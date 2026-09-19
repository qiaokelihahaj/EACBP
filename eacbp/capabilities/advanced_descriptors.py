"""Shared declarations for optional advanced analysis extensions.

The advanced extension planner used to maintain a second ``METHODS`` table
and a second set of audit requirements.  Keeping those values beside the
implementation makes a new extension easy to wire incorrectly.  This module
is intentionally small: it contains discovery metadata and a plan factory,
while execution remains in the real capability classes and auditing remains
independent in ``eacbp.auditor``.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from pydantic import BaseModel, ConfigDict

from eacbp.capabilities.base import CapabilityDescriptor
from eacbp.schemas.artifact import ArtifactType


class AdvancedExtensionParameters(BaseModel):
    """Permissive model for extension parameters.

    Advanced methods have method-specific validation at execution time.  The
    declaration still validates that the contract is a mapping, and
    ``extra='allow'`` is deliberate: orchestration metadata and method
    options must survive descriptor normalisation and resume fingerprinting.
    """

    model_config = ConfigDict(extra="allow")


def _extension_uri(
    scheme: str,
    study_id: str,
    name: str,
    target_branch: Optional[str] = None,
) -> str:
    """Build a namespaced extension URI when a multi-target branch is set."""

    branch = str(target_branch or "").strip().strip("/")
    path = f"{branch}/{name}" if branch else name
    return f"{scheme}://{study_id}/{path}/v1"


def _extension_plan(capability_name: str, **context: Any) -> Dict[str, Any]:
    """Return the output/operation skeleton consumed by ``extend_plan``.

    The factory is intentionally data-only and side-effect free.  It receives
    keyword context so callers can pass target branch information without
    changing the legacy planner signature.
    """

    study_id = str(context.get("study_id") or context.get("sid") or "study")
    target_branch = context.get("target_branch")
    outputs = {
        "functional_activity": [_extension_uri("table", study_id, "functional_activity", target_branch)],
        "donor_sensitivity": [
            _extension_uri("table", study_id, "donor_sensitivity", target_branch),
            _extension_uri("table", study_id, "donor_sensitivity_summary", target_branch),
        ],
        "doublet_detection": [
            _extension_uri("adata", study_id, "doublet_checked", target_branch),
            _extension_uri("table", study_id, "doublet_scores", target_branch),
        ],
        "cell_annotation": [_extension_uri("adata", study_id, "reference_annotated", target_branch)],
        "background_removal": [
            _extension_uri("adata", study_id, "background_corrected", target_branch),
            _extension_uri("json", study_id, "cellbender_report", target_branch),
        ],
        "liana_communication": [
            _extension_uri("table", study_id, "liana_communication", target_branch),
            _extension_uri("table", study_id, "liana_condition_comparison", target_branch),
        ],
    }
    try:
        return {"expected_outputs": list(outputs[capability_name])}
    except KeyError as exc:
        raise ValueError(f"No advanced extension plan is declared for {capability_name!r}") from exc


def _descriptor(
    capability_name: str,
    method: str,
    *,
    input_types,
    output_types,
    required_audit_ids,
    method_aliases=(),
    description: str,
) -> CapabilityDescriptor:
    return CapabilityDescriptor(
        capability_name=capability_name,
        method=method,
        parameter_model=AdvancedExtensionParameters,
        input_types=input_types,
        output_types=output_types,
        required_audit_ids=required_audit_ids,
        method_aliases=method_aliases,
        plan_factory=lambda **context: _extension_plan(capability_name, **context),
        description=description,
    )


ADVANCED_EXTENSION_DESCRIPTORS: Dict[str, CapabilityDescriptor] = {
    "functional_activity": _descriptor(
        "functional_activity",
        "decoupler_ulm_v2",
        input_types=(ArtifactType.ANNDATA, ArtifactType.TABLE),
        output_types=(ArtifactType.TABLE,),
        required_audit_ids=("advanced_statistics_integrity",),
        method_aliases=("decoupler_functional_v1",),
        description="Donor-level decoupler ULM activity and condition comparison.",
    ),
    "donor_sensitivity": _descriptor(
        "donor_sensitivity",
        "pydeseq2_leave_one_donor_out_v1",
        input_types=(ArtifactType.ANNDATA,),
        output_types=(ArtifactType.TABLE,),
        required_audit_ids=("advanced_statistics_integrity",),
        method_aliases=("donor_loo_pydeseq2_v1",),
        description="Leave-one-donor-out PyDESeq2 sensitivity analysis.",
    ),
    "doublet_detection": _descriptor(
        "doublet_detection",
        "scanpy_scrublet_v1",
        input_types=(ArtifactType.ANNDATA,),
        output_types=(ArtifactType.ANNDATA, ArtifactType.TABLE, ArtifactType.JSON),
        required_audit_ids=("advanced_qc_integrity",),
        method_aliases=("scrublet", "scrublet_v1"),
        description="Scanpy Scrublet doublet scoring with optional filtering.",
    ),
    "cell_annotation": _descriptor(
        "cell_annotation",
        "celltypist_local_v1",
        input_types=(ArtifactType.ANNDATA,),
        output_types=(ArtifactType.ANNDATA, ArtifactType.TABLE, ArtifactType.JSON),
        required_audit_ids=("advanced_qc_integrity",),
        method_aliases=("celltypist",),
        description="Local CellTypist annotation with conflict provenance.",
    ),
    "background_removal": _descriptor(
        "background_removal",
        "cellbender_cli_v1",
        input_types=(ArtifactType.ANNDATA,),
        output_types=(ArtifactType.ANNDATA, ArtifactType.JSON),
        required_audit_ids=("advanced_qc_integrity",),
        method_aliases=("cellbender",),
        description="Explicit CellBender external background-removal adapter.",
    ),
    "liana_communication": _descriptor(
        "liana_communication",
        "liana_rank_aggregate_v1",
        input_types=(ArtifactType.ANNDATA, ArtifactType.SPATIAL_DATA),
        output_types=(ArtifactType.TABLE,),
        required_audit_ids=("liana_integrity",),
        description="Donor-condition LIANA ligand/receptor rank aggregation.",
    ),
}


def descriptor_for_capability(capability_name: str, method: Optional[str] = None) -> Optional[CapabilityDescriptor]:
    """Resolve a shared advanced descriptor by capability and method/alias."""

    descriptor = ADVANCED_EXTENSION_DESCRIPTORS.get(capability_name)
    if descriptor is None:
        return None
    if method is None or method == descriptor.method or method in descriptor.method_aliases:
        return descriptor
    return None


def advanced_extension_names() -> tuple[str, ...]:
    """Return the declared extension names in deterministic order."""

    return tuple(ADVANCED_EXTENSION_DESCRIPTORS)


def advanced_method_ids() -> frozenset[str]:
    """Return canonical advanced method IDs for independent auditors."""

    return frozenset(descriptor.method for descriptor in ADVANCED_EXTENSION_DESCRIPTORS.values())


__all__ = [
    "AdvancedExtensionParameters",
    "ADVANCED_EXTENSION_DESCRIPTORS",
    "descriptor_for_capability",
    "advanced_extension_names",
    "advanced_method_ids",
]
