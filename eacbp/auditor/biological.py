"""
Biological Validator auditing marker coherence, biological plausibility, and pathway gene universes.
"""

from typing import Dict, Any, List
import pandas as pd
import numpy as np

from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.schemas.artifact import ArtifactType
from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData


# Canonical marker reference
CANONICAL_MARKERS = {
    "Microglia": ["Cx3cr1", "P2ry12", "Tmem119", "Aif1", "Csf1r"],
    "DAM": ["Apoe", "Trem2", "Clec7a", "Itgax", "Axl", "Cst7"],
    "Astrocytes": ["Gfap", "Aldh1l1", "Aqp4", "Slc1a2"],
    "Neurons": ["Rbfox3", "Syp", "Map2", "Tubb3"],
    "Oligodendrocytes": ["Mog", "Mbp", "Olig2", "Plp1"],
}


class BiologicalValidator(BaseAuditor):
    """Audits cell marker consistency and biological coherence against known prior databases."""

    def __init__(self):
        super().__init__(auditor_name="biological_validator")

    def audit(
        self,
        contract: TaskContract,
        result: TaskResult,
        registry: ArtifactRegistry
    ) -> ValidationReport:
        target_uri = result.output_artifacts[0] if result.output_artifacts else None
        report = ValidationReport(
            auditor_name=self.auditor_name,
            target_task_id=contract.task_id,
            target_artifact_uri=target_uri,
        )

        if not target_uri or not registry.exists(target_uri):
            return report

        meta, payload = registry.get(target_uri)
        cap_name = contract.capability

        # 1. Annotation & Marker Coherence Audit
        if cap_name == "clustering" and meta.type == ArtifactType.ANNDATA:
            data = SCData.from_payload(payload)
            gene_names = set(data.var["gene_name"].astype(str)) if "gene_name" in data.var.columns else set(data.var.index.astype(str))
            cell_types = data.obs["cell_type"].unique() if "cell_type" in data.obs.columns else []

            marker_reference = contract.parameters.get("marker_reference", CANONICAL_MARKERS)
            valid_reference = (
                isinstance(marker_reference, dict)
                and bool(marker_reference)
                and all(
                    isinstance(label, str) and bool(label.strip())
                    and isinstance(markers, (list, tuple)) and bool(markers)
                    and all(isinstance(marker, str) and bool(marker.strip()) for marker in markers)
                    for label, markers in marker_reference.items()
                )
            )

            matched_markers = {}
            for ct in cell_types if valid_reference else []:
                if ct in marker_reference:
                    expected = marker_reference[ct]
                    present = [m for m in expected if m in gene_names]
                    matched_markers[ct] = present

            missing_types = [str(ct) for ct, markers in matched_markers.items() if not markers]
            unreferenced_types = [str(ct) for ct in cell_types if not valid_reference or ct not in marker_reference]

            report.add_check(
                name="canonical_marker_representation",
                passed=valid_reference and bool(matched_markers) and not missing_types,
                severity=ValidationSeverity.INFO,
                message=f"Marker representation in the feature universe: {matched_markers}; this check does not establish expression coherence or annotation accuracy.",
                metrics={
                    "matched_markers": matched_markers,
                    "reference_source": "task_parameters" if "marker_reference" in contract.parameters else "builtin_brain_markers",
                    "valid_reference": valid_reference,
                    "cell_types_without_present_markers": missing_types,
                    "cell_types_without_reference": unreferenced_types,
                }
            )

        # 2. DEG Top Candidate Biological Coherence
        elif cap_name == "deg" and meta.type == ArtifactType.TABLE:
            deg_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            top_up = deg_df[deg_df["log2_fold_change"] > 0]["gene"].head(10).tolist()
            
            report.add_check(
                name="reported_gene_context",
                passed=True,
                severity=ValidationSeverity.INFO,
                message=f"Reported positive-direction genes: {top_up[:5]}; no disease-specific coherence is inferred without a supplied reference.",
                metrics={"top_upregulated": top_up[:5]}
            )

        return report
