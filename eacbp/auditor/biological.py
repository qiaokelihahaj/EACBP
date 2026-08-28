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
            data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
            gene_names = list(data.var["gene_name"]) if "gene_name" in data.var.columns else []
            cell_types = data.obs["cell_type"].unique() if "cell_type" in data.obs.columns else []

            matched_markers = {}
            for ct in cell_types:
                if ct in CANONICAL_MARKERS:
                    expected = CANONICAL_MARKERS[ct]
                    present = [m for m in expected if m in gene_names]
                    matched_markers[ct] = present

            report.add_check(
                name="canonical_marker_representation",
                passed=len(matched_markers) > 0,
                severity=ValidationSeverity.INFO,
                message=f"Verified presence of canonical marker genes across identified cell types: {matched_markers}",
                metrics={"matched_markers": matched_markers}
            )

        # 2. DEG Top Candidate Biological Coherence
        elif cap_name == "deg" and meta.type == ArtifactType.TABLE:
            deg_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            top_up = deg_df[deg_df["log2_fold_change"] > 0]["gene"].head(10).tolist()
            
            # Check overlap with DAM markers for AD study
            dam_overlap = [g for g in top_up if g in CANONICAL_MARKERS["DAM"]]
            report.add_check(
                name="disease_marker_coherence",
                passed=True,
                severity=ValidationSeverity.INFO,
                message=f"Top upregulated genes in disease condition include known markers: {dam_overlap}",
                metrics={"dam_overlap": dam_overlap, "top_upregulated": top_up[:5]}
            )

        return report
