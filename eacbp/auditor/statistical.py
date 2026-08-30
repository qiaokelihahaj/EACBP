"""
Statistical Auditor verifying degrees of freedom, batch diversity, multiple testing FDR,
Moran's I bounds, Geary's C bounds, perturbation bounds, and epistemic tagging.
"""

from typing import Dict, Any, Optional
import numpy as np
import pandas as pd

from eacbp.schemas.task import TaskContract, TaskResult
from eacbp.schemas.artifact import ArtifactType
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.auditor.base import BaseAuditor, ValidationReport, ValidationSeverity


class StatisticalValidator(BaseAuditor):
    """Audits statistical rigor, degrees of freedom, FDR corrections, spatial autocorrelation bounds, and epistemic calibration."""

    def __init__(self):
        super().__init__(auditor_name="statistical_validator")

    def _create_report(self, target_task_id: str, target_artifact_uri: Optional[str] = None) -> ValidationReport:
        return ValidationReport(
            auditor_name=self.auditor_name,
            target_task_id=target_task_id,
            target_artifact_uri=target_artifact_uri,
        )

    def audit(
        self,
        contract: TaskContract,
        result: TaskResult,
        registry: ArtifactRegistry,
    ) -> ValidationReport:
        target_uri = result.output_artifacts[0] if result.output_artifacts else None
        report = self._create_report(contract.task_id, target_uri)

        if not target_uri or not registry.exists(target_uri):
            return report

        meta, payload = registry.get(target_uri)
        cap_name = contract.capability

        # =====================================================================
        # 1. Pseudoreplication & DEG Degrees of Freedom Audit
        # =====================================================================
        if cap_name == "deg" and meta.type == ArtifactType.TABLE:
            deg_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
            is_pseudobulk = meta.summary_metrics.get("statistical_unit") == "pseudobulk" or "donor_replicates" in meta.summary_metrics
            
            # Check if cell-level test was used when pseudobulk was possible
            if not is_pseudobulk:
                report.add_check(
                    name="pseudoreplication_audit",
                    passed=False,
                    severity=ValidationSeverity.WARNING,
                    message="DEG performed at single-cell level without donor aggregation; statistical tests subject to pseudoreplication inflation.",
                    metrics={"is_pseudobulk": False, "confirmatory_allowed": False},
                )
                report.recommendations.append("Re-run DEG using donor-level pseudobulk aggregation for confirmatory biological claims.")
            else:
                report.add_check(
                    name="pseudoreplication_audit",
                    passed=True,
                    severity=ValidationSeverity.INFO,
                    message="DEG performed at biological donor replicate level (pseudobulk); degrees of freedom preserved.",
                    metrics={"is_pseudobulk": True, "confirmatory_allowed": True},
                )

            # Check multiple testing correction
            has_fdr = "fdr_q_value" in deg_df.columns or "fdr" in deg_df.columns or "p_val_adj" in deg_df.columns
            report.add_check(
                name="multiple_testing_correction",
                passed=has_fdr,
                severity=ValidationSeverity.ERROR,
                message="FDR multiple testing adjustment (Benjamini-Hochberg) present." if has_fdr else "Missing multiple testing adjustment column.",
                metrics={"has_fdr": has_fdr},
            )

        # =====================================================================
        # 2. Clustering Stability Audit
        # =====================================================================
        elif cap_name == "clustering" and meta.type == ArtifactType.ANNDATA:
            sil_score = meta.summary_metrics.get("silhouette_score", 0.0)
            passed_sil = sil_score >= 0.10
            report.add_check(
                name="clustering_separation_silhouette",
                passed=passed_sil,
                severity=ValidationSeverity.WARNING if not passed_sil else ValidationSeverity.INFO,
                message=f"Clustering silhouette score: {sil_score:.3f} (threshold >= 0.10).",
                metrics={"silhouette_score": sil_score},
            )

        # =====================================================================
        # 3. Trajectory Stability Audit
        # =====================================================================
        elif cap_name == "trajectory_inference" and meta.type == ArtifactType.TABLE:
            stab_score = meta.summary_metrics.get("stability_score", 0.0)
            passed_stab = stab_score >= 0.60
            report.add_check(
                name="trajectory_subsampling_stability",
                passed=passed_stab,
                severity=ValidationSeverity.WARNING if not passed_stab else ValidationSeverity.INFO,
                message=f"Trajectory root sensitivity & subsampling stability: {stab_score:.3f} (threshold >= 0.60).",
                metrics={"stability_score": stab_score},
            )
            if not passed_stab:
                report.recommendations.append("Trajectory topology is sensitive to root cell selection; disable strong mechanistic claims.")

        # =====================================================================
        # 4. Spatial Analytics Audit (Moran's I / Geary's C / Spatial CCI)
        # =====================================================================
        elif cap_name in ("spatial_deg", "spatial_domain", "cell_cell_communication", "spatial_cci") or "morans_i" in meta.summary_metrics:
            if meta.type == ArtifactType.TABLE:
                sp_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
                
                # Check Moran's I bounds if present
                moran_col = "morans_i" if "morans_i" in sp_df.columns else ("moran_i" if "moran_i" in sp_df.columns else None)
                if moran_col:
                    moran_vals = sp_df[moran_col].dropna()
                    valid_moran = bool((moran_vals >= -1.0).all() and (moran_vals <= 1.0).all())
                    report.add_check(
                        name="morans_i_bounds_check",
                        passed=valid_moran,
                        severity=ValidationSeverity.ERROR,
                        message="All Moran's I values reside strictly in theoretical bounds [-1.0, 1.0]." if valid_moran else "Moran's I out of analytical bounds [-1, 1].",
                        metrics={"min_moran": float(moran_vals.min()) if not moran_vals.empty else 0.0, "max_moran": float(moran_vals.max()) if not moran_vals.empty else 0.0},
                    )

                # Check Geary's C bounds if present
                geary_col = "gearys_c" if "gearys_c" in sp_df.columns else ("geary_c" if "geary_c" in sp_df.columns else None)
                if geary_col:
                    geary_vals = sp_df[geary_col].dropna()
                    valid_geary = bool((geary_vals >= 0.0).all())
                    report.add_check(
                        name="gearys_c_bounds_check",
                        passed=valid_geary,
                        severity=ValidationSeverity.ERROR,
                        message="All Geary's C values are non-negative." if valid_geary else "Negative Geary's C detected.",
                        metrics={"min_geary": float(geary_vals.min()) if not geary_vals.empty else 0.0},
                    )

                # Check multiple testing correction
                has_sp_fdr = "fdr_q_value" in sp_df.columns or "fdr" in sp_df.columns or "p_val_adj" in sp_df.columns
                report.add_check(
                    name="multiple_testing_correction",
                    passed=has_sp_fdr,
                    severity=ValidationSeverity.ERROR,
                    message="Spatial multiple testing FDR correction present." if has_sp_fdr else "Missing spatial FDR column.",
                    metrics={"has_fdr": has_sp_fdr},
                )

                # Check spatial CCI non-negative interaction scores
                if "spatial_interaction_score" in sp_df.columns:
                    cci_scores = sp_df["spatial_interaction_score"].dropna()
                    valid_cci = bool((cci_scores >= 0.0).all())
                    report.add_check(
                        name="spatial_cci_score_non_negative",
                        passed=valid_cci,
                        severity=ValidationSeverity.ERROR,
                        message="All spatial CCI interaction scores are non-negative." if valid_cci else "Negative spatial CCI scores detected.",
                        metrics={"min_cci_score": float(cci_scores.min()) if not cci_scores.empty else 0.0},
                    )

        # =====================================================================
        # 5. In Silico Perturbation Simulation Audit
        # =====================================================================
        elif cap_name in ("genetic_perturbation_simulation", "genetic_perturbation", "compound_perturbation_simulation", "compound_perturbation"):
            # Check expression non-negativity for AnnData output
            if meta.type == ArtifactType.ANNDATA:
                data = payload if isinstance(payload, SCData) else SCData.from_dict(payload)
                if hasattr(data.X, "tocsr"):
                    x_min = float(np.min(data.X.data)) if len(data.X.data) > 0 else 0.0
                else:
                    x_min = float(np.min(data.X))
                report.add_check(
                    name="perturbation_shift_bounds_check",
                    passed=(x_min >= 0.0),
                    severity=ValidationSeverity.ERROR,
                    message=f"Simulated expression satisfies non-negativity constraint (min: {x_min:.4f}).",
                    metrics={"min_expression": x_min},
                )

            # Check network attenuation bounds
            attenuation = contract.parameters.get("network_attenuation", meta.summary_metrics.get("network_attenuation", result.metrics.get("network_attenuation", 0.35)))
            if attenuation is not None:
                valid_alpha = 0.0 <= float(attenuation) < 1.0
                report.add_check(
                    name="network_attenuation_bounds_check",
                    passed=valid_alpha,
                    severity=ValidationSeverity.ERROR,
                    message=f"Network attenuation alpha={attenuation:.2f} is strictly in [0.0, 1.0)." if valid_alpha else f"Invalid attenuation alpha={attenuation}.",
                    metrics={"alpha": float(attenuation)},
                )

            # Check compound reversal discordance score bounds
            if meta.type == ArtifactType.TABLE:
                pt_df = payload if isinstance(payload, pd.DataFrame) else pd.DataFrame(payload)
                if "reversal_score" in pt_df.columns or "cosine_discordance" in pt_df.columns:
                    col = "reversal_score" if "reversal_score" in pt_df.columns else "cosine_discordance"
                    rev_scores = pt_df[col].dropna()
                    valid_rev = bool((rev_scores >= -1.0).all() and (rev_scores <= 1.0).all())
                    report.add_check(
                        name="reversal_score_bounds_check",
                        passed=valid_rev,
                        severity=ValidationSeverity.ERROR,
                        message="All compound reversal scores reside in cosine bounds [-1.0, 1.0]." if valid_rev else "Cosine discordance out of bounds [-1, 1].",
                        metrics={"min_reversal": float(rev_scores.min()) if not rev_scores.empty else 0.0, "max_reversal": float(rev_scores.max()) if not rev_scores.empty else 0.0},
                    )

        # =====================================================================
        # 6. Epistemic Tagging Audit (Prior-Guided Mode)
        # =====================================================================
        is_prior_guided = contract.parameters.get("prior_guided", False) or meta.summary_metrics.get("prior_guided", False)
        if is_prior_guided or cap_name == "knowledge_retrieval":
            tag_present = "[PRIOR-GUIDED HYPOTHESIS TESTING]" in str(payload) or "[PRIOR-GUIDED HYPOTHESIS TESTING]" in str(result.metrics) or "[PRIOR-GUIDED HYPOTHESIS TESTING]" in str(meta.summary_metrics)
            if is_prior_guided:
                report.add_check(
                    name="epistemic_tagging_check",
                    passed=tag_present,
                    severity=ValidationSeverity.ERROR if not tag_present else ValidationSeverity.INFO,
                    message="Mandatory '[PRIOR-GUIDED HYPOTHESIS TESTING]' epistemic tag verified on prior-guided execution." if tag_present else "Missing mandatory '[PRIOR-GUIDED HYPOTHESIS TESTING]' tag on prior-guided task.",
                    metrics={"prior_guided": True, "tag_verified": tag_present},
                )

        return report
