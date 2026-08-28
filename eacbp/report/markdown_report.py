"""
Scientific Report Generator producing 4-tier calibrated reports with interactive sentence-to-artifact provenance links,
multimodal DAG provenance, independent auditor sign-offs, and mandatory epistemic tagging.
"""

from typing import List, Dict, Any, Optional
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.evidence import ClaimNode, LanguageTier
from eacbp.evidence.graph import EvidenceGraph
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.report.provenance import SentenceProvenanceTracker
from eacbp.auditor.base import ValidationReport


class ScientificReportGenerator:
    """Generates rigorous scientific manuscripts with strict 4-Tier Language and clickable provenance cards."""

    def __init__(
        self,
        manifest: StudyManifest,
        evidence_graph: EvidenceGraph,
        artifact_registry: ArtifactRegistry,
        task_history: List[Any],
        audit_reports: Optional[List[ValidationReport]] = None,
    ):
        self.manifest = manifest
        self.evidence_graph = evidence_graph
        self.artifact_registry = artifact_registry
        self.task_history = task_history
        self.audit_reports = audit_reports or []
        self.tracker = SentenceProvenanceTracker(evidence_graph, artifact_registry, task_history)

    def generate_markdown(self) -> str:
        manifest = self.manifest
        lines = []

        is_prior_guided = (
            manifest.analysis_policy.prior_guided_analysis
            or bool(manifest.hypotheses.user_provided)
            or any("[PRIOR-GUIDED" in c.statement for c in self.evidence_graph.claim_nodes.values())
        )

        # Title & Study Header
        lines.append(f"# Scientific Study Report: {manifest.title}")
        lines.append(
            f"**Study ID**: `{manifest.study_id}` | "
            f"**Species**: *{manifest.biological_design.species}* | "
            f"**Tissue**: *{manifest.biological_design.tissue}* | "
            f"**Disease**: *{manifest.biological_design.disease}*"
        )
        lines.append("")
        lines.append("> [!NOTE]")
        lines.append(
            "> **EACBP Epistemic Governance**: All claims are strictly categorized across 4 calibrated language tiers "
            "and backed by verified computational artifacts, statistical audits, and multimodal biological evidence nodes."
        )
        lines.append("")

        # Prior-Guided Epistemic Callout Banner
        if is_prior_guided:
            lines.append("> [!WARNING]")
            lines.append("> **[PRIOR-GUIDED HYPOTHESIS TESTING]**")
            lines.append(
                "> This study includes prior-guided hypothesis evaluation testing targeted biological axes. "
                "Claims tagged with prior guidance reflect confirmatory hypothesis evaluation and must not be conflated with unbiased exploratory discovery."
            )
            lines.append("")

        # Section 1: Study Design & Experimental Audit
        lines.append("## 1. Study Design & Experimental Audit")
        lines.append(
            f"- **Biological Replication Units**: {manifest.experimental_design.total_samples} "
            f"{manifest.experimental_design.biological_unit}s ({manifest.experimental_design.donor_replicates_per_condition})"
        )
        lines.append(f"- **Batches**: {', '.join(manifest.experimental_design.batches) if manifest.experimental_design.batches else 'None'}")
        lines.append(f"- **Modalities**: {', '.join(manifest.data.modalities)}")
        lines.append(
            f"- **Analysis Policy**: Discovery Mode: `{manifest.analysis_policy.discovery_mode}`, "
            f"Prior-Guided Analysis: `{manifest.analysis_policy.prior_guided_analysis}`, "
            f"Pseudobulk Preferred: `{manifest.analysis_policy.prefer_pseudobulk}`"
        )
        if manifest.hypotheses.user_provided:
            lines.append(f"- **Prior Hypotheses Evaluated**: {', '.join(manifest.hypotheses.user_provided)}")
        lines.append("")

        # Section 2: Computational Task DAG Execution Summary
        lines.append("## 2. Computational Task DAG Execution Summary")
        lines.append("| Task ID | Capability | Method Used | Status | Execution Time | Output Artifacts |")
        lines.append("| :--- | :--- | :--- | :--- | :--- | :--- |")
        for task_res in self.task_history:
            t_id = getattr(task_res, "task_id", "unknown")
            cap = getattr(task_res, "capability", "unknown")
            meth = getattr(task_res, "method_used", "unknown")
            status = getattr(task_res, "status", "unknown")
            status_val = status.value if hasattr(status, "value") else str(status)
            exec_time = f"{getattr(task_res, 'execution_time_sec', 0.0):.2f}s"
            outs = "<br>".join([f"`{u}`" for u in getattr(task_res, "output_artifacts", [])])
            lines.append(f"| `{t_id}` | `{cap}` | `{meth}` | `{status_val}` | {exec_time} | {outs} |")
        lines.append("")

        # Section 3: Evidence-Grounded Multimodal Claims
        lines.append("## 3. Evidence-Grounded Scientific Claims & Multimodal Calibration")
        lines.append("")
        for claim in self.evidence_graph.claim_nodes.values():
            card = self.tracker.resolve_claim_provenance_card(claim.claim_id)
            conf = claim.confidence
            
            lines.append(f"### Claim `{claim.claim_id}`: {claim.statement}")
            lines.append(f"- **Language Tier**: `{claim.language_tier.value}`")
            lines.append(f"- **Causal Status**: `{claim.causal_status}`")
            lines.append(
                f"- **Calibrated Confidence**: Association: `{conf.association:.2f}` | "
                f"Mechanistic: `{conf.mechanistic:.2f}` | "
                f"Causal: `{conf.causal:.2f}` | "
                f"**Overall: `{conf.overall:.2f}`**"
            )
            lines.append("")
            lines.append("#### Supporting Evidence & Sentence-Level Provenance Trace:")
            for ev in card["evidence_chain"]:
                task_info = ev["source_task"]
                lines.append(f"- **Evidence `{ev['evidence_id']}`** (`{ev['evidence_type']}`, `{ev['strength']}`, score: `{ev['score']:.2f}`): {ev['summary']}")
                lines.append(f"  - *Generated by Task*: `{task_info['task_id']}` (`{task_info['capability']}` -> `{task_info['method']}`)")
                for art in ev["artifacts"]:
                    lines.append(f"  - *Source Artifact*: `{art['uri']}` (SHA-256: `{art['sha256'][:12]}...`)")
                    if art.get("lineage_path"):
                        lines.append(f"  - *Lineage Path*: `{' -> '.join(art['lineage_path'])}`")
            lines.append("")

        # Section 4: Scientific Auditor Sign-offs & Independent Verification
        if self.audit_reports:
            lines.append("## 4. Scientific Auditor Sign-offs & Independent Verification")
            lines.append("| Target Task ID | Target Artifact | Overall Audit Status | Verification Checks Passed |")
            lines.append("| :--- | :--- | :--- | :--- |")
            for rep in self.audit_reports:
                t_id = rep.target_task_id
                t_art = rep.target_artifact_uri or "N/A"
                status_str = "PASSED" if rep.overall_passed else "FLAGGED"
                pass_count = sum(1 for c in rep.checks if c.passed)
                total_count = len(rep.checks)
                lines.append(f"| `{t_id}` | `{t_art}` | `{status_str}` | {pass_count}/{total_count} checks passed |")
            lines.append("")

        # Section 5: Registered Artifacts & Lineage DAG
        lines.append("## 5. Registered Artifacts & Lineage DAG")
        lines.append("| Artifact URI | Type | Operation | Task ID | Size (Bytes) | SHA-256 Checksum |")
        lines.append("| :--- | :--- | :--- | :--- | :--- | :--- |")
        for uri, meta in self.artifact_registry.registry.items():
            lines.append(
                f"| `{meta.uri}` | `{meta.type.value}` | `{meta.operation}` | "
                f"`{meta.created_by_task}` | {meta.size_bytes} | `{meta.sha256_hash[:16]}...` |"
            )
        lines.append("")

        # Section 6: Lineage Mermaid Diagram
        lines.append("## 6. Computational Lineage Graph")
        lines.append("```mermaid")
        lines.append(self.artifact_registry.lineage.to_mermaid())
        lines.append("```")
        lines.append("")

        # Section 7: Evidence-to-Claim Mermaid Graph
        lines.append("## 7. Evidence-to-Claim DAG")
        lines.append("```mermaid")
        lines.append(self.evidence_graph.to_mermaid())
        lines.append("```")
        lines.append("")

        return "\n".join(lines)
