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
from eacbp.orchestrator.dag import target_branch_slug_map


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
            "> **EACBP Evidence Policy**: Claims restate admitted evidence; failed audits are excluded. "
            "Confidence scores are heuristic summaries, not calibrated probabilities. Local curated knowledge is unverified context."
        )
        lines.append("")

        simulated = any(c.is_simulated for c in self.evidence_graph.claim_nodes.values()) or any(
            a.summary_metrics.get("is_simulated", False)
            for a in self.artifact_registry.list_artifacts(study_id=manifest.study_id))
        if simulated:
            lines.extend(["> [!WARNING]", "> **SIMULATED DATA — demonstration outputs; not measurements from a real biological study.**", ""])
        rejected = [t for t in self.task_history if getattr(getattr(t, "status", None), "value", "") not in ("success",)]
        if rejected:
            lines.extend(["> [!WARNING]", "> **INCOMPLETE STUDY**: Some tasks failed audit, execution, or were blocked. Only admitted upstream observations are reported.", ""])

        # Prior-Guided Epistemic Callout Banner
        if is_prior_guided:
            lines.append("> [!WARNING]")
            lines.append("> **[PRIOR-GUIDED HYPOTHESIS TESTING]**")
            lines.append(
                "> This study includes prior-guided hypothesis evaluation testing targeted biological axes. "
                "Prior guidance records the question being explored; it does not establish confirmation or independently verified knowledge."
            )
            lines.append("")

        # Section 1: Study Design & Experimental Audit
        lines.append("## 1. Study Design & Experimental Audit")
        observed_design = next((t.metrics for t in self.task_history if getattr(t, "capability", None) == "dataset_audit" and getattr(t.status, "value", "") == "success"), {})
        samples = observed_design.get("total_samples", manifest.experimental_design.total_samples)
        donors = observed_design.get("donor_replicates", manifest.experimental_design.donor_replicates_per_condition)
        batches = observed_design.get("batches", manifest.experimental_design.batches)
        lines.append(f"- **Biological Replication Units**: {samples or 'unknown'} ({donors or 'unknown'})")
        lines.append(f"- **Batches**: {', '.join(map(str, batches)) if batches else 'unknown'}")
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
        decisions = observed_design.get("planning_decisions", [])
        if decisions:
            lines.append("### Analysis branches omitted after metadata audit")
            lines.extend(f"- `{item['task_id']}`: {item['reason']}" for item in decisions)
            lines.append("")
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

        # Multi-target runs have one shared preprocessing prefix and one
        # independent target branch for population-dependent statistics and
        # consumers.  Keep this grouping explicit in the report so a reader
        # can follow a target's provenance without inferring it from URI text.
        targets = list(manifest.biological_design.target_cell_types)
        if len(targets) > 1:
            slug_by_target = target_branch_slug_map(targets)
            lines.extend([
                "### Target branch execution groups",
                "",
                "Preprocessing is shared through annotation. Subset, abundance, DEG, trajectory, and their consumers are isolated per target. Multiple-testing correction remains specific to each analysis within a target; there is no pooled correction across analyses or targets.",
                "",
            ])
            for target in targets:
                slug = slug_by_target[target]
                branch_tasks = [
                    task for task in self.task_history
                    if self._task_target(task, slug_by_target) == (target, slug)
                ]
                lines.extend([f"#### Target: `{target}` (branch `{slug}`)", ""])
                if not branch_tasks:
                    lines.append("No target-specific task result was recorded.")
                    lines.append("")
                    continue
                lines.extend(["| Task ID | Capability | Status | Target provenance | Outputs |", "| :--- | :--- | :--- | :--- | :--- |"])
                for task in branch_tasks:
                    status = getattr(task, "status", "unknown")
                    status_val = status.value if hasattr(status, "value") else str(status)
                    outputs = "<br>".join(f"`{uri}`" for uri in getattr(task, "output_artifacts", []))
                    lines.append(
                        f"| `{getattr(task, 'task_id', 'unknown')}` | `{getattr(task, 'capability', 'unknown')}` | `{status_val}` | `{target}` / `{slug}` | {outputs} |"
                    )
                lines.append("")

        # Section 3: Evidence-Grounded Multimodal Claims
        lines.extend(self._advanced_analysis_details())
        lines.append("## 3. Evidence-Grounded Scientific Claims & Multimodal Calibration")
        lines.append("")
        for claim in self.evidence_graph.claim_nodes.values():
            card = self.tracker.resolve_claim_provenance_card(claim.claim_id)
            conf = claim.confidence
            
            lines.append(f"### Claim `{claim.claim_id}`: {claim.statement}")
            lines.append(f"- **Language Tier**: `{claim.language_tier.value}`")
            lines.append(f"- **Causal Status**: `{claim.causal_status}`")
            lines.append(
                f"- **Heuristic Confidence (not a probability)**: Association: `{conf.association:.2f}` | "
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
                source_result = next((item for item in self.task_history if getattr(item, "task_id", None) == task_info.get("task_id")), None)
                target_info = self._task_target(source_result, target_branch_slug_map(list(manifest.biological_design.target_cell_types))) if source_result is not None and len(manifest.biological_design.target_cell_types) > 1 else None
                if target_info is not None:
                    lines.append(f"  - *Target branch*: `{target_info[0]}` (`{target_info[1]}`); FDR family: `target_cell_type:{target_info[0]}`")
                for art in ev["artifacts"]:
                    lines.append(f"  - *Source Artifact*: `{art['uri']}` (SHA-256: `{art['sha256'][:12]}...`)")
                    if art.get("lineage_path"):
                        lines.append(f"  - *Lineage Path*: `{' -> '.join(art['lineage_path'])}`")
                        if len(art["ancestor_artifacts"]) > len(art["lineage_path"]):
                            lines.append(f"  - *All contributing artifacts*: `{', '.join(art['ancestor_artifacts'])}`")
            lines.append("")

        lines.extend([
            "### Evidence dependence and validation scope",
            "Analyses derived from shared input artifacts are dependent evidence. DEG, inferred functional activity and communication from the same matrix do not establish independent replication or experimental causality.",
            "Confidence dimensions conservatively collapse shared data origins; sensitivity analyses do not contribute additional independent support. Distinct artifact roots alone do not prove independent cohorts.",
            "For multi-target studies, multiple-testing correction remains specific to each analysis within a target. FDR values are not pooled across analyses or target branches.",
            "",
            "| Evidence | Root input artifacts | Validation scope |",
            "| :--- | :--- | :--- |",
        ])
        for node in self.evidence_graph.evidence_nodes.values():
            origins = ", ".join(f"`{uri}`" for uri in node.data_origin_uris) or "Unresolved; independence not established"
            lines.append(f"| `{node.evidence_id}` | {origins} | `{node.validation_scope}` |")
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
        study_artifacts = self.artifact_registry.list_artifacts(study_id=manifest.study_id)
        for meta in study_artifacts:
            lines.append(
                f"| `{meta.uri}` | `{meta.type.value}` | `{meta.operation}` | "
                f"`{meta.created_by_task}` | {meta.size_bytes} | `{meta.sha256_hash[:16]}...` |"
            )
        lines.append("")

        # Section 6: Lineage Mermaid Diagram
        lines.append("## 6. Computational Lineage Graph")
        lines.append("```mermaid")
        from eacbp.artifact.lineage import LineageGraph
        relevant_lineage = LineageGraph()
        selected_uris = {a.uri for a in study_artifacts}
        for uri in list(selected_uris):
            selected_uris.update(self.artifact_registry.lineage.get_ancestors(uri))
        for uri in sorted(selected_uris):
            if self.artifact_registry.exists(uri):
                relevant_lineage.add_artifact(self.artifact_registry.get_metadata(uri))
        lines.append(relevant_lineage.to_mermaid())
        lines.append("```")
        lines.append("")

        # Section 7: Evidence-to-Claim Mermaid Graph
        lines.append("## 7. Evidence-to-Claim DAG")
        lines.append("```mermaid")
        lines.append(self.evidence_graph.to_mermaid())
        lines.append("```")
        lines.append("")

        return "\n".join(lines)

    def _advanced_analysis_details(self):
        """Include audited full-result context, not just significant claims."""
        import pandas as pd
        advanced = {"functional_activity", "donor_sensitivity", "doublet_detection",
                    "cell_annotation", "background_removal", "liana_communication"}
        selected = [t for t in self.task_history if t.capability in advanced or
                    (t.capability == "deg" and "pydeseq2" in t.method_used)]
        if not selected:
            return []
        lines = ["### Advanced analysis results and limitations", "",
                 "Full registered tables preserve nonsignificant and unestimated outcomes. Previews below are limited to 20 rows; missing estimates do not imply a zero effect.", ""]
        audited = {r.target_task_id for r in self.audit_reports if r.overall_passed and not r.stop_rule_triggered}
        columns = ["gene", "source", "target", "pathway", "activity", "condition_a", "condition_b",
                   "log2_fold_change", "lfc_se", "ci_low", "ci_high", "fdr_q_value", "activity_effect_condition_a_vs_b",
                   "n_donors_condition_a", "n_donors_condition_b", "status", "reason",
                   "direction_consistency", "min_log2_fold_change", "max_log2_fold_change",
                   "sender_cell_type", "receiver_cell_type", "ligand", "receptor",
                   "magnitude_rank", "specificity_rank", "donor_id", "condition",
                   "n_donors_a", "n_donors_b", "comparison_status",
                   "comparison_p_value_magnitude", "comparison_fdr_magnitude",
                   "comparison_p_value_specificity", "comparison_fdr_specificity",
                   "estimated_coverage", "significance_retention", "summary_status",
                   "n_estimated_loo", "baseline_log2_fold_change", "baseline_fdr_q_value",
                   "baseline_significant_fdr05", "n_significant_loo", "median_log2_fold_change",
                   "fit_status", "status_reason", "skip_reason", "n_requested_fits",
                   "n_successful_fits", "complete_fit_coverage",
                   "scientific_robustness_claim_supported"]

        metric_keys = {
            "doublet_detection": (
                "n_cells_before", "n_cells_after", "n_doublets_marked", "n_doublets_filtered",
                "retention_rate", "threshold", "threshold_by_batch", "counts_source",
                "filter_doublets",
            ),
            "cell_annotation": (
                "n_cells_before", "n_cells_after", "n_annotations", "n_unknown_annotations",
                "n_label_conflicts", "annotation_key", "use_as_cell_type",
                "existing_cell_type_preserved", "model_metadata",
            ),
            "background_removal": (
                "n_cells_before", "n_cells_after", "n_genes", "returncode", "quality_status",
                "quality_warnings", "scientific_quality_certified", "report_available",
                "metrics_available", "log_available", "corrected_counts_layer",
            ),
        }

        def cell(value):
            if value is None or (isinstance(value, (float, int)) and pd.isna(value)):
                return "not estimated"
            if isinstance(value, float):
                return f"{value:.4g}"
            return str(value).replace("|", "\\|").replace("\n", " ")

        def append_metric_table(task):
            keys = metric_keys.get(task.capability, ())
            present = [key for key in keys if key in task.metrics]
            if not present:
                return
            lines.extend(["", "Recorded task metrics:", "", "| Metric | Recorded value |", "| :--- | :--- |"])
            for key in present:
                lines.append(f"| `{key}` | {cell(task.metrics[key])} |")

        def append_liana_uncertainty(payload):
            fdr_columns = [
                key for key in ("comparison_fdr_magnitude", "comparison_fdr_specificity")
                if key in payload
            ]
            if not fdr_columns and "comparison_status" not in payload:
                return
            non_significant = 0
            not_estimated = 0
            for key in fdr_columns:
                values = pd.to_numeric(payload[key], errors="coerce")
                non_significant += int((values >= 0.05).sum())
                not_estimated += int(values.isna().sum())
            if "comparison_status" in payload:
                statuses = payload["comparison_status"].astype(str).str.casefold()
                not_estimated = max(not_estimated, int((statuses != "tested").sum()))
            estimated_entries = non_significant + not_estimated
            if fdr_columns:
                lines.append(
                    f"LIANA comparison FDR entries: {estimated_entries + int(sum(pd.to_numeric(payload[key], errors='coerce').lt(0.05).sum() for key in fdr_columns))}; "
                    f"non-significant (FDR >= 0.05): {non_significant}; not estimated: {not_estimated}."
                )
            else:
                lines.append(f"LIANA condition comparisons: non-significant status was not estimated; not estimated: {not_estimated}.")

        slug_by_target = target_branch_slug_map(list(self.manifest.biological_design.target_cell_types)) if len(self.manifest.biological_design.target_cell_types) > 1 else {}
        for task in selected:
            target_info = self._task_target(task, slug_by_target) if slug_by_target else None
            suffix = f" — Target: `{target_info[0]}` (branch `{target_info[1]}`)" if target_info else ""
            lines.extend([f"#### `{task.task_id}` — `{task.method_used}`{suffix}", ""])
            status_value = getattr(getattr(task, "status", None), "value", str(getattr(task, "status", "")))
            if status_value != "success" or task.task_id not in audited:
                lines.extend([f"Results not admitted: {cell(task.error_message or 'required audit did not pass')}", ""])
                continue
            for key in ("status", "reason", "skip_reason", "statistical_unit", "design_formula", "design", "interpretation_limit", "limitations"):
                if key in task.metrics:
                    lines.append(f"- {key}: {cell(task.metrics[key])}")
            append_metric_table(task)
            if task.capability == "doublet_detection":
                lines.append("Doublet scores and labels are QC observations; they do not add biological mechanism or causal support.")
            elif task.capability == "cell_annotation":
                lines.append("Unknown and conflicting CellTypist labels are retained as uncertainty; annotation does not add biological mechanism or causal support.")
            elif task.capability == "background_removal":
                lines.append("CellBender convergence and run warnings are retained when reported; exit code 0 records execution completion only and does not certify biological quality.")
            if task.capability == "donor_sensitivity" and len(task.output_artifacts) > 1:
                lines.append(f"LOO summary output: `{task.output_artifacts[1]}` (second registered output).")
            if task.capability == "liana_communication":
                lines.extend(["LIANA ranks are not FDR. Condition comparison FDR is BH-adjusted separately for each rank score. Non-significant and not-estimated comparisons are retained; inferred ligand–receptor communication does not establish signaling or causality.", ""])
            for uri in task.output_artifacts:
                meta, payload = self.artifact_registry.get(uri)
                lines.extend(["", f"Source: `{uri}`", ""])
                if not isinstance(payload, pd.DataFrame):
                    continue
                lines.append(f"Total rows: {len(payload)}.")
                if task.capability == "liana_communication":
                    append_liana_uncertainty(payload)
                if "fdr_q_value" in payload:
                    q = pd.to_numeric(payload["fdr_q_value"], errors="coerce")
                    lines.append(f"FDR < 0.05: {int(q.between(0, .05, inclusive='left').sum())}; unestimated FDR: {int(q.isna().sum())}.")
                if payload.empty:
                    lines.append("No eligible results were returned; no positive finding is inferred.")
                    continue
                cols = [c for c in columns if c in payload]
                if not cols:
                    cols = list(payload.columns[:8])
                lines.extend(["", "| " + " | ".join(cols) + " |", "| " + " | ".join(":---" for _ in cols) + " |"])
                for row in payload[cols].head(20).itertuples(index=False, name=None):
                    lines.append("| " + " | ".join(cell(v) for v in row) + " |")
            lines.append("")
        return lines

    @staticmethod
    def _task_target(task, slug_by_target):
        """Resolve target provenance attached by the orchestrator or encoded in a branch ID/URI."""

        if task is None or not slug_by_target:
            return None
        metrics = getattr(task, "metrics", {}) or {}
        target = metrics.get("target_cell_type")
        slug = metrics.get("target_branch")
        if target is not None and slug is not None:
            return str(target), str(slug)
        task_id = str(getattr(task, "task_id", ""))
        for candidate, candidate_slug in slug_by_target.items():
            if task_id.endswith(f"_{candidate_slug}"):
                return candidate, candidate_slug
            for uri in getattr(task, "output_artifacts", []) or []:
                if f"/{candidate_slug}/" in str(uri):
                    return candidate, candidate_slug
        return None
