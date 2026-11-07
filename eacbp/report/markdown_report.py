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
        run_summary: Optional[Dict[str, Any]] = None,
    ):
        self.manifest = manifest
        self.evidence_graph = evidence_graph
        self.artifact_registry = artifact_registry
        self.task_history = task_history
        self.audit_reports = audit_reports or []
        self.run_summary = run_summary or {}
        self.tracker = SentenceProvenanceTracker(evidence_graph, artifact_registry, task_history)

    def generate_markdown(self) -> str:
        from eacbp.evidence.extraction import validate_evidence_row_bindings
        validate_evidence_row_bindings(self.evidence_graph.evidence_nodes.values(),
                                       self.artifact_registry, self.task_history)
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

        lines.extend(self._autonomous_request_details())

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
        lines.extend(self._scientific_result_details())
        lines.extend(self._qc_disposition_details())
        lines.extend(self._qc_review_details())
        lines.extend(self._clustering_execution_details())
        lines.extend(self._cluster_marker_details())
        lines.extend(self._rna_velocity_details())
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
        lines.extend(self._fit_diagnostics_details())
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
            if card.get("inference_contract_ids"):
                lines.append("Inference contracts: " + ", ".join(f"`{value}`" for value in card["inference_contract_ids"]))
            for heading, chain in (
                ("Supporting Evidence & Sentence-Level Provenance Trace", card["evidence_chain"]),
                ("Contradicting Evidence & Sentence-Level Provenance Trace", card.get("contradiction_evidence_chain", [])),
            ):
                lines.append(f"#### {heading}:")
                if not chain:
                    lines.append("No admitted evidence is attached in this category; this does not establish that contrary observations are absent.")
                for ev in chain:
                    lines.extend(self._evidence_trace_lines(ev, manifest))
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

    def _evidence_trace_lines(self, ev, manifest):
        import json

        task_info = ev["source_task"]
        lines = [
            f"- **Evidence `{ev['evidence_id']}`** (`{ev['evidence_type']}`, `{ev['strength']}`, score: `{ev['score']:.2f}`): {ev['summary']}",
            f"  - *Generated by Task*: `{task_info['task_id']}` (`{task_info['capability']}` -> `{task_info['method']}`)",
            f"  - *Source row verification*: `{ev.get('row_verification', 'unavailable')}`",
        ]
        selector = ev.get("row_selector")
        if selector:
            lines.append(f"  - *Unique source row*: `{selector['artifact_uri']}`; key `{json.dumps(selector['key'], sort_keys=True, ensure_ascii=False)}`")
            lines.append(f"  - *Bound row metrics*: `{json.dumps(ev.get('metrics', {}), sort_keys=True, ensure_ascii=False)}`")
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
        return lines

    def _fit_diagnostics_details(self):
        batches = []
        for task in self.task_history:
            if getattr(task, "capability", None) not in {"deg", "donor_sensitivity"}:
                continue
            reports = [report for report in self.audit_reports if report.target_task_id == task.task_id]
            status = getattr(task, "status", None)
            if getattr(status, "value", status) != "success" or not reports or any(not report.overall_passed or report.stop_rule_triggered for report in reports):
                continue
            for uri in getattr(task, "output_artifacts", []):
                meta = self.artifact_registry.get_metadata(uri)
                if meta.operation != "pydeseq2_fit_diagnostics":
                    continue
                _, payload = self.artifact_registry.get(uri)
                batches.append((uri, meta, payload))
        if not batches:
            return []
        lines = ["### Per-fit model diagnostics", "",
                 "Fit completion records execution. Warnings and dispersion-trend fallback do not by themselves prove failure or convergence. Worker-process stderr is outside the recorded Python-warning scope.",
                 "", "| Task / left-out donor | Fit ID | Fit status | Actual dispersion trend | Samples / genes | Recorded warnings | Receipt artifact |",
                 "| :--- | :--- | :--- | :--- | :--- | ---: | :--- |"]
        warning_lines = []
        for uri, meta, payload in batches:
            if not payload["fits"]:
                lines.append(f"| {meta.created_by_task} | unavailable | not_attempted | unavailable | unavailable | 0 | `{uri}` |")
                warning_lines.append(f"- Task `{meta.created_by_task}` not attempted: {payload.get('not_attempted_reason') or 'Reason unavailable'}")
            for fit in payload["fits"]:
                frame = fit["counts"]["frame"]
                lines.append(f"| {fit['task_id']} / {fit['left_out_donor'] or 'full model'} | `{fit['fit_id']}` | {fit['status']} | {fit['actual_dispersion_trend'] or 'unknown'} | {len(frame['index'])} / {len(frame['columns'])} | {len(fit['warnings'])} | `{uri}` |")
                if fit.get("error"):
                    warning_lines.append(f"- Fit `{fit['fit_id']}` failed: `{fit['error']['type']}`: {fit['error']['message']}")
                if fit["warnings"]:
                    for warning in fit["warnings"]:
                        message = warning["message"].replace("\n", " ").replace("|", "\\|")
                        warning_lines.append(f"- Fit `{fit['fit_id']}` ({fit['task_id']}, left-out donor `{fit['left_out_donor']}`): `{warning['category']}`: {message}")
        if warning_lines:
            lines.extend(["", *warning_lines])
        lines.extend(["", "The registered JSON retains each fit's counts, sample identities, design matrix, contrast, effective options, library versions, error and individual receipt digest.", ""])
        return lines

    def _autonomous_request_details(self):
        """Expose the persisted question contract independently of task success."""
        summary = self.run_summary
        if "analysis_request_status" not in summary:
            return []
        lines = ["### Scientific question completion", "",
                 f"- Pipeline execution: `{summary.get('status', 'unknown')}`",
                 f"- Requested analysis: `{summary['analysis_request_status']}`",
                 f"- Autonomy: `{summary.get('autonomy_kind', 'unknown')}`; policy `{summary.get('policy_version', 'unknown')}`",
                 "- Completion records required audited analyses; it does not certify agreement with an external reproduction benchmark."]
        checks = summary.get("requirement_completion_checks", [])
        if checks:
            lines.extend(["", f"Frozen requirement contract: `{summary.get('requirement_contract_id', 'unknown')}`", "",
                          "| Requested dimension | Population | Contrast | Status | Admitted tasks |",
                          "| :--- | :--- | :--- | :--- | :--- |"])
            def cell(value):
                return str(value).replace("|", "\\|").replace("\n", " ").replace("\r", " ")
            for check in checks:
                contrast = check.get("contrast")
                label = (f"{contrast.get('tested')} versus {contrast.get('reference')}"
                         if contrast and contrast.get("resolved") else "unresolved" if contrast else "—")
                row = (check.get("dimension"), check.get("target_cell_type") or "all",
                       label, check.get("status"), ", ".join(check.get("task_ids", [])) or "—")
                lines.append("| " + " | ".join(cell(value) for value in row) + " |")
            lines.append("")
        for gap in summary.get("scientific_gaps", []):
            lines.append(f"- Unanswered requirement: {gap}")
        for limitation in summary.get("analysis_limitations", []):
            lines.append(f"- Limitation: {limitation}")
        return lines + [""]

    def _qc_disposition_details(self):
        """Explain audited filtering reasons without equating low RNA with damage."""
        admitted = {report.target_task_id for report in self.audit_reports
                    if report.overall_passed and not report.stop_rule_triggered}
        tasks = [task for task in self.task_history
                 if task.capability == "qc" and task.task_id in admitted
                 and getattr(task.status, "value", task.status) == "success"
                 and task.metrics.get("cell_disposition_schema_version") == 1]
        if not tasks:
            return []
        lines = ["### QC cell dispositions", "",
                 "The registered QC artifact records every input cell, its measured QC quantities "
                 "and filtering flags. Reasons explain the executed thresholds; they do not "
                 "independently diagnose damaged cells.", ""]
        for task in tasks:
            lines.extend([f"Task: `{task.task_id}`; registered QC output: "
                          + ", ".join(f"`{uri}`" for uri in task.output_artifacts), "",
                          "| Disposition | Input cells |", "| :--- | ---: |"])
            for reason, count in sorted(task.metrics.get("cell_disposition_counts", {}).items()):
                lines.append(f"| {reason} | {count} |")
            warning = task.metrics.get("interpretation_warning")
            if warning:
                lines.extend(["", str(warning)])
            lines.append("")
        return lines

    def _qc_review_details(self):
        """Expose candidate counts only with a matching durable admission."""
        tasks = [task for task in self.task_history if task.capability == "qc_review"]
        if not tasks:
            return []
        records = self.artifact_registry.list_audit_records()
        lines = ["### Excluded-cell QC review", "",
                 "Excluded cells remain separate from the primary analysis. Raw counts and filtering "
                 "reasons are preserved for inspection; low RNA does not independently diagnose damage.", ""]
        for task in tasks:
            identity = ArtifactRegistry._result_identity(task.model_dump(mode="python"))
            matching = [record for record in records if record.task_id == task.task_id
                        and record.result_fingerprint == identity and record.overall_passed
                        and not record.stop_rule_triggered and record.status.value == "passed"
                        and record.output_artifacts == task.output_artifacts]
            if getattr(task.status, "value", task.status) != "success" or len(matching) != 1:
                lines.extend([f"Task `{task.task_id}`: review unavailable; no unique matching durable audit admission.", ""])
                continue
            record = matching[0]
            outputs = self.artifact_registry.get_audited_many(record.signature, record.contract)
            metadata, summary = next((meta, payload) for meta, payload in outputs
                                     if meta.operation == "raw_count_qc_review_summary")
            lines.extend([f"Task: `{task.task_id}`; prepared-input SHA-256: `{summary['source_sha256']}`.",
                          f"Primary retained cells: {summary['retained_cells']}; excluded candidates: {summary['candidate_cells']}; "
                          f"preserved input features: {summary['n_genes']}.", "",
                          "| Disposition | Input cells |", "| :--- | ---: |"])
            for reason, count in summary["cell_disposition_counts"].items():
                lines.append(f"| {reason} | {count} |")
            lines.extend(["", f"Mitochondrial metric source: `{summary['mitochondrial_metric_source']}`.",
                          str(summary["interpretation_warning"]), ""])
            for meta, _ in outputs:
                path = str(meta.storage_path).replace("\\", "/")
                lines.append(f"- `{meta.uri}`: [registered file](<{path}>); SHA-256 `{meta.sha256_hash}`.")
            features = next(payload for meta, payload in outputs if meta.operation == "qc_review_feature_counts")
            lines.extend(["", "Raw detection descriptions below use the excluded cells themselves. "
                          "They are not cell-type assignments, differential-expression tests or evidence of viability.", "",
                          "| Filtering reason | Gene | Detected cells | Cells in group | Detection fraction | Raw count sum |",
                          "| :--- | :--- | ---: | ---: | ---: | ---: |"])
            groups = features[features["reason"].ne("retained") & features["n_cells"].gt(0)]
            for reason, rows in groups.groupby("reason", sort=True, observed=True):
                ranked = rows.sort_values(["detection_fraction", "total_counts", "gene_index"], ascending=[False, False, True]).head(5)
                for _, row in ranked.iterrows():
                    gene = str(row["gene"]).replace("|", "\\|").replace("\n", " ")
                    lines.append(f"| {reason} | {gene} | {int(row['n_detected'])} | {int(row['n_cells'])} | "
                                 f"{row['detection_fraction']:.4g} | {int(row['total_counts'])} |")
            lines.append("")
        return lines

    def _cluster_marker_details(self):
        """Show reference-free feature rankings only from admitted evidence."""
        from eacbp.schemas.evidence import EvidenceType

        nodes = [node for node in self.evidence_graph.evidence_nodes.values()
                 if node.type == EvidenceType.CLUSTER_MARKERS and node.audit_passed]
        if not nodes:
            return []
        lines = ["### Descriptive cluster markers", "",
                 "These genes characterize groups found in this dataset. Effects are differences "
                 "of mean log1p library-normalized expression. The same cells define the clusters "
                 "and their rankings; these are exploratory descriptions, not donor-level tests "
                 "or independently established cell-type identities.", ""]
        for node in nodes:
            lines.extend([f"Evidence: `{node.evidence_id}`; registered table: "
                          + ", ".join(f"`{uri}`" for uri in node.source_artifact_uris), "",
                          "| Cluster | Gene | Rank | Mean log1p expression difference | Detection fraction difference |",
                          "| :--- | :--- | ---: | ---: | ---: |"])
            for cluster, markers in sorted(node.metrics.get("top_markers_by_cluster", {}).items()):
                for row in markers[:5]:
                    label = str(row["gene"]).replace("|", "\\|").replace("\n", " ")
                    lines.append(f"| {cluster} | {label} | {row['rank']} | "
                                 f"{row['mean_log1p_difference']:.4g} | {row['detection_difference']:.4g} |")
            lines.append("")
        return lines

    def _clustering_execution_details(self):
        """Present effective settings from admitted payloads and their contracts."""
        tasks = [task for task in self.task_history
                 if task.capability == "clustering" and task.method_used == "scanpy_leiden_umap_v1"]
        if not tasks:
            return []
        from eacbp.capabilities.sc_data import SCData
        from eacbp.auditor.clustering import validate_clustering_execution_receipt
        records = self.artifact_registry.list_audit_records()
        lines = ["### Actual clustering execution", "",
                 "Execution reconstruction checks the recorded representation, graph and partition. "
                 "It does not establish cell-type truth or biological independence. UMAP is a display; "
                 "two Leiden iterations do not certify optimization convergence.", ""]
        for task in tasks:
            identity = ArtifactRegistry._result_identity(task.model_dump(mode="python"))
            matching = [record for record in records if record.task_id == task.task_id
                        and record.result_fingerprint == identity and record.overall_passed
                        and not record.stop_rule_triggered and record.status.value == "passed"
                        and record.output_artifacts == task.output_artifacts]
            if task.status.value != "success" or len(matching) != 1:
                lines.extend([f"Task `{task.task_id}`: effective settings unavailable; no unique matching durable admission.", ""])
                continue
            record = matching[0]
            outputs = self.artifact_registry.get_audited_many(record.signature, record.contract)
            data = SCData.from_payload(outputs[0][1])
            if "execution_receipt" not in data.uns.get("clustering", {}):
                lines.extend([f"Task `{task.task_id}`: effective settings unavailable; historical computation has no execution receipt.", ""])
                continue
            receipt = validate_clustering_execution_receipt(record.contract, task, self.artifact_registry,
                                                            recompute_partition=False)
            representation, neighbors, leiden = (receipt[key] for key in ("representation", "neighbors", "leiden"))
            lines.extend([f"Task: `{task.task_id}`; input SHA-256: `{receipt['source_sha256']}`.", "",
                          "| Recorded quantity | Effective value |", "| :--- | :--- |",
                          f"| Representation algorithm | `{representation['algorithm']}` |",
                          f"| Representation shape / features | `{representation['shape']}` / `{representation['feature_selection']}` |",
                          f"| Actual representation dimensions | {representation['actual_n_components']} |",
                          f"| Requested / actual neighbors | {neighbors['requested_n_neighbors']} / {neighbors['actual_n_neighbors']} |",
                          f"| Graph representation / metric | `{neighbors['use_rep']}` / `{neighbors['metric']}` |",
                          f"| Leiden resolution / seed | {leiden['resolution']} / {leiden['random_state']} |",
                          f"| Leiden flavor / iterations | `{leiden['flavor']}` / {leiden['n_iterations']} |",
                          f"| Clusters / partition SHA-256 | {leiden['n_clusters']} / `{leiden['partition_sha256']}` |", "",
                          "Representation settings are inherited when an upstream embedding exists; "
                          "a requested clustering PCA dimension is not an additional representation fit. "
                          "This presentation validates stored receipt consistency without refitting the graph.", ""])
        return lines

    def _rna_velocity_details(self):
        from eacbp.schemas.evidence import EvidenceType

        nodes = [node for node in self.evidence_graph.evidence_nodes.values()
                 if node.type == EvidenceType.RNA_VELOCITY and node.audit_passed]
        if not nodes:
            return []
        lines = ["### Model-inferred RNA velocity", "",
                 "Stochastic RNA velocity estimates kinetic direction from spliced and unspliced RNA. "
                 "Pseudotime has no measured time units. The audit checks preserved counts, neighborhoods, "
                 "RNA normalization, moments and summary consistency; it does not verify biological kinetic "
                 "assumptions or independently refit the model. These observations do not establish causal fate.", ""]
        for node in nodes:
            lines.extend([f"Evidence: `{node.evidence_id}`; outputs: "
                          + ", ".join(f"`{uri}`" for uri in node.source_artifact_uris), "",
                          "| Quantity | Recorded value |", "| :--- | :--- |"])
            for key in ("n_cells", "n_genes", "n_velocity_genes", "direction_source", "root_source",
                        "mean_velocity_confidence", "pseudotime_min", "pseudotime_max"):
                lines.append(f"| `{key}` | {node.metrics.get(key, 'not recorded')} |")
            lines.append("")
        return lines

    def _scientific_result_details(self):
        """Keep valid inconclusive estimates and failed preconditions visible."""
        selected = [task for task in self.task_history if getattr(task, "scientific_result", None) is not None]
        if not selected:
            return []
        audited = {report.target_task_id for report in self.audit_reports
                   if report.overall_passed and not report.stop_rule_triggered}
        lines = ["### Scientific result status and inference contracts", "",
                 "Non-significance does not establish absence of an effect or practical equivalence. "
                 "An unknown assumption is not a verified assumption. Failed or unaudited results below "
                 "describe execution limits and are not admitted scientific evidence.", ""]

        def cell(value):
            return str(getattr(value, "value", value)).replace("|", "\\|").replace("\n", " ")

        for task in selected:
            outcome = task.scientific_result
            contract = getattr(task, "inference_contract", None)
            admitted = task.task_id in audited and getattr(task.status, "value", task.status) == "success"
            lines.extend([f"#### `{task.task_id}` — {cell(outcome.status)}", "",
                          f"- Audit admission: {'passed' if admitted else 'not admitted'}",
                          f"- Inference contract: `{getattr(task, 'inference_contract_id', None) or 'not recorded'}`",
                          f"- Estimand: {cell(getattr(contract, 'estimand', 'not resolved'))}",
                          f"- Features tested: {outcome.n_features_tested}; supported: {outcome.n_features_supported}",
                          f"- Result: {cell(outcome.summary)}", ""])
            if contract is not None:
                lines.extend([f"- Contract version: {contract.contract_version}",
                              f"- Units: independent `{cell(contract.independent_unit)}`; observed `{cell(contract.observation_unit)}`",
                              f"- Design: `{cell(contract.design_formula)}`; counts: `{cell(contract.counts_source)}`",
                              f"- Alpha: {contract.alpha:g}; confidence level: {contract.confidence_level:g}; "
                              f"FDR: {cell(contract.fdr_method)} / {cell(contract.fdr_family)}", ""])
            if outcome.assumptions:
                lines.extend(["| Assumption | Status | Scope | Evidence / limitation |", "| :--- | :--- | :--- | :--- |"])
                for assessment in outcome.assumptions:
                    lines.append(f"| {cell(assessment.name)} | {cell(assessment.status)} | {cell(assessment.scope)} | "
                                 f"{cell(assessment.evidence)}; {cell(assessment.reason)} |")
                lines.append("")
        return lines

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
                   "effect_definition", "alpha", "significant_at_alpha",
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
                   "baseline_significant_at_alpha", "baseline_significant_fdr05", "n_significant_loo", "median_log2_fold_change",
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
            for key in ("status", "reason", "skip_reason", "statistical_unit", "design_formula", "design", "contrast_spec", "inference_settings", "alpha", "effect_definition", "interpretation_limit", "limitations"):
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
                    alpha = float(task.metrics.get("alpha", 0.05))
                    lines.append(f"FDR < {alpha:g}: {int(q.between(0, alpha, inclusive='left').sum())}; unestimated FDR: {int(q.isna().sum())}.")
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
