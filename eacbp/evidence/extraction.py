"""Convert audited outputs to narrowly worded evidence, without disease templates."""
import math
import pandas as pd
from eacbp.schemas.evidence import EvidenceNode, EvidenceType, EvidenceStrength, EvidencePolarity


def finite_number(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def extract_evidence(contract, result, report, registry):
    if not report.overall_passed or report.stop_rule_triggered:
        return []
    uris = result.output_artifacts
    if not uris or any(not registry.exists(uri) for uri in uris):
        return []
    # Reading also verifies artifact integrity.
    outputs = [registry.get(uri) for uri in uris]
    cap = contract.capability
    metrics = result.metrics
    evidence = []

    # Follow every parent, rather than just a shortest lineage path. Different
    # computational branches from the same matrix are not independent cohorts.
    roots, seen = set(), set()
    pending = list(uris)
    while pending:
        uri = pending.pop()
        if uri in seen:
            continue
        seen.add(uri)
        metadata = registry.get_metadata(uri)
        if metadata.parent_uris:
            pending.extend(metadata.parent_uris)
        else:
            roots.add(uri)

    # Keep the simulation marker attached to the observation even when a
    # caller supplies it only on an upstream root.  ArtifactRegistry normally
    # inherits this marker, but checking the roots makes the extraction
    # boundary explicit and robust for independently constructed registries.
    simulated = bool(metrics.get("is_simulated")) or any(
        m.summary_metrics.get("is_simulated", False) for m, _ in outputs
    ) or any(
        registry.get_metadata(uri).summary_metrics.get("is_simulated", False)
        for uri in roots
    )

    def add(
        kind,
        summary,
        values=None,
        score=0.5,
        context=None,
        suffix=None,
        polarity=EvidencePolarity.SUPPORTING,
        strength=EvidenceStrength.MODERATE,
    ):
        values = dict(values or {})
        values["capability"] = cap
        evidence.append(EvidenceNode(
            evidence_id=f"E_{contract.task_id}_{suffix or len(evidence)}",
            type=kind, summary=summary, metrics=values, polarity=polarity,
            biological_context=context or {}, score=max(0, min(1, score)),
            strength=strength,
            source_task_id=contract.task_id, source_artifact_uris=uris,
            audit_passed=True, is_simulated=simulated,
            data_origin_uris=sorted(roots),
        ))

    if cap == "knowledge_retrieval":
        for i, raw in enumerate(metrics.get("evidence_nodes", [])):
            node = EvidenceNode.model_validate(raw)
            node.evidence_id = f"E_{contract.task_id}_{i}"
            node.source_task_id = contract.task_id
            node.source_artifact_uris = uris
            node.audit_passed = True
            node.is_simulated = simulated
            node.source_verified = False
            node.data_origin_uris = sorted(roots)
            node.metrics.update({"source_mode": "local_curated_unverified", "capability": cap})
            node.summary = "Local curated reference (not independently verified): " + node.summary
            evidence.append(node)
        return evidence

    table = next((p for _, p in outputs if isinstance(p, pd.DataFrame)), None)
    if cap == "liana_communication":
        add(EvidenceType.CELL_COMMUNICATION,
            f"LIANA inferred {metrics.get('n_unique_interactions', 0)} ligand–receptor interactions across {len(metrics.get('donor_condition_groups', []))} donor-condition groups. Ranks are not FDR; inferred communication does not establish signaling or causality.",
            metrics, score=0)
        return evidence

    if cap == "doublet_detection":
        observed = []
        for label, key in (
            ("cells before", "n_cells_before"),
            ("cells after", "n_cells_after"),
            ("doublets marked", "n_doublets_marked"),
            ("doublets filtered", "n_doublets_filtered"),
            ("retention rate", "retention_rate"),
            ("threshold", "threshold"),
        ):
            if key in metrics:
                observed.append(f"{label}={metrics[key]}")
        detail = ", ".join(observed)
        summary = (
            "Audited doublet detection recorded "
            + (detail if detail else "a QC output with no additional count metrics reported")
            + ". Doublet scores and labels are observations for QC; they do not establish a biological mechanism or causality."
        )
        add(
            EvidenceType.QC_METRICS,
            summary,
            metrics,
            score=0,
            context={"qc_observation": "doublet_detection"},
            polarity=EvidencePolarity.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
        )
        return evidence

    if cap == "cell_annotation":
        observed = []
        for label, key in (
            ("cells before", "n_cells_before"),
            ("cells after", "n_cells_after"),
            ("annotations", "n_annotations"),
            ("unknown annotations", "n_unknown_annotations"),
            ("label conflicts", "n_label_conflicts"),
        ):
            if key in metrics:
                observed.append(f"{label}={metrics[key]}")
        detail = ", ".join(observed)
        summary = (
            "Audited CellTypist annotation recorded "
            + (detail if detail else "an annotation output with no additional count metrics reported")
            + ". Unknown and conflicting labels remain explicit uncertainty; annotation does not establish a biological mechanism or causality."
        )
        add(
            EvidenceType.CELL_ANNOTATION,
            summary,
            metrics,
            score=0,
            context={"annotation_observation": "cell_annotation"},
            polarity=EvidencePolarity.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
        )
        return evidence

    if cap == "background_removal":
        observed = []
        for label, key in (
            ("cells before", "n_cells_before"),
            ("cells after", "n_cells_after"),
            ("genes", "n_genes"),
            ("return code", "returncode"),
            ("quality status", "quality_status"),
            ("scientific quality certified", "scientific_quality_certified"),
        ):
            if key in metrics:
                observed.append(f"{label}={metrics[key]}")
        warning_values = metrics.get("quality_warnings")
        if isinstance(warning_values, (list, tuple)):
            warnings = [str(value) for value in warning_values if str(value).strip()]
            if warnings:
                observed.append("retained quality warnings=" + "; ".join(warnings[:8]))
            elif "quality_warnings" in metrics:
                observed.append("retained quality warnings=none recorded")
        elif "quality_warnings" in metrics:
            observed.append(f"retained quality warnings={warning_values}")
        detail = ", ".join(observed)
        summary = (
            "Audited CellBender background-removal output recorded "
            + (detail if detail else "a corrected-count output with no additional quality metrics reported")
            + "."
        )
        if metrics.get("returncode") == 0:
            summary += " Exit code 0 records execution completion only and does not certify biological quality."
        else:
            summary += " Convergence and run diagnostics remain execution evidence and do not certify biological quality."
        summary += " Corrected counts do not establish a biological mechanism or causality."
        add(
            EvidenceType.QC_METRICS,
            summary,
            metrics,
            score=0,
            context={"qc_observation": "background_removal"},
            polarity=EvidencePolarity.NEUTRAL,
            strength=EvidenceStrength.INSUFFICIENT,
        )
        return evidence

    if cap == "functional_activity" and table is not None:
        for _, row in table.iterrows():
            q = finite_number(row.get("fdr_q_value"))
            effect = finite_number(row.get("activity_effect_condition_a_vs_b"))
            if q is None or effect is None or not 0 <= q < .05:
                continue
            source = str(row.get("source"))
            contrast = f"{row.get('condition_a')} versus {row.get('condition_b')}"
            add(EvidenceType.FUNCTIONAL_ACTIVITY,
                f"Inferred activity of {source} differs for {contrast} (donor-level effect={effect:.3g}, FDR={q:.3g}); activity is inferred from the supplied network, not directly measured.",
                row.to_dict(), score=1-q, context={"pathway": source, "contrast": contrast})
        return evidence
    if cap == "donor_sensitivity":
        if metrics.get("skipped"):
            add(EvidenceType.SENSITIVITY_ANALYSIS,
                f"Donor sensitivity was not evaluated: {metrics.get('skip_reason', 'insufficient eligible donors')}. No robustness conclusion is supported.",
                metrics, score=0)
        else:
            add(EvidenceType.SENSITIVITY_ANALYSIS,
                f"Donor leave-one-out analysis completed {metrics.get('n_successful_fits', 0)} of {metrics.get('n_requested_fits', 0)} requested fits. Per-gene effect and direction sensitivity must be inspected; fit completion alone does not establish robustness.",
                metrics, score=0)
        return evidence
    if cap == "deg" and table is not None:
        required = {"gene", "log2_fold_change", "fdr_q_value"}
        if not required.issubset(table.columns):
            return []
        is_pb = bool(metrics.get("is_pseudobulk", False))
        comparisons = outputs[0][0].parameters
        contrast = f"{comparisons.get('condition_a', comparisons.get('cond_ad', 'test'))} versus {comparisons.get('condition_b', comparisons.get('cond_ctrl', 'reference'))}"
        selected = table.copy()
        selected["fdr_q_value"] = pd.to_numeric(selected["fdr_q_value"], errors="coerce")
        selected = selected[selected["fdr_q_value"].between(0, 0.05, inclusive="left")].sort_values("fdr_q_value").head(10)
        for _, row in selected.iterrows():
            fc, q = finite_number(row["log2_fold_change"]), finite_number(row["fdr_q_value"])
            if fc is None or fc == 0 or q is None:
                continue
            gene = str(row["gene"])
            direction = "higher" if fc > 0 else "lower"
            unit = "donor-level" if is_pb else "exploratory cell-level"
            add(EvidenceType.PSEUDOBULK_DEG if is_pb else EvidenceType.CELL_LEVEL_DEG,
                f"{gene} has {direction} expression in {contrast} ({unit}; log2FC={fc:.3g}, FDR={q:.3g}).",
                row.to_dict(), score=1-q, context={"gene": gene, "contrast": contrast, "direction": direction})
        return evidence

    if cap == "differential_abundance":
        for row in metrics.get("abundance_results", []):
            q = finite_number(row.get("fdr_q_value", row.get("fdr")))
            if q is None or not 0 <= q < .05:
                continue
            add(EvidenceType.DIFFERENTIAL_ABUNDANCE,
                f"State {row.get('state')} differs in donor abundance for {row.get('condition_a')} versus {row.get('condition_b')} (log2 ratio={row.get('log2_ratio_condition_a_vs_b', row.get('log2_ratio'))}, FDR={q:.3g}).",
                row, score=1-q, context={"state": row.get("state")})
        return evidence

    if cap == "spatial_deg" and table is not None:
        if not {"gene", "fdr_q_value", "morans_i"}.issubset(table.columns):
            return []
        for _, row in table.sort_values("fdr_q_value").head(10).iterrows():
            q, moran = finite_number(row["fdr_q_value"]), finite_number(row["morans_i"])
            if q is None or moran is None or not 0 <= q < .05:
                continue
            add(EvidenceType.SPATIAL_LOCALIZATION,
                f"{row['gene']} shows spatial autocorrelation (Moran's I={moran:.3g}, FDR={q:.3g}).",
                row.to_dict(), score=1-q, context={"gene": str(row["gene"]), "spatial": True})
        return evidence

    if cap in ("cell_cell_communication", "spatial_cci") and table is not None:
        for _, row in table.head(5).iterrows():
            interaction = finite_number(row.get("spatial_interaction_score"))
            if interaction is None:
                continue
            add(EvidenceType.SPATIAL_LOCALIZATION,
                f"Computed proximity-weighted interaction score for {row.get('ligand')}-{row.get('receptor')} "
                f"between {row.get('sender_cell_type')} and {row.get('receiver_cell_type')}: {interaction:.3g}.",
                row.to_dict(), context={"ligand": row.get("ligand"), "receptor": row.get("receptor")})
        return evidence

    if cap == "trajectory_inference":
        stability = finite_number(metrics.get("stability_score"), 0)
        if stability < .60 or any(not c.passed for c in report.checks if c.check_name == "trajectory_subsampling_stability"):
            return []
        add(EvidenceType.TRAJECTORY_STABILITY,
            f"Root-distance ordering has subsampling rank correlation {stability:.3g}; "
            f"FDR-selected dynamic genes: {', '.join(metrics.get('top_dynamic_genes', [])) or 'none'}.",
            metrics, score=stability)
        return evidence

    if "perturbation" in cap:
        if "genetic" in cap:
            target = metrics.get("target_gene", contract.parameters.get("target_gene"))
            reversion = finite_number(metrics.get("reversion_rate"))
            summary = f"In silico {metrics.get('perturbation_type', contract.parameters.get('perturbation_type', 'perturbation'))} of {target}"
            summary += f" produced state reversion score {reversion:.3g}." if reversion is not None else " produced a counterfactual expression artifact."
        else:
            summary = f"In silico compound scoring produced {len(table) if table is not None else 0} candidate results."
        add(EvidenceType.PERTURBATION, summary, metrics, score=.35)
        return evidence

    kinds = {
        "quantification": EvidenceType.DATASET_AUDIT,
        "dataset_audit": EvidenceType.DATASET_AUDIT,
        "qc": EvidenceType.QC_METRICS,
        "clustering": EvidenceType.CLUSTERING_STABILITY,
        "spatial_domain": EvidenceType.SPATIAL_LOCALIZATION,
        "spacell_microenvironment_analysis": EvidenceType.SPATIAL_LOCALIZATION,
        "gene_function_reasoning": EvidenceType.PATHWAY_ENRICHMENT,
        "chatcell_dialogue_prediction": EvidenceType.TRAJECTORY_STABILITY,
    }
    if cap in kinds:
        if cap == "quantification":
            summary = f"{'Simulation generated' if simulated else 'Quantification produced'} {metrics.get('n_cells_quantified', 0)} cells using {metrics.get('quant_engine', result.method_used)}."
        elif cap == "dataset_audit":
            summary = f"Dataset audit counted {metrics.get('n_cells', 0)} cells; minimum recorded donor replicates per condition: {metrics.get('min_replicates', 0)}."
        elif cap == "qc":
            summary = f"QC retained fraction {metrics.get('retention_rate', 'unreported')} of input cells."
        else:
            summary = f"{result.method_used} produced {len(uris)} audited output artifacts for {cap}; this observation does not establish a biological mechanism."
        add(kinds[cap], summary, metrics)
    return evidence
