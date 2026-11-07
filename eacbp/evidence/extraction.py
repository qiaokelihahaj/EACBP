"""Convert audited outputs to narrowly worded evidence, without disease templates."""
import math
import hashlib
import json
import pandas as pd
from eacbp.schemas.evidence import EvidenceNode, EvidenceRowSelector, EvidenceType, EvidenceStrength, EvidencePolarity


class EvidenceRowBindingError(ValueError):
    """Evidence does not identify and faithfully describe a unique source row."""


_ROW_IDENTITIES = ("contrast_label", "contrast_kind", "contrast_factor", "condition_a", "condition_b",
                   "target_cell_type", "target_branch")
_ROW_PARAMETERS = ("condition_a", "condition_b", "cond_ad", "cond_ctrl", "contrast_spec",
                   "target_cell_type", "target_branch")
_ROW_STATISTICS = {"log2_fold_change", "activity_effect_condition_a_vs_b", "lfc_se", "standard_error",
                   "effect_se", "ci_lower", "ci_upper", "ci_low", "ci_high", "ci95_lower", "ci95_upper",
                   "activity_se", "activity_ci_low", "activity_ci_high", "fdr_q_value",
                   "p_value", "status", "status_reason", "effect_definition", "alpha"}


def _row_json(value, path):
    if isinstance(value, dict):
        return {str(k): _row_json(v, f"{path}.{k}") for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_row_json(v, f"{path}[{i}]") for i, v in enumerate(value)]
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        raise EvidenceRowBindingError(f"non-finite value at {path}")
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise EvidenceRowBindingError(f"unsupported row binding value at {path}")


def _row_digest(selector, metrics, context):
    payload = {"key": selector.key, "artifact_parameters": selector.artifact_parameters,
               "metrics": {k: metrics[k] for k in selector.metric_columns},
               "context": {k: context[k] for k in selector.context_columns}}
    content = json.dumps(_row_json(payload, "row"), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def _same_row_value(left, right):
    # CSV decimal roundtrips can differ in the final binary digit. This is a
    # serialization tolerance, never a scientific effect-size tolerance.
    if isinstance(left, (int, float)) and not isinstance(left, bool) and isinstance(right, (int, float)) and not isinstance(right, bool):
        return math.isfinite(left) and math.isfinite(right) and math.isclose(left, right, rel_tol=1e-12, abs_tol=0.0)
    return type(left) is type(right) and left == right


def _source_target(values, parameters):
    """Resolve a population from source fields; reject conflicting identities."""
    nested = parameters.get("inference_contract") or {}
    if not isinstance(nested, dict):
        raise EvidenceRowBindingError("source inference contract must be an object")
    candidates = [values.get("target_cell_type"), parameters.get("target_cell_type"),
                  nested.get("target_cell_type")]
    candidates = [value for value in candidates if value is not None and value != ""]
    if candidates and any(not _same_row_value(candidates[0], value) for value in candidates[1:]):
        raise EvidenceRowBindingError("conflicting target_cell_type identities in source row/metadata")
    return candidates[0] if candidates else None


def validate_evidence_target_scope(node, source_parameters, task=None, source_values=None):
    """Check a claimed population against frozen source metadata and task identity.

    Legacy aggregate evidence without a row selector remains explicitly
    unavailable. An old untargeted row does not acquire an invented population.
    Existing v3 selectors can be verified using their frozen source metadata
    even when they did not inline the nested inference contract.
    """
    if node.row_selector is None:
        return "unavailable"
    if source_values is None:
        # Offline checks may use only an explicitly bound identity key, never
        # an arbitrary metric added to a node. Live checks pass the real row.
        selector = node.row_selector
        if ("target_cell_type" in selector.key
                and selector.metric_columns.get("target_cell_type") == "target_cell_type"):
            source_values = {"target_cell_type": selector.key["target_cell_type"]}
        else:
            source_values = {}
    source_target = _source_target(source_values, source_parameters)
    claimed = node.biological_context.get("target_cell_type")
    task_targets = []
    if task is not None:
        inference = getattr(task, "inference_contract", None)
        if inference is not None:
            target = inference.get("target_cell_type") if isinstance(inference, dict) else getattr(inference, "target_cell_type", None)
            if target is not None and target != "":
                task_targets.append(target)
        for collection in (getattr(task, "metrics", {}) or {}, getattr(task, "parameters", {}) or {}):
            target = collection.get("target_cell_type")
            if target is not None and target != "":
                task_targets.append(target)
    if source_target is None:
        if (claimed is not None and claimed != "") or task_targets:
            raise EvidenceRowBindingError("target_cell_type is claimed but absent from source row/metadata")
        return "not_claimed"
    if not _same_row_value(claimed, source_target):
        raise EvidenceRowBindingError("evidence target_cell_type differs from source row/metadata")
    if any(not _same_row_value(target, source_target) for target in task_targets):
        raise EvidenceRowBindingError("source target_cell_type differs from source task target")
    return "verified"


def _row_context(values, parameters, evidence_type):
    if evidence_type in {EvidenceType.PSEUDOBULK_DEG, EvidenceType.CELL_LEVEL_DEG}:
        expected = {"gene": str(values["gene"]), "direction": "higher" if float(values["log2_fold_change"]) > 0 else "lower"}
    elif evidence_type == EvidenceType.FUNCTIONAL_ACTIVITY:
        expected = {"pathway": str(values["source"])}
    else:
        raise EvidenceRowBindingError("row binding is only implemented for DEG and functional activity")
    contrast = values.get("contrast_label") or parameters.get("contrast_spec", {}).get("label")
    if not contrast:
        tested = values.get("condition_a", parameters.get("condition_a", parameters.get("cond_ad", "test")))
        reference = values.get("condition_b", parameters.get("condition_b", parameters.get("cond_ctrl", "reference")))
        contrast = f"{tested} versus {reference}"
    expected["contrast"] = str(contrast)
    for name in ("target_cell_type", "target_branch"):
        value = _source_target(values, parameters) if name == "target_cell_type" else values.get(name, parameters.get(name))
        if value is not None and value != "":
            expected[name] = value
    return expected


def _row_summary(evidence_type, values, parameters, task_metrics, contract_parameters=None):
    """Render native inferential text from the values that are actually cited."""
    context = _row_context(values, parameters, evidence_type)
    alpha = _inference_alpha(values, task_metrics, parameters, contract_parameters)
    if alpha is None:
        raise EvidenceRowBindingError("row evidence has no valid inference threshold")
    q = float(values["fdr_q_value"])
    if not math.isfinite(q) or not 0 <= q < alpha:
        raise EvidenceRowBindingError("native inferential row evidence must pass its recorded FDR threshold")
    if "status" in values and values["status"] != "estimated":
        raise EvidenceRowBindingError("native inferential row evidence requires estimated status")
    contrast = context["contrast"]
    if evidence_type == EvidenceType.FUNCTIONAL_ACTIVITY:
        effect = float(values["activity_effect_condition_a_vs_b"])
        definition = str(values.get("effect_definition") or "donor-level inferred activity difference")
        return (f"Inferred activity of {context['pathway']} differs for {contrast} ({definition}; "
                f"donor-level effect={effect:.3g}, FDR={q:.3g} at alpha={alpha:g}); activity is inferred "
                "from the supplied network, not directly measured.")
    fc = float(values["log2_fold_change"])
    unit = "donor-level" if evidence_type == EvidenceType.PSEUDOBULK_DEG else "exploratory cell-level"
    definition = str(values.get("effect_definition") or "log2 fold change")
    effect_text = (f"has a {context['direction']} {definition} estimate" if values.get("contrast_kind") == "numeric"
                   else f"has {context['direction']} expression for {contrast}")
    return f"{context['gene']} {effect_text} ({unit}; {definition}; log2FC={fc:.3g}, FDR={q:.3g} at alpha={alpha:g})."


def validate_evidence_row_structure(node):
    """Validate frozen bindings without claiming that payload files were read."""
    selector = node.row_selector
    if selector is None:
        return "unavailable"
    if selector.artifact_uri not in node.source_artifact_uris:
        raise EvidenceRowBindingError("row selector artifact is not an evidence source")
    required = ({"source", "activity_effect_condition_a_vs_b", "fdr_q_value"}
                if node.type == EvidenceType.FUNCTIONAL_ACTIVITY else {"gene", "log2_fold_change", "fdr_q_value"})
    if not required.issubset(selector.metric_columns):
        raise EvidenceRowBindingError("row selector omits required scientific fields")
    for name in _ROW_STATISTICS.intersection(node.metrics):
        if name not in selector.metric_columns:
            raise EvidenceRowBindingError(f"row selector omits scientific metric {name}")
    if any(name not in node.metrics for name in selector.metric_columns):
        raise EvidenceRowBindingError("row selector refers to an absent evidence metric")
    if any(name not in node.biological_context for name in selector.context_columns):
        raise EvidenceRowBindingError("row selector refers to absent biological context")
    feature = "source" if node.type == EvidenceType.FUNCTIONAL_ACTIVITY else "gene"
    if feature not in selector.key:
        raise EvidenceRowBindingError("row selector requires a gene or source identity")
    for column, expected in selector.key.items():
        value = _row_json(expected, f"key.{column}")
        if value is None or isinstance(value, (dict, list)):
            raise EvidenceRowBindingError("row selector key must contain non-null scalar values")
        mapped = [name for name, source in selector.metric_columns.items() if source == column]
        if not mapped or any(not _same_row_value(_row_json(node.metrics[name], name), value) for name in mapped):
            raise EvidenceRowBindingError(f"row key differs from evidence metric {column}")
    expected_context = _row_context(node.metrics, selector.artifact_parameters, node.type)
    for name, expected in expected_context.items():
        if not _same_row_value(node.biological_context.get(name), expected):
            raise EvidenceRowBindingError(f"row evidence context differs for {name}")
    if _row_digest(selector, node.metrics, node.biological_context) != selector.row_sha256:
        raise EvidenceRowBindingError("row binding fingerprint differs from evidence values")
    return "structure_checked"


def validate_evidence_row_binding(node, registry, task=None):
    """Independently resolve the unique row and compare every bound value."""
    status = validate_evidence_row_structure(node)
    if status == "unavailable":
        return status
    selector = node.row_selector
    metadata, table = registry.get(selector.artifact_uri)
    if not isinstance(table, pd.DataFrame) or not table.columns.is_unique:
        raise EvidenceRowBindingError("row source must have unique table columns")
    if task is not None and selector.artifact_uri not in task.output_artifacts:
        raise EvidenceRowBindingError("row artifact is not an output of its source task")
    for name, expected in selector.artifact_parameters.items():
        if name not in metadata.parameters or _row_json(metadata.parameters[name], name) != expected:
            raise EvidenceRowBindingError(f"row source artifact parameter differs for {name}")
    for name in _ROW_PARAMETERS:
        if name in metadata.parameters and metadata.parameters[name] is not None and name not in selector.artifact_parameters:
            raise EvidenceRowBindingError(f"row selector omits artifact parameter {name}")
    mask = pd.Series(True, index=table.index)
    for name, value in selector.key.items():
        if name not in table.columns:
            raise EvidenceRowBindingError(f"row source has no key column {name}")
        mask &= table[name].eq(value).fillna(False)
    selected = table.loc[mask]
    if len(selected) != 1:
        raise EvidenceRowBindingError(f"row selector must resolve uniquely; matched {len(selected)} rows")
    row = _evidence_values(selected.iloc[0].to_dict(), getattr(task, "metrics", {}) or {}, metadata.parameters)
    # Empty covariates may also be certified by artifact computation parameters.
    if "covariates" in row and pd.api.types.is_scalar(row["covariates"]) and pd.isna(row["covariates"]) and metadata.parameters.get("covariates") == []:
        row["covariates"] = []
    validate_evidence_target_scope(node, metadata.parameters, task, source_values=row)
    for name in _ROW_STATISTICS.intersection(row):
        if name not in selector.metric_columns:
            raise EvidenceRowBindingError(f"row selector omits source scientific metric {name}")
    for name in _ROW_IDENTITIES:
        if name in row and row[name] is not None and name not in selector.key:
            raise EvidenceRowBindingError(f"row selector omits source identity {name}")
    for name, column in selector.metric_columns.items():
        if column not in row:
            raise EvidenceRowBindingError(f"row source has no metric column {column}")
        expected = _row_json(row[column], f"metrics.{name}")
        actual = _row_json(node.metrics[name], f"metrics.{name}")
        if not _same_row_value(actual, expected):
            raise EvidenceRowBindingError(f"evidence metric differs from source row: {name}")
    for name, column in selector.context_columns.items():
        expected = _row_json(row[column], f"context.{name}") if column in row else None
        if name in {"gene", "pathway"} and expected is not None:
            expected = str(expected)
        if column not in row or not _same_row_value(node.biological_context[name], expected):
            raise EvidenceRowBindingError(f"evidence context differs from source row: {name}")
    for name, expected in _row_context(row, metadata.parameters, node.type).items():
        if not _same_row_value(node.biological_context.get(name), expected):
            raise EvidenceRowBindingError(f"evidence context differs from source row: {name}")
    if node.summary != _row_summary(node.type, row, metadata.parameters, getattr(task, "metrics", {}) or {}):
        raise EvidenceRowBindingError("evidence summary differs from the native source-row summary")
    return "verified"


def validate_evidence_row_bindings(nodes, registry, task_history=()):
    """Check each node while loading a shared table only once per boundary."""
    class TableCache:
        def __init__(self):
            self.values = {}

        def get(self, uri):
            if uri not in self.values:
                self.values[uri] = registry.get(uri)
            return self.values[uri]

    tasks = {task.task_id: task for task in task_history}
    cache = TableCache()
    return {node.evidence_id: validate_evidence_row_binding(node, cache, tasks.get(node.source_task_id))
            for node in nodes}


def _bind_row(node, metadata, row, registry, task):
    feature = "source" if node.type == EvidenceType.FUNCTIONAL_ACTIVITY else "gene"
    values = _evidence_values(row.to_dict(), task.metrics, metadata.parameters)
    # Helper-only selection columns are not source table columns.
    values = {name: value for name, value in values.items() if not name.startswith("__")}
    node.metrics.update(values)
    key = {name: _row_json(values[name], f"key.{name}") for name in (feature,) + _ROW_IDENTITIES
           if name in values and values[name] is not None}
    parameters = {name: _row_json(metadata.parameters[name], name) for name in _ROW_PARAMETERS
                  if name in metadata.parameters and metadata.parameters[name] is not None}
    # New evidence freezes nested population identity. Older selectors may
    # omit this copy; their immutable source metadata is still checked above.
    if isinstance(metadata.parameters.get("inference_contract"), dict):
        parameters["inference_contract"] = _row_json(metadata.parameters["inference_contract"], "inference_contract")
    node.biological_context.update(_row_context(values, parameters, node.type))
    contexts = ({"pathway": "source"} if feature == "source" else {"gene": "gene"})
    contexts.update({name: name for name in ("contrast", "target_cell_type", "target_branch")
                     if (name if name != "contrast" else "contrast_label") in values})
    if "contrast" in contexts:
        contexts["contrast"] = "contrast_label"
    selector = EvidenceRowSelector(artifact_uri=metadata.uri, key=key,
        metric_columns={name: name for name in values}, context_columns=contexts,
        artifact_parameters=parameters, row_sha256="0" * 64)
    selector.row_sha256 = _row_digest(selector, node.metrics, node.biological_context)
    node.row_selector = selector
    validate_evidence_row_binding(node, registry, task)


def finite_number(value, default=None):
    try:
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def _evidence_values(values, task_metrics, artifact_parameters=None):
    """Normalize certified empty descriptors from legacy CSV or typed tables."""
    values = dict(values or {})

    def missing_descriptor(value):
        return pd.api.types.is_scalar(value) and bool(pd.isna(value))

    covariates = task_metrics.get("covariates", (artifact_parameters or {}).get("covariates"))
    empty_covariates = isinstance(values.get("covariates"), str) and values["covariates"] == ""
    if empty_covariates and isinstance(covariates, (list, tuple)) and covariates:
        raise EvidenceRowBindingError("empty source descriptor for declared nonempty covariates")
    if ("covariates" in values and (missing_descriptor(values["covariates"]) or empty_covariates)
            and isinstance(covariates, (list, tuple)) and not covariates):
        values["covariates"] = []
    if (values.get("status") == "estimated" and "status_reason" in values
            and missing_descriptor(values["status_reason"])):
        values["status_reason"] = ""
    if values.get("contrast_kind") == "numeric":
        # Numeric contrasts explicitly have no categorical factor or levels;
        # CSV encodes these nullable descriptors as NaN, not missing estimates.
        for name in ("contrast_factor", "contrast_tested_level", "contrast_reference_level"):
            if name in values and missing_descriptor(values[name]):
                values[name] = None
    for name in ("network_path", "donor_role_limitation"):
        declared = task_metrics if name in task_metrics else (artifact_parameters or {})
        if (name in values and missing_descriptor(values[name]) and name in declared
                and (declared[name] is None or declared[name] == "")):
            values[name] = declared[name]
    # Statistical NaN/Inf values remain unchanged and fail strict snapshots.
    return values


def _inference_alpha(row=None, task_metrics=None, artifact_parameters=None, contract_parameters=None):
    """Read the threshold used by an audited result, with legacy defaults."""
    for source in (row, task_metrics, artifact_parameters, contract_parameters):
        if source is None:
            continue
        for key in ("alpha", "significance_alpha"):
            if key not in source:
                continue
            value = finite_number(source.get(key))
            if value is None or not 0.0 < value < 1.0:
                return None
            return value
    return 0.05


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
    inference_id = getattr(result, "inference_contract_id", None)
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
        source_row=None,
        source_metadata=None,
    ):
        values = _evidence_values(values, metrics)
        values["capability"] = cap
        node = EvidenceNode(
            evidence_id=f"E_{contract.task_id}_{suffix or len(evidence)}",
            type=kind, summary=summary, metrics=values, polarity=polarity,
            biological_context=context or {}, score=max(0, min(1, score)),
            strength=strength,
            source_task_id=contract.task_id, source_artifact_uris=uris,
            inference_contract_id=inference_id,
            audit_passed=True, is_simulated=simulated,
            data_origin_uris=sorted(roots),
        )
        if source_row is not None:
            _bind_row(node, source_metadata, source_row, registry, result)
        evidence.append(node)

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

    table_output = next(((m, p) for m, p in outputs if isinstance(p, pd.DataFrame)), (None, None))
    table_metadata, table = table_output
    scientific_result = getattr(result, "scientific_result", None)
    if inference_id and scientific_result is not None and cap in {"deg", "functional_activity"}:
        scientific_status = getattr(scientific_result.status, "value", scientific_result.status)
        if scientific_status in {"estimated_supported", "estimated_inconclusive", "not_estimable"}:
            add(
                EvidenceType.STATISTICAL_RESULT,
                scientific_result.summary,
                {"scientific_status": scientific_status,
                 "n_features_tested": scientific_result.n_features_tested,
                 "n_features_supported": scientific_result.n_features_supported,
                 "alpha": metrics.get("alpha", contract.parameters.get("alpha", 0.05))},
                score=0, suffix="result_summary", polarity=EvidencePolarity.NEUTRAL,
                strength=EvidenceStrength.INSUFFICIENT,
            )
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
        artifact_parameters = table_metadata.parameters
        for _, row in table.iterrows():
            q = finite_number(row.get("fdr_q_value"))
            effect = finite_number(row.get("activity_effect_condition_a_vs_b"))
            alpha = _inference_alpha(row, metrics, artifact_parameters, contract.parameters)
            if alpha is None or q is None or effect is None or not 0 <= q < alpha:
                continue
            source = str(row.get("source"))
            contrast = str(row.get("contrast_label") or f"{row.get('condition_a')} versus {row.get('condition_b')}")
            add(
                EvidenceType.FUNCTIONAL_ACTIVITY,
                _row_summary(EvidenceType.FUNCTIONAL_ACTIVITY, row, artifact_parameters, metrics, contract.parameters),
                row.to_dict(), score=1-q, context={"pathway": source, "contrast": contrast},
                source_row=row, source_metadata=table_metadata,
            )
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
        artifact_parameters = table_metadata.parameters
        selected = table.copy()
        selected["fdr_q_value"] = pd.to_numeric(selected["fdr_q_value"], errors="coerce")
        if "alpha" in selected:
            selected["__inference_alpha"] = pd.to_numeric(selected["alpha"], errors="coerce")
        else:
            selected["__inference_alpha"] = _inference_alpha(
                None, metrics, artifact_parameters, contract.parameters
            )
        selected = selected[
            selected["fdr_q_value"].notna()
            & selected["fdr_q_value"].ge(0.0)
            & selected["__inference_alpha"].between(0.0, 1.0, inclusive="neither")
            & selected["fdr_q_value"].lt(selected["__inference_alpha"])
        ].sort_values("fdr_q_value").head(10)
        for _, row in selected.iterrows():
            fc, q = finite_number(row["log2_fold_change"]), finite_number(row["fdr_q_value"])
            if fc is None or fc == 0 or q is None:
                continue
            alpha = finite_number(row.get("__inference_alpha"), 0.05)
            gene = str(row["gene"])
            direction = "higher" if fc > 0 else "lower"
            legacy_tested = artifact_parameters.get(
                "condition_a", artifact_parameters.get("cond_ad", "test")
            )
            legacy_reference = artifact_parameters.get(
                "condition_b", artifact_parameters.get("cond_ctrl", "reference")
            )
            contrast = str(
                row.get("contrast_label")
                or artifact_parameters.get("contrast_spec", {}).get("label")
                or f"{row.get('condition_a', legacy_tested)} versus "
                f"{row.get('condition_b', legacy_reference)}"
            )
            evidence_type = EvidenceType.PSEUDOBULK_DEG if is_pb else EvidenceType.CELL_LEVEL_DEG
            add(evidence_type,
                _row_summary(evidence_type, row, artifact_parameters, metrics, contract.parameters),
                {name: value for name, value in row.items() if not name.startswith("__")},
                score=1-q, context={"gene": gene, "contrast": contrast, "direction": direction},
                source_row=row, source_metadata=table_metadata)
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
