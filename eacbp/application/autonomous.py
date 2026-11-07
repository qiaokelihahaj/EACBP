"""Rule-driven analysis from a scientific question and observed input.

This boundary accepts no workflow configuration or reference answer. Decisions
are recorded before the shared study lifecycle imports and prepares the input.
It is a deterministic scientific policy, not an LLM reasoning service.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import math
from pathlib import Path
import re
import time
from typing import Any, Mapping

import numpy as np
import pandas as pd

from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.biological_units import DONOR_ALIASES, resolve_donor_role
from eacbp.application.analysis_requirements import freeze_question_requirements
from eacbp.schemas.study import (
    BiologicalDesign, Constraints, ExperimentalDesign, StudyManifest,
)


POLICY_VERSION = "observed_single_cell_v3"
_ALIASES = {
    "condition": ("condition", "condition_id", "treatment", "stim", "label"),
    "donor": DONOR_ALIASES,
    "sample": ("sample_id", "sample", "specimen_id", "specimen"),
    "batch": ("batch", "batch_id", "library_id", "library"),
    "library": ("library_id", "library"),
    "section": ("section_id", "section", "slide_id", "slide", "slice_id", "slice"),
    "measured_time": ("time", "timepoint", "time_point", "collection_time", "developmental_time", "age", "embryonic_day"),
    "cell_type": ("cell_type", "celltype", "cell_type_label"),
}
_CONTROL_LEVELS = {
    "ctrl", "control", "con", "untreated", "unstim", "unstimulated", "baseline",
    "vehicle", "mock", "wt", "false", "0", "normal", "healthy",
}
_DERIVED_UNS = {"iroot", "paga", "neighbors", "pca", "umap", "rank_genes_groups", "log1p", "clustering", "integration", "normalization", "qc"}
_KINETIC_LAYERS = ("spliced", "unspliced")


def _derived_observation_field(name: str) -> bool:
    """Recognize common analysis exports independently of dataset or software.

    Compact names also cover camel-case and punctuation variants emitted by
    different single-cell formats. Experimental design columns and feature
    identifiers are preserved; quality metrics are recomputed from counts.
    """
    compact = re.sub(r"[^a-z0-9]", "", str(name).casefold())
    if any(token in compact for token in ("cluster", "annotation", "pseudotime", "leiden", "louvain", "velocity", "snnres", "spatialdomain")):
        return True
    if re.fullmatch(r"(?:x)?(?:tsne|umap|pca|diffmap|embedding)(?:[0-9]+|[xyz])?", compact) or re.fullmatch(r"pc[0-9]+", compact):
        return True
    if compact.startswith(("ncount", "nfeature", "totalcounts", "ngenesbycounts", "pctcounts", "percentmito", "sscore", "g2mscore", "cellcyclescore")):
        return True
    if compact.startswith(("rootcell", "rootscore", "rootprob", "endcell", "endpoint", "endscore", "endprob", "terminalstate", "terminalprob")):
        return True
    return compact in {
        "qcpassed", "phase", "latenttime", "initialsizespliced", "initialsizeunspliced",
        "iroot", "rootcell", "rootcells", "rootscore", "rootscores", "rootprobability",
        "endpoint", "endpoints", "endcell", "endcells", "endscore", "endscores",
        "terminalstate", "terminalstates", "terminalprobability",
    }


def _derived_unstructured_field(name: str) -> bool:
    normalized = str(name).casefold()
    compact = re.sub(r"[^a-z0-9]", "", normalized)
    return (normalized in _DERIVED_UNS or normalized.endswith("_colors")
            or any(token in compact for token in ("cluster", "annotation", "pseudotime", "leiden", "louvain", "velocity", "snnres"))
            or compact.startswith(("spatialconnectivities", "spatialdistances", "spatialneighbors", "spatialdomain"))
            or compact.startswith(("tsne", "umap", "pca", "neighbors", "embedding", "diffmap", "highlyvariable")))


def _complete_labels(values: pd.Series) -> bool:
    normalized = values.astype(str).str.strip()
    return bool(not values.isna().any() and normalized.ne("").all()
                and not normalized.str.casefold().isin({"nan", "none", "null", "<na>"}).any())


def _resolve_column(obs: pd.DataFrame, role: str) -> str | None:
    if role == "donor":
        return resolve_donor_role(obs, {}).column
    candidates = [name for name in _ALIASES[role] if name in obs]
    # A generic label column describes a condition only when its values carry
    # recognizable control semantics. It may otherwise be a reference label.
    if role == "condition" and "label" in candidates:
        values = {str(value).strip().casefold() for value in obs["label"].dropna().unique()}
        if len(values) != 2 or not values.intersection(_CONTROL_LEVELS):
            candidates.remove("label")
    if not candidates:
        return None
    if role == "batch" and any(name in {"batch", "batch_id"} for name in candidates):
        candidates = [name for name in candidates if name in {"batch", "batch_id"}]
    selected = candidates[0]
    for other in candidates[1:]:
        left = obs[selected].astype("string").str.strip()
        right = obs[other].astype("string").str.strip()
        if not left.equals(right):
            raise ValueError(f"Ambiguous {role} metadata: {selected!r} and {other!r} disagree")
    return selected


def _matrix_values(matrix) -> np.ndarray:
    return np.asarray(matrix.data if hasattr(matrix, "tocsr") else matrix)


def _integer_counts(matrix) -> bool:
    values = _matrix_values(matrix)
    if values.dtype.kind not in "biuf" or not np.isfinite(values).all() or np.any(values < 0):
        return False
    if values.size and np.longdouble(values.max()) >= np.longdouble(2) ** 63:
        return False
    return bool(values.dtype.kind in "biu" or np.equal(values, np.floor(values)).all())


def _counts_source(payload: SCData) -> str:
    if "counts" in payload.layers:
        if not _integer_counts(payload.layers["counts"]):
            raise ValueError("Named counts layer is invalid: finite nonnegative integer counts are required")
        return "layers.counts"
    if _integer_counts(payload.X):
        return "X"
    raw = payload.raw
    if raw is not None and _integer_counts(raw["X"]):
        return "raw.X"
    raise ValueError("Unsupported expression representation: no finite nonnegative integer count matrix was observed; normalized or debatched X cannot be declared raw counts")


def _genes_source(var: pd.DataFrame) -> str:
    for name in ("gene_name", "gene_symbol", "gene_symbols", "symbol", "feature_name"):
        if name in var:
            if not _complete_labels(var[name]):
                raise ValueError(f"Gene symbols in var[{name!r}] are incomplete")
            if var[name].astype(str).str.strip().duplicated().any():
                raise ValueError(f"Gene symbols in var[{name!r}] are duplicated; an explicit gene-identity resolution is required")
            return f"var.{name}"
    if not var.index.is_unique or any(not str(value).strip() for value in var.index):
        raise ValueError("Unique nonempty feature IDs or gene symbols are required")
    return "var.index"


def _uses_observed_groups(question: str) -> bool:
    text = question.casefold()
    context = bool(re.search(r"\b(?:existing|observed|provided|annotated)\b|已有|既有|已注释|观测|原有", text))
    groups = any(word in text for word in ("cell type", "cell-type", "population", "细胞群", "细胞类型", "细胞亚群"))
    return context and groups


def _contrast(levels: list[str]) -> dict[str, Any] | None:
    if len(levels) != 2:
        return None
    controls = [value for value in levels if value.strip().casefold() in _CONTROL_LEVELS]
    reference = controls[0] if len(controls) == 1 else levels[0]
    tested = next(value for value in levels if value != reference)
    return {"tested": tested, "reference": reference,
            "reason": "One observed level has control semantics" if len(controls) == 1
            else "No unique control semantics; deterministic lexical orientation, without a causal interpretation"}


def _target_eligibility(obs: pd.DataFrame, columns, contrast) -> list[dict[str, Any]]:
    cell_col, donor_col, condition_col = (columns.get(key) for key in ("cell_type", "donor", "condition"))
    if not contrast or not cell_col or not donor_col or not condition_col:
        return []
    result = []
    for target in sorted(map(str, obs[cell_col].dropna().unique())):
        subset = obs.loc[obs[cell_col].astype(str) == target]
        complete = all(_complete_labels(subset[name]) for name in (donor_col, condition_col))
        groups = {condition: sorted(set(subset.loc[subset[condition_col].astype(str) == condition, donor_col].dropna().astype(str)))
                  for condition in (contrast["tested"], contrast["reference"])}
        donors_a, donors_b = map(set, groups.values())
        overlap = donors_a & donors_b
        paired = bool(overlap)
        eligible = complete and min(len(donors_a), len(donors_b)) >= 2 and (not paired or donors_a == donors_b)
        reason = "Complete paired donor sets" if eligible and paired else "Disjoint independent donor sets" if eligible else "Incomplete metadata, fewer than two biological donors per condition, or unmatched paired donors"
        result.append({"cell_type": target, "n_cells": len(subset), "donors_by_condition": groups,
                       "paired": paired, "eligible": eligible, "reason": reason})
    return result


def _velocity_layers(payload: SCData, counts_source: str) -> dict[str, Any]:
    evidence = {}
    reasons = []
    for name in _KINETIC_LAYERS:
        value = payload.layers.get(name)
        valid = value is not None and tuple(value.shape) == payload.shape and _integer_counts(value)
        evidence[name] = {"present": value is not None, "shape": list(value.shape) if value is not None else None,
                          "validated_integer_counts": valid}
        if not valid:
            reasons.append(f"layers.{name} is absent, misaligned or not finite nonnegative integer counts")
    if counts_source == "raw.X":
        reasons.append("The selected raw feature space is not the input space aligned to kinetic layers")
    return {"available": not reasons, "layers": evidence,
            "reason": "; ".join(reasons) if reasons else "Aligned finite nonnegative integer spliced and unspliced matrices are observed"}


def _observed_spatial_coordinates(payload: SCData) -> tuple[dict[str, Any], np.ndarray | None]:
    """Identify measured coordinates independently of derived input embeddings.

    An observed coordinate system supplies locations, not previously inferred
    tissue domains or neighborhood graphs. Multiple sections remain separate
    observed systems; this receipt does not align or connect those systems.
    """
    candidates = [(f"obsm.{name}", payload.obsm[name]) for name in ("spatial", "spatial_coords") if name in payload.obsm]
    for cols in (("x_coord", "y_coord"), ("x", "y")):
        if all(name in payload.obs for name in cols):
            z = "z_coord" if cols[0] == "x_coord" else "z"
            names = [*cols, *([z] if z in payload.obs else [])]
            candidates.append(("obs." + ",".join(names), payload.obs[names].to_numpy()))
    if not candidates:
        return {"available": False, "reason": "Measured spatial coordinates were not observed"}, None
    source, values = candidates[0]
    try:
        coordinates = np.asarray(values, dtype=np.float64)
        if coordinates.ndim != 2 or coordinates.shape[0] != payload.n_obs or coordinates.shape[1] not in (2, 3):
            raise ValueError("coordinates must have aligned cells and two or three dimensions")
        if not np.isfinite(coordinates).all() or not np.any(np.var(coordinates, axis=0) > 0):
            raise ValueError("coordinates must be finite and have nonzero spatial variation")
        for alternate, values in candidates[1:]:
            other = np.asarray(values, dtype=np.float64)
            if other.shape != coordinates.shape or not np.array_equal(other, coordinates):
                raise ValueError(f"measured coordinate sources {source!r} and {alternate!r} disagree")
    except (TypeError, ValueError) as exc:
        return {"available": False, "source": source, "reason": str(exc)}, None
    identity = hashlib.sha256(np.ascontiguousarray(coordinates).tobytes()).hexdigest()
    return {"available": True, "source": source, "shape": list(coordinates.shape),
            "sha256": identity, "reason": "Aligned finite measured coordinates are observed; no spatial graph or section alignment is reused"}, coordinates


def observe_input(payload: SCData, *, question: str) -> dict[str, Any]:
    """Observe matrix identity and sample design without modifying the input."""
    if min(payload.shape) < 3:
        raise ValueError("Autonomous single-cell analysis requires at least three cells and features")
    if not payload.obs.index.is_unique:
        raise ValueError("Input cell IDs must be unique")
    source = _counts_source(payload)
    var = payload.raw["var"] if source == "raw.X" else payload.var
    gene_source = _genes_source(var)
    columns = {role: _resolve_column(payload.obs, role) for role in _ALIASES}
    normalized_obs = payload.obs.copy(deep=True)
    for original in columns.values():
        if original:
            normalized_obs[original] = normalized_obs[original].astype("string").str.strip()
    condition_col = columns["condition"]
    conditions = sorted(map(str, normalized_obs[condition_col].dropna().unique())) if condition_col else []
    contrast = _contrast(conditions)
    grouping_requested = _uses_observed_groups(question)
    targets = _target_eligibility(normalized_obs, columns, contrast) if grouping_requested else []
    matrix = payload.layers["counts"] if source == "layers.counts" else payload.raw["X"] if source == "raw.X" else payload.X
    genes = var[gene_source[4:]].astype(str).str.strip().to_numpy() if gene_source != "var.index" else var.index.astype(str).str.strip().to_numpy()
    detected = np.asarray((matrix > 0).sum(axis=1)).ravel()
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if "mt" in var:
        if var["mt"].dtype != bool or var["mt"].isna().any():
            raise ValueError("var.mt must contain complete boolean mitochondrial annotations")
        mito = var["mt"].to_numpy(dtype=bool)
    else:
        mito = np.char.startswith(np.char.lower(genes.astype(str)), "mt-")
    mito_percent = None
    if mito.any():
        mito_percent = np.divide(np.asarray(matrix[:, mito].sum(axis=1)).ravel() * 100., totals,
                                 out=np.zeros(len(totals)), where=totals > 0)
    logs = np.log1p(detected)
    lower = np.median(logs) - 3 * 1.4826 * np.median(np.abs(logs - np.median(logs)))
    min_genes = max(1, int(np.floor(np.expm1(lower))))
    max_mito = 100.0
    if mito_percent is not None:
        median = np.median(mito_percent)
        max_mito = min(100., max(5., float(median + 3 * 1.4826 * np.median(np.abs(mito_percent - median)))))
    batch_col = columns["batch"]
    batches = sorted(map(str, normalized_obs[batch_col].dropna().unique())) if batch_col else []
    batch_complete = bool(batch_col and _complete_labels(normalized_obs[batch_col]))
    confounded = False
    if batch_col and condition_col and len(conditions) > 1:
        design = pd.crosstab(normalized_obs[batch_col], normalized_obs[condition_col])
        confounded = bool((design.gt(0).sum(axis=1) == 1).all())
    spatial, _ = _observed_spatial_coordinates(payload)
    section_col = columns["section"]
    return {
        "n_cells": payload.n_obs, "n_features": int(matrix.shape[1]), "counts_source": source,
        "counts_identity_basis": "Observed finite nonnegative integer values; this check does not establish experimental provenance",
        "gene_name_source": gene_source, "obs_column_map": columns, "conditions": conditions,
        "condition_metadata_complete": bool(condition_col and _complete_labels(payload.obs[condition_col])),
        "contrast": contrast, "use_observed_cell_types": grouping_requested, "targets": targets,
        "observed_cell_types": sorted(map(str, normalized_obs[columns["cell_type"]].dropna().unique())) if grouping_requested and columns["cell_type"] else [],
        "rna_velocity": _velocity_layers(payload, source),
        "spatial": spatial, "sections": sorted(map(str, normalized_obs[section_col].dropna().unique())) if section_col else [],
        "section_metadata_complete": bool(section_col and _complete_labels(normalized_obs[section_col])),
        "batches": batches, "batch_metadata_complete": batch_complete, "batch_condition_confounded": confounded,
        "qc": {"min_genes": min_genes, "max_mito_pct": max_mito, "n_mitochondrial_genes": int(mito.sum()),
               "rule": "Lower 3 scaled MAD on log1p detected genes; upper 3 scaled MAD on mitochondrial percent with 5% floor; unavailable mitochondrial symbols disable that filter"},
    }


def prepare_input(payload: SCData, settings: Mapping[str, Any]) -> tuple[SCData, dict[str, Any]]:
    """Prepare one in-memory input while quarantining prior analysis outputs.

    ``settings`` is the internally generated autonomous_analysis receipt. No
    caller-selected parameters are accepted by the public analyze boundary.
    """
    observation = observe_input(payload, question=str(settings["question"]))
    expected = settings.get("observations")
    if expected is not None and observation != expected:
        raise ValueError("Input observations changed between autonomous planning and import")
    source = observation["counts_source"]
    if source == "raw.X":
        prepared = SCData(payload.raw["X"], payload.obs, payload.raw["var"], uns=payload.uns)
    else:
        prepared = payload.copy()
        prepared.X = (payload.layers["counts"] if source == "layers.counts" else payload.X).copy()
    prepared.layers = {"counts": prepared.X.copy()}
    if observation["rna_velocity"]["available"]:
        prepared.layers.update({name: payload.layers[name].copy() for name in _KINETIC_LAYERS})
    gene_source = observation["gene_name_source"]
    prepared.var["gene_name"] = prepared.var.index.astype(str).str.strip() if gene_source == "var.index" else prepared.var[gene_source[4:]].astype(str).str.strip().to_numpy()
    removed = [f"layers.{name}" for name in payload.layers if name not in prepared.layers]
    removed.extend(f"varm.{name}" for name in payload.varm)
    removed.extend(f"varp.{name}" for name in payload.varp)
    if payload.raw is not None:
        removed.append("raw")
    for name in list(prepared.obs):
        cell_label = name in _ALIASES["cell_type"]
        if _derived_observation_field(name) or (cell_label and not observation["use_observed_cell_types"]):
            removed.append(f"obs.{name}")
            del prepared.obs[name]
    for name in list(prepared.var):
        compact = re.sub(r"[^a-z0-9]", "", str(name).casefold())
        if str(name).casefold() in {"variance_log1p", "means", "dispersions", "dispersions_norm"} or compact.startswith("highlyvariable") or str(name).casefold().startswith(("velocity_", "fit_")):
            removed.append(f"var.{name}")
            del prepared.var[name]
    _, coordinates = _observed_spatial_coordinates(payload)
    retained_coordinates = [name for name in ("spatial", "spatial_coords") if name in payload.obsm] if coordinates is not None else []
    removed.extend(f"obsm.{name}" for name in payload.obsm if name not in retained_coordinates)
    removed.extend(f"obsp.{name}" for name in payload.obsp)
    prepared.obsm, prepared.obsp = {}, {}
    if coordinates is not None:
        prepared.obsm.update({name: np.asarray(payload.obsm[name]).copy() for name in retained_coordinates})
        prepared.obsm.setdefault("spatial", coordinates.copy())
    prepared.varm, prepared.varp = {}, {}
    for name in list(prepared.uns):
        if _derived_unstructured_field(name):
            removed.append(f"uns.{name}")
            del prepared.uns[name]
    for canonical, original in observation["obs_column_map"].items():
        if original and (canonical != "cell_type" or observation["use_observed_cell_types"]):
            prepared.obs[canonical] = payload.obs[original].astype("string").str.strip().to_numpy()
    prepared.raw = None
    receipt = {"policy_version": POLICY_VERSION, "counts_source": source, "gene_name_source": gene_source,
               "obs_column_map": observation["obs_column_map"], "removed_analysis_fields": sorted(set(removed)),
               "removal_reason": "Prior embeddings, graphs, labels and analysis choices cannot supply the autonomous exploration answer",
               "cell_type_source": observation["obs_column_map"]["cell_type"] if observation["use_observed_cell_types"] else None,
               "cell_type_interpretation": "Input observational grouping requested by the scientific question; not autonomous annotation" if observation["use_observed_cell_types"] else "No input cell-type labels used",
               "retained_kinetic_layers": [name for name in _KINETIC_LAYERS if name in prepared.layers],
               "retained_measured_coordinate_fields": sorted(prepared.obsm),
               "source_cells_and_features_preserved": source != "raw.X", "raw_feature_space_used": source == "raw.X"}
    prepared.uns["autonomous_input_preparation"] = deepcopy(receipt)
    return prepared, receipt


def analyze_study(*, question: str, data: str | Path, species: str, tissue: str,
                  run_dir: str | Path, threads: int = 1, max_runtime_hours: float = 12.0) -> dict[str, Any]:
    """Analyze through native run_study using only question/data/context inputs."""
    from eacbp.application import study_service
    from eacbp.application.qc_sensitivity import freeze_qc_sensitivity_plan
    from eacbp.application.qc_sensitivity_execution import execute_qc_sensitivity
    from eacbp.schemas.runtime import RunConfig

    started_at = time.monotonic()

    for label, value in (("question", question), ("species", species), ("tissue", tissue)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{label} must be a nonempty string")
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError("threads must be a positive integer")
    if isinstance(max_runtime_hours, bool) or not math.isfinite(max_runtime_hours) or max_runtime_hours <= 0:
        raise ValueError("max_runtime_hours must be finite and positive")
    source_path = Path(data).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"data must be an existing h5ad file: {source_path}")
    if Path(run_dir).expanduser().resolve().exists():
        raise ValueError("run-dir already exists; choose a new isolated directory")
    before = study_service._hash_file(source_path)
    payload = study_service._import_h5ad(source_path)
    observation = observe_input(payload, question=question)
    after = study_service._hash_file(source_path)
    if before != after:
        raise ValueError("Input h5ad changed during autonomous observation")
    requirement_contract = freeze_question_requirements(question=question, observations=observation, source_sha256=before)
    eligible = [target for target in observation["targets"] if target["eligible"]]
    has_contrast = observation["contrast"] is not None and observation["condition_metadata_complete"]
    comparison_requested = bool(re.search(r"\b(?:compare|comparison|response|effect|differential|change|changes)\b|对比|比较|响应|差异|影响|处理|刺激", question, re.IGNORECASE))
    response = has_contrast and comparison_requested and observation["use_observed_cell_types"] and bool(eligible)
    gaps = []
    limitations = []
    if comparison_requested and not has_contrast:
        gaps.append("The comparison question lacks a complete, unambiguous two-condition experimental design; no condition contrast was invented")
    if has_contrast and comparison_requested and not observation["use_observed_cell_types"]:
        gaps.append("A cell-population contrast requires explicit question context permitting input observed cell-type groups; no mixed-cell replacement estimand was used")
    if has_contrast and comparison_requested and observation["use_observed_cell_types"] and not eligible:
        gaps.append("No observed cell type has complete metadata and at least two qualified biological donors per condition")
    if comparison_requested and observation["use_observed_cell_types"]:
        gaps.extend(f"Requested population {target['cell_type']!r} is not estimable: {target['reason']}"
                    for target in observation["targets"] if not target["eligible"])
    trajectory = bool(re.search(r"trajectory|differentiation|development|lineage|pseudotime|rna\s*velocity|RNA速度|轨迹|发育|分化|谱系|拟时", question, re.IGNORECASE))
    velocity = trajectory and observation["rna_velocity"]["available"]
    if velocity:
        limitations.append("Stochastic RNA velocity infers kinetic direction under model assumptions; developmental time is not observed, and direction does not establish causality")
    elif trajectory:
        gaps.append("No independently justified root cell was provided or inferred; directed DPT cannot be estimated from expression alone")
        gaps.append("RNA velocity is not estimable: " + observation["rna_velocity"]["reason"])
    identity = hashlib.sha256((before + question + species + tissue).encode("utf-8")).hexdigest()[:16]
    manifest = StudyManifest(
        study_id=f"autonomous_{identity}", title=question.strip(),
        biological_design=BiologicalDesign(species=species.strip(), tissue=tissue.strip(),
                                          conditions=observation["conditions"],
                                          target_cell_types=[target["cell_type"] for target in eligible] if response else []),
        experimental_design=ExperimentalDesign(biological_unit="donor" if observation["obs_column_map"]["donor"] else "unknown",
                                              total_samples=len(set(payload.obs[observation["obs_column_map"]["donor"]].dropna().astype(str))) if observation["obs_column_map"]["donor"] else 0,
                                              batches=observation["batches"]),
        constraints=Constraints(max_runtime_hours=float(max_runtime_hours), gpu_allowed=False, min_biological_replicates=2),
    )
    integration = "harmonypy_v1" if len(observation["batches"]) > 1 and observation["batch_metadata_complete"] and not observation["batch_condition_confounded"] else "no_correction_v1"
    decisions = [
        {"decision": "analysis_intent", "value": "condition_comparison" if response else "exploration", "reason": "Observed design and scientific question, without dataset-name dispatch"},
        {"decision": "input_preparation", "value": observation["counts_source"], "reason": "Select a validated integer matrix; quarantine prior analysis outputs"},
        {"decision": "quality_control", "value": observation["qc"], "reason": "Robust distribution policy evaluated on this input before reference evaluation"},
        {"decision": "integration", "value": integration, "reason": "Correct only complete multiple technical batches without complete condition confounding"},
        {"decision": "population_targets", "value": eligible if response else [], "reason": "Every eligible observed population, without manually selecting a cell type"},
        {"decision": "statistics", "value": "pydeseq2_pseudobulk_v1 + donor_sensitivity" if response else "omitted", "reason": "Count-based biological-donor inference and automatic donor robustness check" if response else "No qualified requested population-level contrast"},
        {"decision": "trajectory", "value": "rna_velocity" if velocity else "omitted", "reason": "Requested kinetic direction has validated spliced/unspliced inputs; infer stochastic direction without claiming observed time" if velocity else "No independent root evidence; no stored input root is reused"},
        {"decision": "population_markers", "value": "cluster_markers", "reason": "Native marker discovery on newly computed clusters; no expected gene list or reference annotation is supplied"},
    ]
    settings = {"policy_version": POLICY_VERSION, "autonomy_kind": "deterministic_rule_driven", "question": question,
                "source_sha256": before, "observations": observation, "decisions": decisions, "scientific_gaps": gaps,
                "requirement_contract": requirement_contract,
                "limitations": limitations,
                "intent": "condition_comparison" if response else "exploration",
                "input_preparation": {key: observation[key] for key in ("counts_source", "gene_name_source", "obs_column_map")}}
    config = {"mode": "real", "method_profile": "standard", "threads": threads, "autonomous_analysis": settings,
              "analysis_extensions": {"cluster_markers": True, "qc_review": True},
              "method_overrides": {"integration": integration},
              "capability_parameters": {"qc": {key: observation["qc"][key] for key in ("min_genes", "max_mito_pct")},
                                        "normalization": {"n_top_genes": min(2000, observation["n_features"])},
                                        "integration": {"n_components": min(30, payload.n_obs - 1, observation["n_features"] - 1)}}}
    if velocity:
        config["analysis_extensions"]["rna_velocity"] = True
    if response:
        contrast = observation["contrast"]
        config["advanced_analysis"] = True
        config["analysis_extensions"]["donor_sensitivity"] = True
        config["capability_parameters"]["deg"] = {"condition_col": "condition", "donor_col": "donor", "condition_a": contrast["tested"], "condition_b": contrast["reference"], "n_cpus": threads, "min_donors": 2}
        config["target_parameters"] = {target["cell_type"]: {"deg": {"paired": target["paired"]}, "donor_sensitivity": {"paired": target["paired"]}} for target in eligible}
        # Single-target plans do not expand target_parameters in the legacy DAG.
        if len(eligible) == 1:
            config["capability_parameters"]["deg"]["paired"] = eligible[0]["paired"]
            config["capability_parameters"]["donor_sensitivity"] = {"paired": eligible[0]["paired"]}
            # A single-target plan retains legacy URI names, but its admitted
            # results still need the exact requested population identity.
            for capability in ("subset_cells", "deg", "donor_sensitivity"):
                config["capability_parameters"].setdefault(capability, {})["target_cell_type"] = eligible[0]["cell_type"]
    # Freeze the validated runtime configuration, including defaults, before
    # computing anything. Branches can change only the two declared QC values.
    config = RunConfig.from_mapping(study_service.normalise_config_paths(config, base_dir=Path.cwd())).to_mapping()
    sensitivity_plan = freeze_qc_sensitivity_plan(payload=payload, observations=observation,
                                                source_sha256=before, baseline_config=config)
    config["autonomous_analysis"]["qc_sensitivity_plan"] = sensitivity_plan
    group_column = sensitivity_plan["group_retention"]["source_column"]
    observed_groups = payload.obs[group_column].copy() if group_column else None
    del payload
    summary = study_service.run_study(manifest=manifest, config=config, data=source_path, run_dir=run_dir)
    return execute_qc_sensitivity(plan=sensitivity_plan, manifest=manifest, config=config,
                                 data=source_path, run_dir=Path(run_dir).expanduser().resolve(),
                                 baseline_summary=summary, observed_groups=observed_groups,
                                 started_at=started_at)
