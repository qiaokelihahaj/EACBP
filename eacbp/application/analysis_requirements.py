"""Frozen scientific-question scope and post-execution coverage checks.

This policy records requests, not scientific answers or workflow configuration.
An unavailable method never makes a requested dimension disappear. Callers must
supply only durably admitted results to the coverage evaluator; a successful
task alone is not independent admission. Legacy runs without this contract must
retain their original completion policy rather than acquiring a new contract.
"""

from __future__ import annotations

from copy import deepcopy
import hashlib
import itertools
import json
import re
from typing import Any, Iterable, Mapping


REQUIREMENT_SCHEMA_VERSION = 1
_CONTROLS = {"ctrl", "control", "con", "untreated", "unstim", "unstimulated", "baseline",
             "vehicle", "mock", "wt", "normal", "healthy"}
_PATTERNS = {
    "comparison": r"\b(?:compare|comparison|response|effect|differential|changes?)\b|对比|比较|响应|差异|影响|处理|刺激",
    "abundance": r"\b(?:abundance|proportions?|composition|frequency|frequencies)\b|丰度|比例|组成|占比",
    "pathway": r"\b(?:pathways?|programs?|gene\s+sets?|enrichment)\b|通路|程序|基因集|富集",
    "regulator": r"\b(?:regulators?|regulatory|regulons?|transcription\s+factors?|TFs?)\b|调控|转录因子",
    "communication": r"\b(?:communication|signaling|signalling|ligand|receptor|crosstalk)\b|通讯|通信|配体|受体|细胞间信号",
    "spatial": r"\b(?:spatial|niches?|microenvironments?|tissue\s+domains?)\b|空间|生态位|微环境|组织域",
    "replication": r"\b(?:replicated|replicates?|sections?|slides?|donors?)\b|重复|切片|供体",
    "multi_condition": r"\b(?:multi[ -]condition|multiple\s+conditions?|several\s+conditions?|factorial|interaction\s+effects?|condition\s+interactions?|(?:joint|omnibus)\s+(?:models?|tests?|effects?))\b|多条件|多个条件|条件交互|处理交互|交互效应|析因|联合模型|整体检验",
    "trajectory": r"trajectory|trajectories|differentiation|development|lineage|pseudotime|rna\s*velocity|RNA速度|轨迹|发育|分化|谱系|拟时",
    "branching": r"\b(?:branching|branches|bifurcation|bifurcations)\b|分支|分叉",
    "fate": r"\b(?:fates?|terminal\s+states?|absorption\s+probabilities)\b|命运|终末状态|终末细胞|吸收概率",
    "dynamic_kinetics": r"\b(?:dynamical|dynamic\s+kinetics?|non[ -]steady[ -]state|latent\s+time|recover\s+dynamics)\b|动态动力学|动力学模型|非稳态|潜在时间",
    "robustness": r"\b(?:robustness|robust|sensitivity|leave[ -]one[ -]out|reproducibility|stability)\b|稳健|鲁棒|敏感性|留一|稳定性|可重复",
}


def _fingerprint(value: Mapping[str, Any]) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
                                     separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def requested_dimensions(question: str) -> list[str]:
    """Recognize the supported request vocabulary without dataset dispatch.

    This bounded parser is intentionally not a claim that arbitrary prose has
    been understood. Dimensions independently accumulate in compound questions.
    """
    if not isinstance(question, str) or not question.strip():
        raise ValueError("question must be a nonempty string")
    return [name for name, pattern in _PATTERNS.items() if re.search(pattern, question, re.IGNORECASE)]


def _requested_contrasts(question: str, observations: Mapping[str, Any]) -> list[dict[str, Any]]:
    levels = sorted(set(map(str, observations.get("conditions", []))))
    if not observations.get("condition_metadata_complete"):
        return [{"tested": None, "reference": None, "resolved": False,
                 "reason": "Complete observed condition metadata are required"}]
    binary = observations.get("contrast")
    if len(levels) == 2 and binary:
        return [{**deepcopy(binary), "resolved": True}]
    # More than two observed levels do not authorize an arbitrary pair or an
    # arbitrary omnibus model. Explicit all-pairs/control requests name their
    # estimands using only the observed condition universe.
    if len(levels) > 2 and re.search(r"all\s+pairwise|pairwise\s+comparisons?|两两|所有.*配对比较", question, re.IGNORECASE):
        return [{"tested": b, "reference": a, "resolved": True,
                 "reason": "Explicit all-pairwise request; lexical orientation is noncausal"}
                for a, b in itertools.combinations(levels, 2)]
    controls = [level for level in levels if level.casefold() in _CONTROLS]
    if len(levels) > 2 and len(controls) == 1 and re.search(
            r"(?:each|every|all)\b.*(?:versus|vs\.?|against|relative\s+to).*\b(?:control|baseline)|(?:各|每|所有).*对照", question, re.IGNORECASE):
        return [{"tested": level, "reference": controls[0], "resolved": True,
                 "reason": "Explicit comparison of every observed condition to the uniquely observed control"}
                for level in levels if level != controls[0]]
    return [{"tested": None, "reference": None, "resolved": False,
             "reason": "The question does not resolve a contrast from the observed condition universe"}]


def freeze_question_requirements(*, question: str, observations: Mapping[str, Any],
                                 source_sha256: str | None = None) -> dict[str, Any]:
    """Freeze all requested populations and dimensions before task eligibility.

    Capabilities listed here express required outputs. Their presence does not
    imply that autonomous planning or an installed method supports the request.
    """
    dimensions = requested_dimensions(question)
    conditions = sorted(set(map(str, observations.get("conditions", []))))
    populations = sorted(set(map(str, observations.get("observed_cell_types", []))))
    if not populations and observations.get("use_observed_cell_types"):
        populations = sorted({str(item["cell_type"]) for item in observations.get("targets", [])})
    comparison = "comparison" in dimensions
    contrasts = _requested_contrasts(question, observations) if comparison else []
    targets = {str(item["cell_type"]): item for item in observations.get("targets", [])}
    requirements: list[dict[str, Any]] = []

    def add(dimension, capabilities, *, target=None, contrast=None, prerequisites=(), binding=None):
        requirement = {"dimension": dimension, "capabilities": list(capabilities),
                       "target_cell_type": target, "contrast": deepcopy(contrast),
                       "prerequisite_gaps": list(prerequisites), "binding": dict(binding or {})}
        requirement["requirement_id"] = "R_" + _fingerprint(requirement)[:20]
        requirements.append(requirement)

    for capability in ("dataset_audit", "qc", "normalization", "integration", "clustering", "cluster_markers"):
        add("exploration", [capability])

    population_scopes = populations if observations.get("use_observed_cell_types") and populations else [None]
    for contrast in contrasts:
        for population in population_scopes:
            gaps = []
            if population is None:
                gaps.append("Population comparison requires question-authorized observed population labels")
            if not observations.get("obs_column_map", {}).get("donor"):
                gaps.append("Biological donor metadata are unavailable; sample/library identifiers are not donor identities")
            if not contrast["resolved"]:
                gaps.append(contrast["reason"])
            # The native binary eligibility observation applies only to this
            # exact contrast; new multi-condition models must audit their own
            # observed donor design rather than reuse a binary qualification.
            if population in targets and not targets[population].get("eligible"):
                gaps.append(str(targets[population].get("reason", "Population donor design is not estimable")))
            add("population_expression", ["subset_cells", "deg"], target=population,
                contrast=contrast, prerequisites=gaps)
            add("donor_robustness", ["donor_sensitivity"], target=population,
                contrast=contrast, prerequisites=gaps)

    for dimension, capabilities, binding in (
            ("abundance", ["differential_abundance"], {}),
            ("pathway", ["functional_activity"], {"network_kind": "pathway"}),
            ("regulator", ["functional_activity"], {"network_kind": "regulon"}),
            ("communication", ["liana_communication"], {})):
        if dimension not in dimensions:
            continue
        scopes = population_scopes if dimension in {"pathway", "regulator"} else [None]
        for contrast in contrasts or [None]:
            for population in scopes:
                gaps = [contrast["reason"]] if contrast is not None and not contrast["resolved"] else []
                add(dimension, capabilities, target=population, contrast=contrast,
                    prerequisites=gaps, binding=binding)

    if "spatial" in dimensions:
        spatial = observations.get("spatial", {})
        gaps = [] if spatial.get("available") else [str(spatial.get("reason", "Measured spatial coordinates are unavailable"))]
        add("spatial", ["spatial_domain", "spatial_deg"], prerequisites=gaps)
        if "replication" in dimensions:
            section = observations.get("sections", [])
            replication_gaps = gaps + ([] if len(section) > 1 and observations.get("section_metadata_complete") else ["Replicated spatial inference requires complete observed section identities"])
            add("replicated_spatial", ["replicated_spatial_inference"], prerequisites=replication_gaps)
    joint_model_requested = bool(re.search(r"\bfactorial\b|interaction\s+effects?|condition\s+interactions?|(?:joint|omnibus)\s+(?:models?|tests?|effects?)|条件交互|处理交互|交互效应|析因|联合模型|整体检验", question, re.IGNORECASE))
    if joint_model_requested or ("multi_condition" in dimensions or comparison and len(conditions) > 2) and not (contrasts and all(item["resolved"] for item in contrasts)):
        # A dedicated requirement prevents independent pairwise DEG results
        # from standing in for a requested joint/factorial conclusion.
        add("multi_condition", ["multi_condition_inference"],
            prerequisites=[] if len(conditions) >= (2 if joint_model_requested else 3) else ["The requested model lacks enough observed condition levels"])
    if "trajectory" in dimensions or "branching" in dimensions or "fate" in dimensions or "dynamic_kinetics" in dimensions:
        kinetic = observations.get("rna_velocity", {})
        gaps = [] if kinetic.get("available") else [str(kinetic.get("reason", "Aligned spliced/unspliced counts are unavailable"))]
        add("kinetic_direction", ["rna_velocity"], prerequisites=gaps)
    for dimension, capability in (("branching", "branching_inference"), ("fate", "fate_mapping"),
                                   ("dynamic_kinetics", "dynamical_rna_velocity")):
        if dimension in dimensions:
            add(dimension, [capability])
    if "robustness" in dimensions and not comparison:
        add("analysis_robustness", ["analysis_robustness"])

    contract = {
        "schema_version": REQUIREMENT_SCHEMA_VERSION, "question": question,
        "source_sha256": source_sha256, "requested_dimensions": dimensions,
        "condition_universe": conditions, "population_universe": populations,
        "contrast_universe": contrasts,
        "observed_metadata_roles": deepcopy(observations.get("obs_column_map", {})),
        "population_design_observations": deepcopy(observations.get("targets", [])),
        "requirements": requirements,
        "policy_limitations": ["A deterministic bounded request vocabulary is used; recorded requirements do not establish scientific capability."]}
    contract["contract_id"] = "Q_" + _fingerprint(contract)
    return contract


def freeze_plan_requirements(*, purpose: str, tasks, source_sha256: str) -> dict[str, Any]:
    """Freeze an explicit configured request before runtime branch pruning.

    A task identity binds its scope to the immutable admitted contract. Missing
    primary capabilities remain requirements, rather than disappearing because
    the planner cannot currently supply them.
    """
    from eacbp.orchestrator.research_scope import PURPOSE_CAPABILITIES

    requirements = []
    for task in tasks:
        parameters = task.parameters
        contrast = ({"tested": parameters["condition_a"], "reference": parameters["condition_b"], "resolved": True}
                    if parameters.get("condition_a") is not None and parameters.get("condition_b") is not None else None)
        requirement = {
            "dimension": task.capability, "capabilities": [task.capability],
            "task_id": task.task_id,
            "input_artifacts": list(task.input_artifacts),
            "target_cell_type": task.parameters.get("target_cell_type", task.parameters.get("cell_type") if task.capability == "subset_cells" else None),
            "contrast": contrast, "prerequisite_gaps": [],
            "binding": {key: parameters[key] for key in (
                "condition_col", "donor_col", "cell_type", "target_cell_type", "paired", "design_formula",
                "counts_layer", "root_cell_id", "terminal_states", "network_sha256", "model_sha256",
            ) if key in parameters},
        }
        requirement["requirement_id"] = "R_" + _fingerprint(requirement)[:20]
        requirements.append(requirement)
    present = {task.capability for task in tasks}
    for capability in sorted(PURPOSE_CAPABILITIES[purpose] - present):
        requirement = {"dimension": capability, "capabilities": [capability],
                       "target_cell_type": None, "contrast": None, "binding": {},
                       "prerequisite_gaps": ["The selected purpose has no configured task"]}
        requirement["requirement_id"] = "R_" + _fingerprint(requirement)[:20]
        requirements.append(requirement)
    contract = {"schema_version": REQUIREMENT_SCHEMA_VERSION, "question": purpose,
                "source_sha256": source_sha256, "requirements": requirements,
                "requested_dimensions": [purpose], "policy_limitations": [
                    "Completion records independently admitted requested tasks, not biological truth."]}
    contract["contract_id"] = "Q_" + _fingerprint(contract)
    return contract


def _value(obj: Any, key: str, default=None):
    return obj.get(key, default) if isinstance(obj, Mapping) else getattr(obj, key, default)


def _status(value: Any) -> str | None:
    return getattr(value, "value", value)


def _scope_matches(requirement, result, contract):
    if requirement.get("task_id") is not None:
        if not requirement["task_id"] == _value(result, "task_id") == _value(contract, "task_id"):
            return False
        if requirement.get("input_artifacts") != list(_value(contract, "input_artifacts", [])):
            return False
    parameters = _value(contract, "parameters", {}) or {}
    metrics = _value(result, "metrics", {}) or {}
    target = requirement["target_cell_type"]
    actual_target = parameters.get("cell_type") if _value(result, "capability") == "subset_cells" else parameters.get("target_cell_type")
    actual_target = actual_target or metrics.get("target_cell_type") or metrics.get("subset_cell_type")
    if target is not None and actual_target != target:
        return False
    if target is None and actual_target is not None:
        return False
    contrast = requirement["contrast"]
    if contrast is not None and _value(result, "capability") != "subset_cells":
        actual = (parameters.get("condition_a", metrics.get("condition_a")),
                  parameters.get("condition_b", metrics.get("condition_b")))
        if actual != (contrast["tested"], contrast["reference"]):
            return False
    return all(parameters.get(key, metrics.get(key)) == value for key, value in requirement["binding"].items())


def evaluate_requirement_coverage(requirement_contract: Mapping[str, Any], admitted_results: Iterable[tuple[Any, Any]], *,
                                  task_history: Iterable[Any] = (), pipeline_status: str = "success",
                                  available_capabilities: Iterable[str] | None = None) -> dict[str, Any]:
    """Require exact scientific scope and successful independently admitted outputs.

    ``admitted_results`` pairs a result with its recorded audit contract. The
    caller verifies durable admission, fingerprints and artifact identity first.
    Unsupported/unestimable requirements remain partial, while an audited finite
    estimate with no supported features may satisfy a requested estimate.
    """
    contract = deepcopy(dict(requirement_contract))
    identity = contract.pop("contract_id", None)
    if contract.get("schema_version") != REQUIREMENT_SCHEMA_VERSION or identity != "Q_" + _fingerprint(contract):
        raise ValueError("Scientific requirement contract identity or version is invalid")
    admitted = list(admitted_results)
    history = list(task_history)
    available = None if available_capabilities is None else set(available_capabilities)
    population_outputs: dict[str, set[str]] = {}
    for result, recorded in admitted:
        if _value(result, "capability") == "subset_cells" and _status(_value(result, "status")) == "success":
            parameters = _value(recorded, "parameters", {}) or {}
            metrics = _value(result, "metrics", {}) or {}
            target = parameters.get("cell_type") or metrics.get("subset_cell_type")
            if target is not None:
                population_outputs.setdefault(str(target), set()).update(_value(result, "output_artifacts", []))
    checks, gaps = [], []
    for requirement in contract["requirements"]:
        task_ids, scientific_statuses, missing = [], [], []
        state = "not_estimable" if requirement["prerequisite_gaps"] else "passed"
        for capability in requirement["capabilities"]:
            matches = []
            for result, recorded in admitted:
                if _value(result, "capability") != capability or _status(_value(result, "status")) != "success":
                    continue
                if not _value(result, "method_used") or not _value(result, "output_artifacts") or not _scope_matches(requirement, result, recorded):
                    continue
                scientific = _value(result, "scientific_result")
                scientific_status = _status(_value(scientific, "status"))
                if scientific_status:
                    scientific_statuses.append(scientific_status)
                if capability in {"deg", "donor_sensitivity", "functional_activity", "differential_abundance", "replicated_spatial_inference", "multi_condition_inference"}:
                    if scientific_status not in {"estimated_supported", "estimated_inconclusive"} or _value(scientific, "n_features_tested", 0) < 1:
                        continue
                target = requirement["target_cell_type"]
                if target is not None and capability in {"deg", "donor_sensitivity", "functional_activity"} and not population_outputs.get(target, set()).intersection(_value(result, "input_artifacts", [])):
                    continue
                metrics = _value(result, "metrics", {}) or {}
                if capability == "donor_sensitivity" and (metrics.get("skipped") or metrics.get("all_requested_fits_completed") is not True):
                    continue
                if capability == "liana_communication" and requirement["contrast"] is not None and metrics.get("comparison_status") != "tested":
                    continue
                if not requirement.get("task_id") and capability in {"replicated_spatial_inference", "multi_condition_inference", "branching_inference", "fate_mapping", "dynamical_rna_velocity", "analysis_robustness"} and metrics.get("requirement_dimension") != requirement["dimension"]:
                    continue
                matches.append(_value(result, "task_id"))
            if matches:
                task_ids.extend(matches)
            else:
                missing.append(capability)
        if missing and state == "passed":
            attempted = [_status(_value(result, "status")) for result in history if _value(result, "capability") in missing]
            state = "assumptions_failed" if "assumptions_failed" in scientific_statuses else "not_estimable" if "not_estimable" in scientific_statuses else "unsupported" if available is not None and any(capability not in available for capability in missing) else "failed" if "failed" in attempted else "blocked" if "blocked" in attempted else "missing_or_unadmitted"
        check = {"requirement_id": requirement["requirement_id"], "dimension": requirement["dimension"],
                 "target_cell_type": requirement["target_cell_type"], "contrast": requirement["contrast"],
                 "status": state, "task_ids": sorted(set(task_ids)),
                 "scientific_statuses": sorted(set(scientific_statuses)), "missing_capabilities": missing,
                 "prerequisite_gaps": list(requirement["prerequisite_gaps"])}
        checks.append(check)
        if state != "passed":
            scope = f" for population {requirement['target_cell_type']!r}" if requirement["target_cell_type"] is not None else ""
            contrast = requirement["contrast"]
            if contrast and contrast["resolved"]:
                scope += f" ({contrast['tested']} versus {contrast['reference']})"
            reason = "; ".join(requirement["prerequisite_gaps"]) or ", ".join(missing)
            gaps.append(f"Requested {requirement['dimension']}{scope} is {state}: {reason}")
    return {"analysis_request_status": "complete" if pipeline_status == "success" and not gaps else "partial",
            "requirement_contract_id": identity, "requirement_completion_checks": checks, "scientific_gaps": gaps}
