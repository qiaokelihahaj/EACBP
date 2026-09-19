"""Shared plan assembly and a read-only, pre-audit plan preview."""

from copy import deepcopy

from pydantic import ValidationError

from eacbp.capabilities import create_default_capability_registry
from eacbp.capabilities.advanced_descriptors import ADVANCED_EXTENSION_DESCRIPTORS
from eacbp.orchestrator.dag import ComputationalDAGPlanner
from eacbp.orchestrator.router import CapabilityRouter
from eacbp.schemas.runtime import RunConfig


def build_study_tasks(manifest, state, capability_registry):
    """Assemble the same built-in and registered extension DAG for all callers."""
    planner_input = deepcopy(state)
    if "analysis_extensions" in planner_input:
        planner_input["analysis_extensions"] = {
            name: settings for name, settings in planner_input["analysis_extensions"].items()
            if name in ADVANCED_EXTENSION_DESCRIPTORS
        }
    tasks = ComputationalDAGPlanner.build_study_plan(manifest, planner_input)
    tasks.extend(capability_registry.plan_extensions(manifest, state, tasks))
    producers = {}
    for task in tasks:
        for uri in task.expected_outputs:
            if uri in producers:
                raise ValueError(f"Multiple task outputs use the same artifact URI: {uri}")
            producers[uri] = task.task_id
    for task in tasks:
        task.depends_on = list(dict.fromkeys([
            *task.depends_on,
            *(producers[uri] for uri in task.input_artifacts
              if uri in producers and producers[uri] != task.task_id),
        ]))
    return ComputationalDAGPlanner.order_tasks(tasks)


def resolve_task_contract(task, manifest, state, capability_registry, artifact_registry=None, *, router=None):
    """Resolve method, operation bounds and descriptor parameters before hashing."""
    router = router if router is not None else CapabilityRouter(capability_registry)
    task.parameters["study_id"] = manifest.study_id
    task.method = router.resolve_method(task.capability, manifest, state)
    implementation = capability_registry.get(task.capability, task.method)
    operations = getattr(implementation, "contract_operations", None)
    resolver = getattr(implementation, "resolve_contract_operations", None)
    if callable(resolver):
        operations = resolver(task.parameters)
    if operations is not None:
        task.allowed_operations = list(operations)
    return capability_registry.prepare_contract(task, artifact_registry)


def preview_study_plan(manifest, config=None, capability_registry=None):
    """Preview contracts without creating a registry, loading matrices or computing.

    Explicit external resource files may be read and hashed by extension planners.
    No task capability is executed. Dataset audit may later change methods,
    contrasts or omit branches; this document is not a resolved execution receipt.
    """
    state = RunConfig.from_mapping(config).to_mapping()
    state.setdefault("method_profile", "baseline" if state.get("mode") == "demo" else "standard")
    registry = capability_registry if capability_registry is not None else create_default_capability_registry()
    preview = {
        "schema_version": 1,
        "study_id": manifest.study_id,
        "phase": "before_dataset_audit",
        "method_profile": state["method_profile"],
        "tasks": [],
        "external_inputs": [],
        "missing_parameters": [],
        "errors": [],
        "notices": [
            "Methods and contrasts are provisional until dataset audit; unsupported branches may be omitted.",
            "Input matrices, observed cell IDs, optional runtime dependencies and available resources are not validated by this preview.",
            "Runtime-hour, GPU and intermediate-retention constraints are not enforced by the current executor.",
        ],
    }
    if "method_profile" not in RunConfig.from_mapping(config).to_mapping():
        preview["notices"].append(
            "Default profile assumes the supplied mode; run-time metadata identifying simulated input can select baseline instead."
        )
    try:
        tasks = build_study_tasks(manifest.model_copy(deep=True), state, registry)
    except (ValueError, KeyError, OSError) as exc:
        preview["errors"].append({"task_id": None, "message": str(exc)})
        preview["valid"] = False
        return preview
    outputs = {uri for task in tasks for uri in task.expected_outputs}
    preview["external_inputs"] = sorted({uri for task in tasks for uri in task.input_artifacts if uri not in outputs})
    for task in tasks:
        try:
            resolve_task_contract(task, manifest, state, registry)
        except (ValueError, KeyError, TypeError) as exc:
            preview["errors"].append({"task_id": task.task_id, "message": str(exc)})
            if isinstance(exc, ValidationError):
                for error in exc.errors(include_url=False):
                    if error["type"] == "missing":
                        preview["missing_parameters"].append({
                            "task_id": task.task_id, "parameter": ".".join(map(str, error["loc"])),
                            "effect": "task_cannot_execute",
                        })
        if task.method == "scanpy_dpt_v1" and not task.parameters.get("root_cell_id"):
            preview["missing_parameters"].append({
                "task_id": task.task_id, "parameter": "root_cell_id",
                "effect": "trajectory_and_consumers_omitted_after_dataset_audit",
            })
        if task.capability == "quantification" and task.parameters.get("mode") != "demo":
            keys = ("genome_dir", "whitelist_path") if task.method == "starsolo_v1" else ("index_path", "t2g_path")
            for key in keys:
                if not task.parameters.get(key):
                    preview["missing_parameters"].append({
                        "task_id": task.task_id, "parameter": key, "effect": "task_cannot_execute",
                    })
        record = task.model_dump(mode="json")
        record["target_branch"] = task.parameters.get("target_branch")
        record["target_cell_type"] = task.parameters.get("target_cell_type")
        preview["tasks"].append(record)
    preview["valid"] = not preview["errors"] and not any(
        item["effect"] == "task_cannot_execute" for item in preview["missing_parameters"]
    )
    return preview
