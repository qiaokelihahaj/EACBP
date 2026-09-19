"""Compose explicitly requested biological analysis extensions."""
from pathlib import Path

from eacbp.schemas.task import TaskContract
from eacbp.orchestrator.resources import pin_resource_files

from eacbp.capabilities.advanced_descriptors import ADVANCED_EXTENSION_DESCRIPTORS

METHODS = {name: descriptor.method for name, descriptor in ADVANCED_EXTENSION_DESCRIPTORS.items()}


def _validate_background_parameters(values):
    """Validate the external CellBender boundary before adding a DAG branch.

    The adapter has the final execution-time checks, but a plan containing an
    implicit input, executable, or destination is already invalid.  Failing at
    planning time keeps the request visible and prevents a later task from
    silently substituting the filtered QC artifact.
    """

    required = ("unfiltered_input_path", "executable", "output_path")
    missing = [key for key in required if values.get(key) is None or not str(values.get(key)).strip()]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(
            "background_removal requires explicit unfiltered_input_path, executable, and output_path; "
            f"missing: {joined}"
        )

    input_path = Path(str(values["unfiltered_input_path"])).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"CellBender unfiltered_input_path does not exist: {input_path}")

    output_path = Path(str(values["output_path"])).expanduser().resolve()
    if output_path.suffix.casefold() != ".h5":
        raise ValueError("CellBender output_path must end with .h5")
    if output_path == input_path:
        raise ValueError("CellBender output_path must differ from unfiltered_input_path")
    if not output_path.parent.is_dir():
        raise FileNotFoundError(
            f"CellBender output_path parent directory does not exist: {output_path.parent}"
        )

    if "run_cwd" in values:
        run_cwd = Path(str(values["run_cwd"])).expanduser().resolve()
        if not run_cwd.is_dir():
            raise FileNotFoundError(f"CellBender run_cwd does not exist: {run_cwd}")


def extend_plan(tasks, manifest, state):
    """Return a DAG whose new tasks preserve the original analysis population."""
    sid = manifest.study_id
    settings = state.get("analysis_extensions", {})
    if not isinstance(settings, dict) or set(settings) - METHODS.keys():
        raise ValueError("analysis_extensions must map supported capability names to true, false, or parameter objects")
    requested = dict(settings)
    if state.get("advanced_analysis"):
        requested.setdefault("donor_sensitivity", True)
    by_cap = {t.capability: t for t in tasks}
    deg = by_cap["deg"]
    if state.get("advanced_analysis"):
        deg.method = "pydeseq2_pseudobulk_v1"
        deg.validation_requirements.append("advanced_statistics_integrity")
        deg.parameters.setdefault("min_donors", max(2, manifest.constraints.min_biological_replicates))
    params = state.get("capability_parameters", {})
    design = dict(deg.parameters)
    audit = by_cap["dataset_audit"]
    for key in ("condition_col", "donor_col"):
        if key in design:
            audit.parameters[key] = design[key]
    for capability, config in requested.items():
        if config is False:
            continue
        if config is not True and not isinstance(config, dict):
            raise ValueError(f"Extension {capability} must be true, false, or a parameter object")
        values = dict(config) if isinstance(config, dict) else {}
        if capability in ("functional_activity", "donor_sensitivity"):
            values = {**design, **values}
        values.update(params.get(capability, {}))
        values.setdefault("species", manifest.biological_design.species)
        values.setdefault("tissue", manifest.biological_design.tissue)
        values.setdefault("random_seed", manifest.reproducibility.random_seed)
        if capability == "background_removal":
            _validate_background_parameters(values)
        values = pin_resource_files(values)
        task_id = "task_ext_" + capability
        output = f"table://{sid}/{capability}/v1"
        input_uri = deg.input_artifacts[0]
        requirements = ["advanced_statistics_integrity"]
        forbidden = ["filter_cells", "recluster", "in_place_mutation"]
        if capability == "background_removal":
            if "qc" not in by_cap or "normalization" not in by_cap:
                raise ValueError(
                    "background_removal requires the base QC and normalization tasks"
                )
            input_uri = by_cap["qc"].input_artifacts[0]
            output = f"adata://{sid}/background_corrected/v1"
            by_cap["qc"].input_artifacts = [output]
            input_layer = by_cap["normalization"].parameters.get("input_layer")
            if input_layer not in (None, "corrected_counts"):
                raise ValueError(
                    "background_removal requires normalization input_layer='corrected_counts'; "
                    f"received {input_layer!r}"
                )
            by_cap["normalization"].parameters["input_layer"] = "corrected_counts"
            requirements = ["advanced_qc_integrity"]
        elif capability == "doublet_detection":
            input_uri = by_cap["qc"].expected_outputs[0]
            output = f"adata://{sid}/doublet_checked/v1"
            by_cap["normalization"].input_artifacts = [output]
            forbidden = ["in_place_mutation"]
            requirements = ["advanced_qc_integrity"]
        elif capability == "cell_annotation":
            input_uri = by_cap["clustering"].expected_outputs[0]
            output = f"adata://{sid}/reference_annotated/v1"
            for existing in tasks:
                existing.input_artifacts = [output if uri == input_uri else uri for uri in existing.input_artifacts]
            forbidden = ["filter_cells", "in_place_mutation"]
            requirements = ["advanced_qc_integrity"]
        elif capability == "liana_communication":
            input_uri = next((t.expected_outputs[0] for t in tasks if t.capability == "cell_annotation"), by_cap["clustering"].expected_outputs[0])
            requirements = ["liana_integrity"]
        expected_outputs = [output]
        if capability == "background_removal":
            expected_outputs.append(f"json://{sid}/cellbender_report/v1")
        task = TaskContract(task_id=task_id, capability=capability, method=METHODS[capability],
                            input_artifacts=[input_uri], expected_outputs=expected_outputs,
                            parameters=values, forbidden_operations=forbidden,
                            validation_requirements=requirements)
        if capability == "doublet_detection":
            task.expected_outputs.append(f"table://{sid}/doublet_scores/v1")
            if not values.get("filter_doublets", False):
                task.forbidden_operations.append("filter_cells")
        if capability == "donor_sensitivity":
            # The primary artifact contains the full-model and donor-specific
            # rows.  A second, stable URI carries the per-gene sensitivity
            # summary (including an explicit skipped artifact when the LOO
            # donor threshold is not met).
            task.expected_outputs.append(f"table://{sid}/donor_sensitivity_summary/v1")
        tasks.append(task)
        if capability == "liana_communication":
            task.expected_outputs.append(f"table://{sid}/liana_condition_comparison/v1")
        descriptor = ADVANCED_EXTENSION_DESCRIPTORS[capability]
        task.validation_requirements = list(descriptor.required_audit_ids)
        task.expected_outputs = descriptor.build_plan(study_id=sid)["expected_outputs"]
    # Annotation can be requested after another extension in JSON ordering.
    annotated = next((t.expected_outputs[0] for t in tasks if t.capability == "cell_annotation"), None)
    if annotated:
        original = by_cap["clustering"].expected_outputs[0]
        for task in tasks:
            if task.capability != "cell_annotation":
                task.input_artifacts = [annotated if uri == original else uri for uri in task.input_artifacts]
    return tasks
