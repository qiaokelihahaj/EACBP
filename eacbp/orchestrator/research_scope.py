"""Explicit research purposes select tasks without changing legacy plans."""

PURPOSE_CAPABILITIES = {
    "overview": {"clustering"},
    "deg": {"deg"},
    "abundance": {"differential_abundance"},
    "trajectory": {"trajectory_inference"},
    "batch": {"integration", "clustering"},
    "annotation": {"cell_annotation"},
}


def select_research_tasks(tasks, state, manifest=None):
    purpose = state.get("research_purpose")
    if purpose is None:
        return tasks
    if purpose not in PURPOSE_CAPABILITIES:
        raise ValueError(f"Unknown research purpose: {purpose!r}")
    if state.get("autonomous_analysis") is not None:
        raise ValueError("research_purpose and autonomous_analysis cannot define competing analysis scopes")
    # Explicit extensions remain intentional additions. Merely configuring a
    # capability's parameters does not request that capability's execution.
    requested = PURPOSE_CAPABILITIES[purpose] | {
        name for name, setting in state.get("analysis_extensions", {}).items() if setting is not False
    }
    if (purpose == "deg" and state.get("advanced_analysis")
            and state.get("analysis_extensions", {}).get("donor_sensitivity") is not False):
        requested.add("donor_sensitivity")
    modalities = set(manifest.data.modalities) if manifest is not None else set()
    # Preserve explicitly enabled non-extension branches and their producers.
    if state.get("include_spatial") or "spatial" in modalities:
        requested.update({"spatial_domain", "spatial_deg", "cell_cell_communication"})
    if state.get("include_cci") or state.get("run_cell_communication") or modalities.intersection({"communication", "cci"}):
        requested.add("cell_cell_communication")
    if (state.get("include_knowledge") or state.get("run_knowledge_engine")
            or manifest is not None and (manifest.analysis_policy.prior_guided_analysis or manifest.hypotheses.user_provided)):
        requested.add("knowledge_retrieval")
    if state.get("cellrank_terminal_states") or any(
        values.get("fate_mapping", {}).get("terminal_states")
        for values in state.get("target_parameters", {}).values()
    ):
        requested.add("fate_mapping")
    adapters = set(state.get("adapters", []))
    for alias, capability in (("spacell", "spacell_microenvironment_analysis"),
                              ("geneagent", "gene_function_reasoning"),
                              ("chatcell", "chatcell_dialogue_prediction")):
        if state.get("include_adapters") or alias in adapters or state.get(f"run_{alias}_adapter"):
            requested.add(capability)
    if state.get("run_perturbation") or state.get("include_perturbation") or state.get("perturbation_targets") or "perturbation" in modalities:
        requested.add("genetic_perturbation_simulation")
    if state.get("run_compound_perturbation"):
        requested.add("compound_perturbation_simulation")
    if state.get("full_e2e"):
        raise ValueError("research_purpose cannot be combined with full_e2e; choose an explicit scope")
    by_id = {task.task_id: task for task in tasks}
    missing = requested - {task.capability for task in tasks}
    if missing:
        raise ValueError(f"Research purpose {purpose!r} requires configured capabilities: {sorted(missing)}")
    producers = {uri: task.task_id for task in tasks for uri in task.expected_outputs}
    keep = {task.task_id for task in tasks if task.capability in requested}
    pending = list(keep)
    while pending:
        task = by_id[pending.pop()]
        dependencies = set(task.depends_on) | {producers[uri] for uri in task.input_artifacts if uri in producers}
        for task_id in dependencies - keep:
            if task_id not in by_id:
                raise ValueError(f"Research plan references missing dependency: {task_id}")
            keep.add(task_id)
            pending.append(task_id)
    return [task for task in tasks if task.task_id in keep]
