"""
Dynamic Computational DAG Generator building auditable task contracts with strict operation bounds.
Dynamically inspects StudyManifest data specifications (spatial coordinates, modalities, prior-guided policy, adapters, perturbations).
"""

from copy import deepcopy
import hashlib
import re
from typing import List, Dict, Any, Optional, Mapping
import networkx as nx
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, RetryPolicy
from eacbp.artifact.uri import ArtifactURI


# These are the analysis capabilities whose input population is the selected
# target cell type.  Shared preprocessing (through annotation) and analyses of
# the complete annotated/spatial matrix remain single tasks.  Keeping the
# population boundary explicit also gives extension planners a stable marker
# when they add target-dependent capabilities later.
_TARGET_BRANCH_CAPABILITIES = frozenset({
    "subset_cells",
    "differential_abundance",
    "deg",
    "trajectory_inference",
    "gene_function_reasoning",
    "chatcell_dialogue_prediction",
    "knowledge_retrieval",
    "genetic_perturbation_simulation",
    "compound_perturbation_simulation",
    "fate_mapping",
    "functional_activity",
    "donor_sensitivity",
})


def _target_slug_base(target: str) -> str:
    """Return the readable, URI-safe stem used for a target branch."""

    value = re.sub(r"[^a-zA-Z0-9]+", "_", str(target).strip().casefold()).strip("_")
    if not value:
        value = "target"
    if value[0].isdigit():
        value = f"target_{value}"
    return value


def target_branch_slug_map(targets: List[str]) -> Dict[str, str]:
    """Map exact requested target names to deterministic branch slugs.

    Target names are preserved in task parameters and provenance.  Slugs only
    identify filesystem/URI and task-id namespaces, so case/punctuation
    variants that would otherwise collide receive a short content hash.  An
    exact duplicate is rejected explicitly because silently running the same
    requested population twice would create misleading duplicate evidence.
    """

    values = [str(target) for target in targets]
    if len(values) != len(set(values)):
        duplicates = sorted({value for value in values if values.count(value) > 1})
        raise ValueError(f"Duplicate target cell types are ambiguous: {duplicates}")
    grouped: Dict[str, List[str]] = {}
    for value in values:
        grouped.setdefault(_target_slug_base(value), []).append(value)
    result: Dict[str, str] = {}
    used = set()
    for value in values:
        stem = _target_slug_base(value)
        if len(grouped[stem]) > 1:
            digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:8]
            slug = f"{stem}_{digest}"
        else:
            slug = stem
        if slug in used:
            # Hash collisions are extraordinarily unlikely, but keep the
            # contract deterministic and fail rather than merging branches.
            raise ValueError(f"Target branch slug collision for {value!r}: {slug!r}")
        used.add(slug)
        result[value] = slug
    return result


def _branch_uri(uri: str, slug: str) -> str:
    """Namespace an artifact URI below its artifact name for one target."""

    parsed = ArtifactURI.parse(uri)
    return ArtifactURI(
        parsed.scheme,
        parsed.study_id,
        f"{parsed.name}/{slug}",
        parsed.version,
    ).to_string()


def _is_target_branch_task(task: TaskContract) -> bool:
    """Identify base tasks that must be replicated per target population."""

    if task.capability in _TARGET_BRANCH_CAPABILITIES:
        return True
    if task.parameters.get("target_cell_type") or task.parameters.get("target_dependent"):
        return True
    return False


def _target_parameter_override(current_state: Mapping[str, Any], target: str, slug: str, capability: str) -> Dict[str, Any]:
    """Resolve optional per-target capability settings without changing legacy state."""

    overrides = current_state.get("target_parameters", {})
    if not isinstance(overrides, Mapping):
        return {}
    raw = overrides.get(target, overrides.get(slug, {}))
    if not isinstance(raw, Mapping):
        return {}
    # Settings are scoped by capability within each target.
    if capability in raw and isinstance(raw.get(capability), Mapping):
        return dict(raw[capability])
    return {}


def _expand_target_branches(tasks: List[TaskContract], targets: List[str], current_state: Mapping[str, Any]) -> List[TaskContract]:
    """Expand only population-dependent tasks into independent target branches.

    The expansion is deliberately URI driven.  Every branch gets a private
    producer namespace and consumers are rewired to that namespace, while
    shared upstream tasks keep their original contracts and are reused by all
    branches.  This makes the normal dependency calculation at the end of the
    planner sufficient to preserve independent failure propagation.
    """

    if len(targets) <= 1:
        return tasks
    slug_by_target = target_branch_slug_map(targets)
    templates = [task for task in tasks if _is_target_branch_task(task)]
    branch_outputs: Dict[str, Dict[str, str]] = {}
    for target in targets:
        slug = slug_by_target[target]
        branch_outputs[target] = {
            uri: _branch_uri(uri, slug)
            for template in templates
            for uri in template.expected_outputs
        }

    expanded: List[TaskContract] = []
    for task in tasks:
        if not _is_target_branch_task(task):
            expanded.append(task)
            continue
        for target in targets:
            slug = slug_by_target[target]
            clone = task.model_copy(deep=True)
            clone.task_id = f"{task.task_id}_{slug}"
            mapping = branch_outputs[target]
            clone.input_artifacts = [mapping.get(uri, uri) for uri in task.input_artifacts]
            clone.expected_outputs = [mapping.get(uri, _branch_uri(uri, slug)) for uri in task.expected_outputs]
            clone.parameters["target_cell_type"] = target
            clone.parameters["target_branch"] = slug
            if clone.capability == "subset_cells":
                clone.parameters["cell_type"] = target
            clone.parameters["fdr_family"] = f"target_cell_type:{target}"
            clone.parameters["target_provenance"] = {
                "target_cell_type": target,
                "target_branch": slug,
                "shared_preprocessing": [
                    "task_000_quant", "task_001_audit", "task_002_qc",
                    "task_003_norm", "task_004_integration", "task_005_clustering",
                ],
                "fdr_scope": "branch_local",
            }
            if clone.capability == "differential_abundance":
                clone.parameters["state_col"] = "microglia_state" if target.casefold() == "microglia" else "sub_state"
            if clone.capability == "knowledge_retrieval":
                clone.parameters["cell_types"] = [target]
            # The root/terminal states are population-specific observations.
            # A single global DPT root must not silently be reused for every
            # subset; callers can provide these through target_parameters.
            override = _target_parameter_override(current_state, target, slug, clone.capability)
            if clone.capability == "trajectory_inference" and current_state.get("method_profile") == "standard":
                if "root_cell_id" not in override:
                    clone.parameters.pop("root_cell_id", None)
                    clone.parameters["target_root_required"] = True
            if clone.capability == "fate_mapping" and "terminal_states" not in override:
                raise ValueError(f"CellRank for {target!r} requires target_parameters terminal_states")
            if override:
                clone.parameters.update(override)
            # Reassert provenance after caller-supplied maps so these fields
            # cannot accidentally be changed into a sibling branch namespace.
            clone.parameters["target_cell_type"] = target
            clone.parameters["target_branch"] = slug
            clone.parameters["fdr_family"] = f"target_cell_type:{target}"
            if clone.capability == "subset_cells":
                if clone.parameters.get("cell_type") != target:
                    raise ValueError(f"Subset cell_type must match its target branch: {target!r}")
            expanded.append(clone)
    return expanded


class ComputationalDAGPlanner:
    """Plans dynamic computational task contracts with explicit allowed/forbidden operation bounds."""

    @staticmethod
    def order_tasks(tasks: List[TaskContract]) -> List[TaskContract]:
        by_id = {task.task_id: task for task in tasks}
        if len(by_id) != len(tasks):
            raise ValueError("Task IDs must be unique within a study.")
        graph = nx.DiGraph()
        graph.add_nodes_from(by_id)
        for task in tasks:
            for parent in task.depends_on:
                if parent not in by_id:
                    raise ValueError(f"Unknown task dependency: {parent}")
                graph.add_edge(parent, task.task_id)
        if not nx.is_directed_acyclic_graph(graph):
            raise ValueError("Task dependencies contain a cycle.")
        return [by_id[task_id] for task_id in nx.topological_sort(graph)]

    @staticmethod
    def build_study_plan(manifest: StudyManifest, current_state: Optional[Dict[str, Any]] = None) -> List[TaskContract]:
        current_state = current_state or {}
        sid = manifest.study_id
        tasks = []

        # =========================================================================
        # Stage 00: FASTQ Upstream Quantification (if raw FASTQ provided)
        # =========================================================================
        has_fastq = (
            manifest.data.has_raw_fastq
            or "FASTQ" in manifest.data.modalities
            or (manifest.data.raw_artifact_uri and manifest.data.raw_artifact_uri.startswith("fastq://"))
            or current_state.get("has_fastq", False)
        )
        if has_fastq:
            quant_in = manifest.data.raw_artifact_uri or f"fastq://{sid}/raw_reads/v1"
            quant_out = f"adata://{sid}/raw/v1"
            tasks.append(TaskContract(
                task_id="task_000_quant",
                capability="quantification",
                method=current_state.get("quant_tool", "kb_python_v1"),
                input_artifacts=[quant_in],
                allowed_operations=["kb_count_alignment", "starsolo_alignment", "sc_quant_demultiplex", "umi_deduplication", "gene_annotation_mapping", "generate_synthetic_counts"],
                forbidden_operations=["filter_cells", "normalize", "recluster"],
                parameters={
                    "mode": current_state.get("mode", "real"),
                    "chemistry": manifest.data.chemistry,
                    "species": manifest.biological_design.species,
                    "target_gene": current_state.get("target_gene", "Kat8"),
                    "random_seed": manifest.reproducibility.random_seed,
                    "threads": current_state.get("threads", 16),
                    "work_dir": current_state.get("work_dir", f"outputs/kb_quant_{sid}"),
                    "index_path": current_state.get("index_path", None),
                    "t2g_path": current_state.get("t2g_path", None),
                    "star_bin": current_state.get("star_bin", "STAR"),
                    "genome_dir": current_state.get("genome_dir"),
                    "whitelist_path": current_state.get("whitelist_path"),
                    "gtf_path": current_state.get("gtf_path"),
                    "num_reads": current_state.get("num_reads", None),
                    "n_cells": current_state.get("n_cells", 1500),
                    "n_genes": current_state.get("n_genes", 400),
                },
                expected_outputs=[quant_out],
                validation_requirements=["finite_expression_check"],
            ))
            base_raw_adata = quant_out
        else:
            base_raw_adata = manifest.data.raw_artifact_uri or f"adata://{sid}/raw/v1"

        # =========================================================================
        # Stage 0: Dataset Audit
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_001_audit",
            capability="dataset_audit",
            method="sc_audit_v1",
            input_artifacts=[base_raw_adata],
            allowed_operations=["audit_metadata", "assess_replication", "assess_batches"],
            forbidden_operations=["filter_cells", "normalize", "recluster"],
            expected_outputs=[f"table://{sid}/dataset_audit/v1"],
            validation_requirements=["sample_count_check", "replicate_check"],
        ))

        # =========================================================================
        # Stage 1: QC Filtering
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_002_qc",
            capability="qc",
            method="sc_qc_v1",
            input_artifacts=[base_raw_adata],
            allowed_operations=["filter_low_quality_cells", "mitochondrial_filtering"],
            forbidden_operations=["normalize", "recluster", "infer_trajectory"],
            parameters={"min_genes": 10, "max_mito_pct": 20.0},
            expected_outputs=[f"adata://{sid}/qc/v1"],
            validation_requirements=["retention_rate_check", "finite_expression_check"],
        ))

        # =========================================================================
        # Stage 2: Normalization & HVG
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_003_norm",
            capability="normalization",
            method="sc_normalize_log1p_v1",
            input_artifacts=[f"adata://{sid}/qc/v1"],
            allowed_operations=["normalize_counts", "log1p_transform", "select_highly_variable_genes"],
            forbidden_operations=["filter_cells", "recluster", "batch_correct"],
            parameters={"target_sum": 10000.0, "n_top_genes": 300},
            expected_outputs=[f"adata://{sid}/normalized/v2"],
            validation_requirements=["finite_expression_check", "hvg_count_check"],
        ))

        # =========================================================================
        # Stage 3: Batch Integration
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_004_integration",
            capability="integration",
            method="batch_mean_centering_v1",
            input_artifacts=[f"adata://{sid}/normalized/v2"],
            allowed_operations=["compute_pca", "batch_mean_centering_pca", "no_correction_baseline", "evaluate_batch_mixing"],
            forbidden_operations=["filter_cells", "normalize", "recluster"],
            parameters={"n_components": 20, "random_seed": manifest.reproducibility.random_seed},
            expected_outputs=[f"adata://{sid}/integrated/v3"],
            validation_requirements=["batch_mixing_audit", "finite_embedding_check"],
        ))

        # =========================================================================
        # Stage 4: Clustering & Annotation
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_005_clustering",
            capability="clustering",
            method="kmeans_marker_embedding_v1",
            input_artifacts=[f"adata://{sid}/integrated/v3"],
            allowed_operations=["pca_first_two_display", "lloyd_kmeans", "marker_score_annotation", "calculate_silhouette"],
            forbidden_operations=["filter_cells", "normalize", "batch_correct"],
            parameters={"k_clusters": 4, "random_seed": manifest.reproducibility.random_seed},
            expected_outputs=[f"adata://{sid}/annotated/v4"],
            validation_requirements=["silhouette_check", "marker_coherence_check"],
        ))

        # =========================================================================
        # Stage 5: Target Subpopulation Subsetting (Microglia)
        # =========================================================================
        target_ct = manifest.biological_design.target_cell_types[0] if manifest.biological_design.target_cell_types else "Microglia"
        tasks.append(TaskContract(
            task_id="task_006_subset",
            capability="subset_cells",
            method="subset_cells_v1",
            input_artifacts=[f"adata://{sid}/annotated/v4"],
            allowed_operations=["subset_cells", "recompute_local_pca", "subcluster_states"],
            forbidden_operations=["normalize", "batch_correct"],
            parameters={"cell_type": target_ct, "obs_key": "cell_type", "random_seed": manifest.reproducibility.random_seed},
            expected_outputs=[f"adata://{sid}/microglia_subset/v5"],
            validation_requirements=["subset_non_empty", "finite_embedding_check"],
        ))

        # =========================================================================
        # Stage 6: Differential Abundance
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_007_abundance",
            capability="differential_abundance",
            method="state_abundance_donor_welch_v1",
            input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
            allowed_operations=["aggregate_donor_proportions", "welch_t_test", "benjamini_hochberg"],
            forbidden_operations=["filter_cells", "recluster", "infer_trajectory"],
            parameters={"state_col": "microglia_state"},
            expected_outputs=[f"table://{sid}/abundance_results/v1"],
            validation_requirements=["fdr_correction_check"],
        ))

        # =========================================================================
        # Stage 7: Differential Expression Analysis (DEG)
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_008_deg",
            capability="deg",
            method="donor_pseudobulk_welch_v1",
            input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
            allowed_operations=["aggregate_counts_by_donor", "library_size_normalize_donor_pseudobulk", "welch_t_test_donor_pseudobulk", "benjamini_hochberg", "compute_cell_level_effect", "welch_t_test_single_cell_exploratory"],
            forbidden_operations=["filter_cells", "recluster", "normalize"],
            parameters={},
            expected_outputs=[f"table://{sid}/deg_results/v1"],
            validation_requirements=["pseudoreplication_audit", "multiple_testing_correction"],
        ))

        # =========================================================================
        # Stage 8: Trajectory Inference & Stability
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_009_trajectory",
            capability="trajectory_inference",
            method="root_distance_pseudotime_v1",
            input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
            allowed_operations=["estimate_root_distance_pseudotime", "test_subsample_stability", "spearman_gene_association", "benjamini_hochberg_correction"],
            forbidden_operations=["filter_cells", "normalize", "batch_correct", "recluster"],
            parameters={"random_seed": manifest.reproducibility.random_seed},
            expected_outputs=[f"table://{sid}/trajectory_results/v1"],
            validation_requirements=["trajectory_subsampling_stability", "dynamic_genes_check"],
        ))

        # =========================================================================
        # Stage 9-11: Spatial Plane (Spatial Domain, Spatial DEG, Spatial CCI)
        # =========================================================================
        has_spatial = (
            manifest.data.has_spatial_coordinates
            or "spatial" in manifest.data.modalities
            or current_state.get("include_spatial", False)
        )

        if has_spatial:
            tasks.append(TaskContract(
                task_id="task_010_spatial_domain",
                capability="spatial_domain",
                method="spatial_domain_knn_v1",
                input_artifacts=[f"adata://{sid}/annotated/v4"],
                allowed_operations=[
                    "validate_spatial_coordinates",
                    "build_spatial_knn",
                    "build_spatial_connectivities",
                    "spatially_smoothed_embedding",
                    "cluster_spatial_domains",
                    "calculate_silhouette",
                    "identify_domains",
                ],
                forbidden_operations=["filter_cells", "normalize", "batch_correct"],
                parameters={"k_neighbors": 6, "n_domains": 4, "smoothing_lambda": 0.3, "random_seed": manifest.reproducibility.random_seed},
                expected_outputs=[f"adata://{sid}/spatial_domains/v6"],
                validation_requirements=["finite_expression_check", "silhouette_check"],
            ))

            tasks.append(TaskContract(
                task_id="task_011_spatial_deg",
                capability="spatial_deg",
                method="spatial_deg_morans_i_v1",
                input_artifacts=[f"adata://{sid}/spatial_domains/v6"],
                allowed_operations=[
                    "build_spatial_connectivities",
                    "calculate_moran_i",
                    "calculate_geary_c",
                    "analytical_significance_test",
                    "benjamini_hochberg_correction",
                    "fdr_correction",
                    "identify_spatial_degs",
                    "calculate_morans_i",
                    "calculate_gearys_c",
                    "extract_spatial_graph",
                ],
                forbidden_operations=["filter_cells", "normalize", "recluster"],
                parameters={"k_neighbors": 6, "min_moran_i": 0.15, "fdr_threshold": 0.05},
                expected_outputs=[f"table://{sid}/spatial_deg/v1"],
                validation_requirements=["multiple_testing_correction", "morans_i_bounds_check"],
            ))

            tasks.append(TaskContract(
                task_id="task_012_spatial_cci",
                capability="cell_cell_communication",
                method="cci_ligand_receptor_v1",
                input_artifacts=[f"adata://{sid}/spatial_domains/v6"],
                allowed_operations=[
                    "load_lr_database",
                    "calculate_spatial_contact_density",
                    "compute_spatial_cci_score",
                    "run_spatial_permutation_test",
                    "extract_spatial_proximity",
                    "match_ligand_receptor_pairs",
                    "score_spatial_interactions",
                    "permutation_testing",
                    "evaluate_cell_cell_communication",
                    "ligand_receptor_cci",
                ],
                forbidden_operations=["filter_cells", "recluster", "normalize"],
                parameters={"fdr_threshold": 0.05, "n_permutations": 200, "random_seed": manifest.reproducibility.random_seed},
                expected_outputs=[f"table://{sid}/spatial_cci/v1"],
                validation_requirements=["fdr_correction_check"],
            ))
        else:
            # Standard non-spatial cell-cell communication
            include_cci = (
                "communication" in manifest.data.modalities
                or "cci" in manifest.data.modalities
                or current_state.get("include_cci", False)
                or current_state.get("run_cell_communication", False)
            )
            if include_cci:
                tasks.append(TaskContract(
                    task_id="task_010_cell_cell_communication",
                    capability="cell_cell_communication",
                    method="cci_ligand_receptor_v1",
                    input_artifacts=[f"adata://{sid}/annotated/v4"],
                    allowed_operations=[
                        "load_lr_database",
                        "calculate_spatial_contact_density",
                        "compute_spatial_cci_score",
                        "run_spatial_permutation_test",
                        "extract_spatial_proximity",
                        "match_ligand_receptor_pairs",
                        "score_spatial_interactions",
                        "permutation_testing",
                        "evaluate_cell_cell_communication",
                        "ligand_receptor_cci",
                    ],
                    forbidden_operations=["filter_cells", "recluster", "normalize"],
                    parameters={"fdr_threshold": 0.05, "n_permutations": 100, "random_seed": manifest.reproducibility.random_seed},
                    expected_outputs=[f"table://{sid}/spatial_cci/v1"],
                    validation_requirements=["fdr_correction_check"],
                ))

        # =========================================================================
        # Stage 12-14: External Agent Adapters (SpaCell, GeneAgent, ChatCell)
        # =========================================================================
        include_adapters = current_state.get("include_adapters", False)
        requested_adapters = current_state.get("adapters", [])

        if has_spatial and (include_adapters or "spacell" in requested_adapters or current_state.get("run_spacell_adapter", False)):
            tasks.append(TaskContract(
                task_id="task_013_spacell_adapter",
                capability="spacell_microenvironment_analysis",
                method="spacell_agent_v1",
                input_artifacts=[f"adata://{sid}/spatial_domains/v6"],
                allowed_operations=[
                    "compute_spatial_neighbors",
                    "spatial_domain_clustering",
                    "analyze_microenvironment",
                    "identify_domains",
                    "evaluate_contact_enrichment",
                    "calculate_morans_i",
                ],
                forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
                parameters={"k_neighbors": 6, "n_clusters": 4, "random_seed": manifest.reproducibility.random_seed},
                expected_outputs=[f"adata://{sid}/spacell_domains/v6_spacell", f"table://{sid}/spacell_niche_metrics/v1", f"json://{sid}/spacell_summary/v1"],
                validation_requirements=["sample_count_check"],
            ))

        if include_adapters or "geneagent" in requested_adapters or current_state.get("run_geneagent_adapter", False):
            tasks.append(TaskContract(
                task_id="task_014_geneagent_adapter",
                capability="gene_function_reasoning",
                method="gene_agent_v1",
                input_artifacts=[f"table://{sid}/deg_results/v1"],
                allowed_operations=[
                    "query_gene_ontology",
                    "map_reactome_pathways",
                    "gene_function_reasoning",
                    "ortholog_lookup",
                ],
                forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
                parameters={"species": manifest.biological_design.species},
                expected_outputs=[f"table://{sid}/gene_agent_pathways/v1", f"json://{sid}/gene_agent_summary/v1"],
                validation_requirements=["fdr_correction_check"],
            ))

        if include_adapters or "chatcell" in requested_adapters or current_state.get("run_chatcell_adapter", False):
            tasks.append(TaskContract(
                task_id="task_015_chatcell_adapter",
                capability="chatcell_dialogue_prediction",
                method="chatcell_agent_v1",
                input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
                allowed_operations=[
                    "query_cell_dialogue",
                    "predict_state_transition",
                    "summarize_phenotype",
                ],
                forbidden_operations=["filter_cells", "recluster", "in_place_mutation"],
                parameters={"target_cell_type": target_ct, "source_state": "control", "target_state": "AD"},
                expected_outputs=[f"table://{sid}/chatcell_state_transitions/v1", f"json://{sid}/chatcell_dialogue/v1"],
                validation_requirements=["finite_expression_check"],
            ))

        # =========================================================================
        # Stage 15: Multi-Source Knowledge Engine
        # =========================================================================
        include_knowledge = (
            current_state.get("run_knowledge_engine", False)
            or current_state.get("include_knowledge", False)
            or manifest.analysis_policy.prior_guided_analysis
            or bool(manifest.hypotheses.user_provided)
            or current_state.get("full_e2e", False)
        )

        if include_knowledge:
            is_prior_guided = manifest.analysis_policy.prior_guided_analysis or bool(manifest.hypotheses.user_provided)
            tasks.append(TaskContract(
                task_id="task_016_knowledge",
                capability="knowledge_retrieval",
                method="knowledge_engine_prior_v1" if is_prior_guided else "knowledge_engine_discovery_v1",
                input_artifacts=[f"table://{sid}/deg_results/v1"],
                allowed_operations=["query_pubmed_literature", "query_gene_ontology_ora", "query_reactome_pathways", "synthesize_knowledge_report"],
                forbidden_operations=["filter_cells", "recluster", "mutate_raw"],
                parameters={
                    "prior_guided": is_prior_guided,
                    "hypotheses": manifest.hypotheses.user_provided,
                    "species": manifest.biological_design.species,
                    "tissue": manifest.biological_design.tissue,
                    "disease": manifest.biological_design.disease,
                    "cell_types": manifest.biological_design.target_cell_types,
                    "target_genes": current_state.get("target_genes", []),
                },
                expected_outputs=[f"table://{sid}/knowledge_evidence/v1", f"json://{sid}/knowledge_report/v1"],
                validation_requirements=["epistemic_tagging_check"],
            ))

        # =========================================================================
        # Stage 16-17: In Silico Perturbation Simulation
        # =========================================================================
        include_perturbation = (
            current_state.get("run_perturbation", False)
            or current_state.get("include_perturbation", False)
            or "perturbation" in manifest.data.modalities
            or current_state.get("full_e2e", False)
            or bool(current_state.get("perturbation_targets"))
        )

        if include_perturbation:
            target_gene = current_state.get("target_gene")
            if not target_gene:
                raise ValueError("Perturbation requires an explicit target_gene in current_state.")
            tasks.append(TaskContract(
                task_id="task_017_genetic_perturbation",
                capability="genetic_perturbation_simulation",
                method="in_silico_crispr_ko_v1",
                input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
                allowed_operations=[
                    "construct_grn_adjacency",
                    "simulate_genetic_perturbation",
                    "propagate_network_shift",
                    "compute_state_reversion",
                    "construct_grn",
                    "compute_grn_propagator",
                    "simulate_knockout",
                    "project_counterfactual_state",
                ],
                forbidden_operations=["recluster", "in_place_mutation"],
                parameters={
                    "target_gene": target_gene,
                    "perturbation_type": current_state.get("perturbation_type", "knockout"),
                    "network_attenuation": 0.35,
                    "delta": 0.95,
                    "random_seed": manifest.reproducibility.random_seed,
                },
                expected_outputs=[f"adata://{sid}/perturbation_{current_state.get('perturbation_type', 'knockout')}_{target_gene.lower()}/v6", f"table://{sid}/perturbation_shift_{target_gene.lower()}/v1"],
                validation_requirements=["perturbation_shift_bounds_check", "matrix_non_empty"],
            ))

            if current_state.get("run_compound_perturbation", False):
                tasks.append(TaskContract(
                    task_id="task_018_compound_perturbation",
                    capability="compound_perturbation_simulation",
                    method="in_silico_compound_response_v1",
                    input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
                    allowed_operations=[
                        "compute_disease_signature",
                        "calculate_cmap_discordance",
                        "simulate_counterfactual_transitions",
                        "compute_transition_matrix",
                        "load_cmap_signatures",
                        "compute_cosine_discordance",
                        "simulate_state_transitions",
                    ],
                    forbidden_operations=["recluster", "in_place_mutation"],
                    parameters={
                        "candidate_compounds": ["Compound_A", "Compound_B", "Compound_C"],
                        "random_seed": manifest.reproducibility.random_seed,
                    },
                    expected_outputs=[f"table://{sid}/compound_reversal_scores/v1", f"adata://{sid}/perturbation_drug_compound_a/v1"],
                    validation_requirements=["state_transition_stochasticity_check"],
                ))

        target_fate_requested = len(manifest.biological_design.target_cell_types) > 1 and any(
            isinstance(settings, Mapping) and settings.get("fate_mapping", {}).get("terminal_states")
            for settings in current_state.get("target_parameters", {}).values()
        )
        if current_state.get("cellrank_terminal_states") or target_fate_requested:
            if current_state.get("method_profile") != "standard":
                raise ValueError("CellRank fate mapping requires method_profile='standard' and an explicit DPT root")
            trajectory = next(t for t in tasks if t.capability == "trajectory_inference")
            trajectory.expected_outputs.append(f"adata://{sid}/diffusion_pseudotime/v1")
            tasks.append(TaskContract(task_id="task_019_cellrank", capability="fate_mapping", method="cellrank_fate_v1",
                input_artifacts=[f"adata://{sid}/diffusion_pseudotime/v1"],
                parameters={"terminal_states": current_state.get("cellrank_terminal_states", {}), "kernel": "pseudotime"},
                expected_outputs=[f"table://{sid}/fate_probabilities/v1"],
                forbidden_operations=["filter_cells", "normalize", "recluster"]))

        # Multiple requested populations share preprocessing but receive
        # independent subset/statistical/consumer branches.  Keep the legacy
        # single-target and broad-study contracts untouched below.
        # A broad study analyzes all annotated cells; never silently select Microglia.
        if not manifest.biological_design.target_cell_types:
            tasks = [task for task in tasks if task.capability != "subset_cells"]
            for task in tasks:
                task.input_artifacts = [uri.replace("/microglia_subset/v5", "/annotated/v4") for uri in task.input_artifacts]
                if task.capability == "differential_abundance":
                    task.parameters["state_col"] = "cluster"
        elif target_ct.lower() != "microglia":
            for task in tasks:
                task.input_artifacts = [uri.replace("/microglia_subset/v5", "/target_subset/v5") for uri in task.input_artifacts]
                task.expected_outputs = [uri.replace("/microglia_subset/v5", "/target_subset/v5") for uri in task.expected_outputs]
                if task.capability == "differential_abundance":
                    task.parameters["state_col"] = "sub_state"
        # Apply explicit per-capability settings to the contract, not only the router.
        for task in tasks:
            if task.capability == "chatcell_dialogue_prediction":
                # This adapter produces a transition table, not an expression matrix.
                task.validation_requirements = ["table_finite_values", "json_valid_payload"]
            global_parameters = current_state.get("capability_parameters", {}).get(task.capability, {})
            if isinstance(global_parameters, Mapping):
                task.parameters.update(global_parameters)
            if task.capability == "deg":
                task.parameters.setdefault("min_replicates", manifest.constraints.min_biological_replicates)
        if current_state.get("advanced_analysis") or current_state.get("analysis_extensions"):
            from eacbp.orchestrator.advanced_plan import extend_plan
            tasks = extend_plan(tasks, manifest, current_state)
        if len(manifest.biological_design.target_cell_types) > 1:
            tasks = _expand_target_branches(tasks, list(manifest.biological_design.target_cell_types), current_state)
        producers = {uri: task.task_id for task in tasks for uri in task.expected_outputs}
        for task in tasks:
            task.depends_on = list(dict.fromkeys(producers[uri] for uri in task.input_artifacts
                                                if uri in producers and producers[uri] != task.task_id))
            if task.capability == "qc":
                task.depends_on.append("task_001_audit")
        return tasks

    @staticmethod
    def adapt_after_audit(tasks, state):
        """Prune unsupported inference branches from audited observed metadata.

        Explicit requests fail with an actionable error; automatic branches are
        omitted and their consumers are omitted transitively.
        """
        if "observed_conditions" not in state:
            return tasks, []
        conditions = state["observed_conditions"]
        contrast_ok = len(conditions) == 2 and state.get("condition_metadata_complete", False)
        reasons = {}
        for task in tasks:
            if task.capability == "trajectory_inference" and state.get("method_profile") == "standard" and not task.parameters.get("root_cell_id"):
                reasons[task.task_id] = "DPT requires an explicit root_cell_id; no root cell was invented"
            if task.capability in {"deg", "differential_abundance", "functional_activity", "donor_sensitivity"}:
                a, b = task.parameters.get("condition_a"), task.parameters.get("condition_b")
                explicit = a is not None or b is not None
                if explicit:
                    if a == b or a not in conditions or b not in conditions:
                        raise ValueError("Explicit contrast must name two distinct observed conditions")
                elif not contrast_ok:
                    reasons[task.task_id] = "No complete, unambiguous two-condition contrast in audited metadata"
                else:
                    task.parameters.update(condition_a=conditions[0], condition_b=conditions[1])
            if task.capability == "differential_abundance" and not state.get("biological_replication_sufficient", False):
                reasons[task.task_id] = "Insufficient independent donors for differential abundance"
            elif task.capability == "differential_abundance" and state.get("paired_donors_observed", False):
                reasons[task.task_id] = "The available independent-donor abundance method cannot model paired donors"
        changed = True
        while changed:
            changed = False
            for task in tasks:
                if task.task_id not in reasons and any(dep in reasons for dep in task.depends_on):
                    reasons[task.task_id] = "Required inference branch was omitted after metadata audit"
                    changed = True
        return [task for task in tasks if task.task_id not in reasons], [
            {"task_id": key, "reason": value} for key, value in reasons.items()]
