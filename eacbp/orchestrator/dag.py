"""
Dynamic Computational DAG Generator building auditable task contracts with strict operation bounds.
Dynamically inspects StudyManifest data specifications (spatial coordinates, modalities, prior-guided policy, adapters, perturbations).
"""

from typing import List, Dict, Any, Optional
from eacbp.schemas.study import StudyManifest
from eacbp.schemas.task import TaskContract, RetryPolicy


class ComputationalDAGPlanner:
    """Plans dynamic computational task contracts with explicit allowed/forbidden operation bounds."""

    @staticmethod
    def build_study_plan(manifest: StudyManifest, current_state: Optional[Dict[str, Any]] = None) -> List[TaskContract]:
        current_state = current_state or {}
        sid = manifest.study_id
        tasks = []

        # =========================================================================
        # Stage 0: Dataset Audit
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_001_audit",
            capability="dataset_audit",
            method="sc_audit_v1",
            input_artifacts=[manifest.data.raw_artifact_uri or f"adata://{sid}/raw/v1"],
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
            input_artifacts=[manifest.data.raw_artifact_uri or f"adata://{sid}/raw/v1"],
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
            method="harmony",
            input_artifacts=[f"adata://{sid}/normalized/v2"],
            allowed_operations=["compute_pca", "harmony_pca_alignment", "no_correction_baseline", "evaluate_batch_mixing"],
            forbidden_operations=["filter_cells", "normalize", "recluster"],
            parameters={"n_components": 20},
            expected_outputs=[f"adata://{sid}/integrated/v3"],
            validation_requirements=["batch_mixing_audit", "finite_embedding_check"],
        ))

        # =========================================================================
        # Stage 4: Clustering & Annotation
        # =========================================================================
        tasks.append(TaskContract(
            task_id="task_005_clustering",
            capability="clustering",
            method="leiden_knn_v1",
            input_artifacts=[f"adata://{sid}/integrated/v3"],
            allowed_operations=["build_neighbor_graph", "find_clusters", "annotate_cell_types", "calculate_silhouette"],
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
            method="state_abundance_v1",
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
            method="deg_pseudobulk_v1",
            input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
            allowed_operations=["pseudobulk_aggregation", "cell_level_mannwhitney", "fdr_correction"],
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
            method="paga_dpt",
            input_artifacts=[f"adata://{sid}/microglia_subset/v5"],
            allowed_operations=["build_neighbor_graph", "estimate_pseudotime", "test_root_sensitivity", "find_dynamic_genes"],
            forbidden_operations=["filter_cells", "normalize", "batch_correct", "recluster"],
            parameters={},
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
                expected_outputs=[f"table://{sid}/spacell_niche_table/v1", f"json://{sid}/spacell_summary/v1"],
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
                expected_outputs=[f"table://{sid}/geneagent_pathways/v1", f"json://{sid}/geneagent_summary/v1"],
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
                expected_outputs=[f"table://{sid}/chatcell_dynamics/v1", f"json://{sid}/chatcell_dialogue/v1"],
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
                    "hypotheses": manifest.hypotheses.user_provided or ["DAM TREM2-APOE axis"],
                    "target_genes": current_state.get("target_genes", ["Trem2", "Apoe", "Clec7a", "Tyrobp"]),
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
            target_gene = current_state.get("target_gene", "Trem2")
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
                    "perturbation_type": "knockout",
                    "network_attenuation": 0.35,
                    "delta": 0.95,
                    "random_seed": manifest.reproducibility.random_seed,
                },
                expected_outputs=[f"adata://{sid}/perturbation_ko_trem2/v6", f"table://{sid}/perturbation_shift_trem2/v1"],
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

        return tasks
