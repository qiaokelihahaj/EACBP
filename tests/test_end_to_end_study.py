"""
End-to-End Tests: Complete execution of simulated AD Mouse Brain Single-Cell Studies.
Validates all 6 Planes, Dual DAGs, Task Contracts, Independent Auditing, Multimodal Evidence, and Provenance.
"""

import pytest
import shutil
from pathlib import Path
import pandas as pd
import numpy as np

from eacbp.schemas.artifact import ArtifactType
from eacbp.schemas.study import (
    StudyManifest,
    BiologicalDesign,
    ExperimentalDesign,
    DataSpec,
    Hypotheses,
    AnalysisPolicy,
    ReproducibilityConfig,
)
from eacbp.schemas.evidence import LanguageTier, ClaimType
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.capabilities.spatial import create_synthetic_spatial_ad_study
from eacbp.report.markdown_report import ScientificReportGenerator


def test_full_ad_mouse_study_pipeline(tmp_path):
    """
    Baseline scRNA-seq pipeline execution (Stages 0-8).
    """
    storage_dir = tmp_path / "artifacts_e2e_scrna"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    # 1. User Prompt & Intent Parsing
    prompt = "分析阿尔茨海默病小鼠脑组织单细胞数据，寻找小胶质细胞时空演化轨迹与关键调控通路。"
    study_id = "AD_mouse_study_001"
    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id=study_id)

    assert manifest.study_id == study_id
    assert manifest.biological_design.species == "mus_musculus"

    # 2. Ingest Synthetic Raw Data (12 mice: 6 AD, 6 Ctrl, 2 batches, realistic markers)
    raw_data = SCData.create_synthetic_ad_study(
        n_cells=600,
        n_genes=300,
        n_ad_mice=6,
        n_ctrl_mice=6,
        random_seed=42,
    )

    raw_uri = f"adata://{study_id}/raw/v1"
    manifest.data.raw_artifact_uri = raw_uri

    registry.register(
        uri_str=raw_uri,
        payload=raw_data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id=study_id,
        created_by_task="task_000_ingest",
        operation="raw_data_ingest",
    )

    # 3. Execute Full Study Workflow
    study_summary = orchestrator.run_study(manifest)

    assert study_summary["study_id"] == study_id
    assert study_summary["tasks_executed"] == 9
    assert study_summary["artifacts_created"] >= 9
    assert study_summary["evidence_nodes_count"] >= 3
    assert study_summary["claims_count"] >= 2

    # 4. Check Invariant 1 & 2: Raw data untouched, versions immutable
    assert registry.exists(raw_uri)
    assert registry.exists(f"adata://{study_id}/qc/v1")
    assert registry.exists(f"adata://{study_id}/normalized/v2")
    assert registry.exists(f"adata://{study_id}/integrated/v3")
    assert registry.exists(f"adata://{study_id}/annotated/v4")
    assert registry.exists(f"adata://{study_id}/microglia_subset/v5")
    assert registry.exists(f"table://{study_id}/deg_results/v1")
    assert registry.exists(f"table://{study_id}/trajectory_results/v1")

    # 5. Check Scientific Audits
    assert len(orchestrator.audit_reports) == 9
    for report in orchestrator.audit_reports:
        assert report.overall_passed is True

    # 6. Check Evidence Graph and Synthesized Claims
    claims = orchestrator.evidence_graph.claim_nodes
    assert "C101_microglia_state_transition" in claims
    assert "C102_dam_marker_expression" in claims

    c101 = claims["C101_microglia_state_transition"]
    assert c101.confidence.association > 0.5
    assert c101.confidence.overall > 0.5
    assert len(c101.support_evidence_ids) >= 2

    # 7. Generate Scientific Report with Provenance Chains
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()

    assert "Scientific Study Report" in report_md
    assert "C101_microglia_state_transition" in report_md
    assert "Computational Lineage Graph" in report_md
    assert "Evidence-to-Claim DAG" in report_md


def test_full_multimodal_spatial_and_perturbation_study_pipeline(tmp_path):
    """
    Comprehensive 18-step E2E study across all 6 planes:
    scRNA + Spatial + SpaCell + GeneAgent + ChatCell + Knowledge Engine (Prior-guided) +
    In Silico CRISPR KO Simulation (Trem2) + Compound Simulation + Multimodal 5-Pillar Claims + Provenance Report.
    """
    storage_dir = tmp_path / "artifacts_e2e_multimodal"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    study_id = "AD_spatial_perturbation_study_001"

    # 1. Study Manifest with multi-modal specifications
    manifest = StudyManifest(
        study_id=study_id,
        title="Multimodal Spatial Single-Cell & In Silico Perturbation Study of AD Microglia",
        biological_design=BiologicalDesign(
            species="mus_musculus",
            tissue="cortex",
            disease="Alzheimer",
            conditions=["AD", "control"],
            target_cell_types=["Microglia"],
        ),
        experimental_design=ExperimentalDesign(
            biological_unit="mouse",
            batches=["batch_1", "batch_2"],
            total_samples=12,
            donor_replicates_per_condition={"AD": 6, "control": 6},
        ),
        data=DataSpec(
            modalities=["scRNA", "spatial", "perturbation"],
            has_spatial_coordinates=True,
            has_rna_velocity=False,
        ),
        hypotheses=Hypotheses(
            user_provided=["DAM TREM2-APOE activation axis in plaque microenvironment"],
        ),
        analysis_policy=AnalysisPolicy(
            discovery_mode=False,
            prior_guided_analysis=True,
            prefer_pseudobulk=True,
        ),
        reproducibility=ReproducibilityConfig(
            random_seed=42,
            save_intermediate_artifacts=True,
        ),
    )

    # 2. Ingest Synthetic Spatial AD Dataset (with plaque foci and spatial markers)
    raw_spatial_data = create_synthetic_spatial_ad_study(
        n_cells=600,
        n_genes=100,
        n_plaques=3,
        random_seed=42,
    )

    raw_uri = f"adata://{study_id}/raw/v1"
    manifest.data.raw_artifact_uri = raw_uri

    registry.register(
        uri_str=raw_uri,
        payload=raw_spatial_data.to_dict(),
        artifact_type=ArtifactType.SPATIAL_DATA,
        study_id=study_id,
        created_by_task="task_000_ingest",
        operation="raw_spatial_data_ingest",
    )

    # 3. Execute Full Multi-Plane Study Workflow
    study_summary = orchestrator.run_study(
        manifest,
        current_state={
            "full_e2e": True,
            "include_spatial": True,
            "include_adapters": True,
            "include_knowledge": True,
            "include_perturbation": True,
            "run_compound_perturbation": True,
            "target_gene": "Trem2",
            "target_genes": ["Trem2", "Apoe", "Clec7a", "Tyrobp"],
        }
    )

    # 4. Verify Task DAG Execution Across All 18 Steps
    assert study_summary["study_id"] == study_id
    assert study_summary["tasks_executed"] == 18, f"Expected 18 tasks, got {study_summary['tasks_executed']}"
    assert study_summary["artifacts_created"] >= 18
    assert study_summary["evidence_nodes_count"] >= 8
    assert study_summary["claims_count"] >= 4

    # 5. Check Invariant 1 (Raw data immutability)
    assert registry.exists(raw_uri)
    raw_meta = registry.get_metadata(raw_uri)
    assert raw_meta.type == ArtifactType.SPATIAL_DATA

    # 6. Check Invariant 2 (Content-addressed versioned artifacts)
    # Check that core artifacts from all stages exist in registry
    registered_uris = list(registry.registry.keys())
    assert any("qc/v1" in u for u in registered_uris)
    assert any("normalized/v2" in u for u in registered_uris)
    assert any("integrated/v3" in u for u in registered_uris)
    assert any("annotated/v4" in u for u in registered_uris)
    assert any("microglia_subset/v5" in u for u in registered_uris)
    assert any("deg_results/v1" in u for u in registered_uris)
    assert any("trajectory_results/v1" in u for u in registered_uris)
    assert any("spatial_domains" in u for u in registered_uris)
    assert any("spatial_deg" in u for u in registered_uris)
    assert any("spatial_cci" in u for u in registered_uris)
    assert any("spacell" in u for u in registered_uris)
    assert any("gene_agent" in u or "geneagent" in u for u in registered_uris)
    assert any("chatcell" in u for u in registered_uris)
    assert any("knowledge_evidence" in u for u in registered_uris)
    assert any("knowledge_report" in u for u in registered_uris)
    assert any("perturbation" in u for u in registered_uris)
    assert any("compound" in u for u in registered_uris)

    # 7. Check Invariant 3: Multi-Tier Calibrated Evidence Across 5 Pillars
    evidence_nodes = orchestrator.evidence_graph.evidence_nodes
    evidence_types = {ev.type for ev in evidence_nodes.values()}
    
    # 5 Pillars Verified:
    # Pillar 1: Spatial localization
    # Pillar 2: Pseudobulk DEG
    # Pillar 3: Literature & Pathway
    # Pillar 4: Trajectory
    # Pillar 5: Perturbation
    from eacbp.schemas.evidence import EvidenceType
    assert EvidenceType.SPATIAL_LOCALIZATION in evidence_types
    assert EvidenceType.PSEUDOBULK_DEG in evidence_types or EvidenceType.CELL_LEVEL_DEG in evidence_types
    assert EvidenceType.LITERATURE_SUPPORT in evidence_types or EvidenceType.PATHWAY_ENRICHMENT in evidence_types
    assert EvidenceType.TRAJECTORY_STABILITY in evidence_types
    assert EvidenceType.PERTURBATION in evidence_types

    # 8. Check Synthesized Multimodal Claims
    claims = orchestrator.evidence_graph.claim_nodes
    assert "C101_microglia_state_transition" in claims
    assert "C102_dam_marker_expression" in claims
    assert "C103_knowledge_pathway_convergence" in claims
    assert "C104_in_silico_perturbation_reversal" in claims

    # Check C102: Spatial Plaque Niche Localization
    c102 = claims["C102_dam_marker_expression"]
    assert "spatial" in c102.statement.lower() or "plaque" in c102.statement.lower()
    assert c102.language_tier == LanguageTier.LEVEL_2_STATISTICAL_INFERENCE
    assert c102.causal_status == "observational"
    assert c102.confidence.association > 0.6
    assert c102.confidence.mechanistic > 0.6
    assert c102.confidence.overall > 0.6

    # Check C103: Prior-guided hypothesis testing tag
    c103 = claims["C103_knowledge_pathway_convergence"]
    assert "[PRIOR-GUIDED HYPOTHESIS TESTING]" in c103.statement
    assert c103.language_tier == LanguageTier.LEVEL_4_HYPOTHESIS

    # Check C104: In Silico Perturbation Simulation
    c104 = claims["C104_in_silico_perturbation_reversal"]
    assert c104.language_tier == LanguageTier.LEVEL_4_HYPOTHESIS
    assert c104.causal_status == "in_silico_perturbed"
    assert c104.confidence.causal <= 0.50, f"In silico causal confidence must not exceed 0.50 ceiling, got {c104.confidence.causal}"
    assert c104.confidence.causal > 0.0

    # 9. Check Invariant 4: Independent Scientific Audits
    assert len(orchestrator.audit_reports) == 18
    for report in orchestrator.audit_reports:
        assert report.overall_passed is True, f"Audit failed for task {report.target_task_id}: {[c.message for c in report.checks if not c.passed]}"

    # 10. Generate Complete Provenance Manuscript Report
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()

    assert "Scientific Study Report" in report_md
    assert "[PRIOR-GUIDED HYPOTHESIS TESTING]" in report_md
    assert "## 1. Study Design & Experimental Audit" in report_md
    assert "## 2. Computational Task DAG Execution Summary" in report_md
    assert "## 3. Evidence-Grounded Scientific Claims & Multimodal Calibration" in report_md
    assert "## 4. Scientific Auditor Sign-offs & Independent Verification" in report_md
    assert "## 5. Registered Artifacts & Lineage DAG" in report_md
    assert "## 6. Computational Lineage Graph" in report_md
    assert "## 7. Evidence-to-Claim DAG" in report_md
    assert "C101_microglia_state_transition" in report_md
    assert "C104_in_silico_perturbation_reversal" in report_md
