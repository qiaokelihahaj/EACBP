'''
Test Suite for Kat8 cKO Autonomous Single-Cell Study Execution in EACBP.
Validates the 6-plane OS on Kat8 conditional knockout epigenetic datasets.
'''

import pytest
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
from eacbp.schemas.evidence import LanguageTier, ClaimType, EvidenceType
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.report.markdown_report import ScientificReportGenerator


def test_kat8_intent_parsing_and_manifest():
    '''Validates that IntentParser correctly translates Kat8 cKO research queries.'''
    prompt = "分析小鼠P12脑组织Kat8条件性敲除(cKO)单细胞数据，探索Kat8调控靶基因与发育轨迹阻滞机制。"
    study_id = "Kat8_P12_cKO_001"
    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id=study_id)

    assert manifest.study_id == study_id
    assert manifest.biological_design.species == "mus_musculus"
    assert manifest.biological_design.conditions == ["cKO", "con"]
    assert "perturbation" in manifest.data.modalities
    assert manifest.analysis_policy.prior_guided_analysis is True
    assert len(manifest.hypotheses.user_provided) >= 1
    assert "Kat8" in manifest.hypotheses.user_provided[0]


def test_full_kat8_cko_study_pipeline(tmp_path):
    '''
    Executes a complete 6-plane autonomous study on Kat8 cKO single-cell data.
    Validates data ingest, QC, normalization, clustering, pseudobulk DEG, trajectory,
    knowledge retrieval, independent auditing, evidence DAG synthesis, and provenance report.
    '''
    storage_dir = tmp_path / "artifacts_kat8_study"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    study_id = "Kat8_P12_cKO_Study_001"

    # 1. Intent Parsing
    prompt = "分析小鼠P12组织Kat8条件性敲除(cKO)单细胞数据，探索Kat8靶基因、H4K16ac表观调控与细胞发育轨迹阻滞。"
    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id=study_id)

    # 2. Ingest synthetic Kat8 dataset (4 cKO vs 4 con mice, 800 cells)
    raw_data = SCData.create_synthetic_kat8_study(
        n_cells=800,
        n_genes=300,
        n_cko_mice=4,
        n_con_mice=4,
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
        operation="raw_kat8_data_ingest",
    )

    # 3. Run Autonomous Study
    study_summary = orchestrator.run_study(
        manifest,
        current_state={
            "target_gene": "Kat8",
            "target_genes": ["Kat8", "Kansl1", "Cdkn1a", "Bax", "Cdk1", "Top2a"],
            "include_knowledge": True,
            "include_perturbation": True,
        }
    )

    # 4. Verify Execution
    assert study_summary["study_id"] == study_id
    assert study_summary["tasks_executed"] >= 8
    assert study_summary["artifacts_created"] >= 8
    assert study_summary["evidence_nodes_count"] >= 3
    assert study_summary["claims_count"] >= 2

    # 5. Check Invariants
    assert registry.exists(raw_uri)
    assert registry.exists(f"adata://{study_id}/qc/v1")
    assert registry.exists(f"adata://{study_id}/normalized/v2")
    assert registry.exists(f"table://{study_id}/deg_results/v1")

    # 6. Verify Audits
    for report in orchestrator.audit_reports:
        assert report.overall_passed is True

    # 7. Check DEG Results contain Kat8
    meta, deg_payload = registry.get(f"table://{study_id}/deg_results/v1")
    deg_df = deg_payload if isinstance(deg_payload, pd.DataFrame) else pd.DataFrame(deg_payload)
    assert "Kat8" in deg_df["gene"].values
    kat8_row = deg_df[deg_df["gene"] == "Kat8"].iloc[0]
    assert kat8_row["log2_fold_change"] < 0 or kat8_row["fdr_q_value"] < 0.05

    # 8. Generate Traceable Markdown Report
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()

    assert "Scientific Study Report" in report_md
    assert study_id in report_md
    assert "Kat8" in report_md or "KAT8" in report_md
    assert "Lineage Graph" in report_md
