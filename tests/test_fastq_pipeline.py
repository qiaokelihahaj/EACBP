'''
Test Suite for EACBP Unified FASTQ-to-Biology Pipeline.
Verifies FASTQQuantificationCapability, non-spatial Cell-Cell Communication (CCI),
and full 6-plane autonomous scientific OS execution from raw FASTQ to final biological claims and report.
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
from eacbp.schemas.task import TaskContract, TaskStatus
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.quantification import FASTQQuantificationCapability
from eacbp.capabilities.spatial.cci import CellCellCommunicationCapability
from eacbp.capabilities.sc_data import SCData
from eacbp.report.markdown_report import ScientificReportGenerator


def test_fastq_quantification_capability(tmp_path):
    storage_dir = tmp_path / "artifacts_quant"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    capability = FASTQQuantificationCapability(implementation_id="sc_quant_v1")

    study_id = "test_quant_study_001"
    fastq_manifest = {
        "samples": {
            "P12K8_con": {"R1": "fake_con_R1.fq.gz", "R2": "fake_con_R2.fq.gz"},
            "P12K8_cKO": {"R1": "fake_cko_R1.fq.gz", "R2": "fake_cko_R2.fq.gz"},
        },
        "chemistry": "10xv3",
        "species": "mus_musculus",
    }

    in_uri = f"fastq://{study_id}/raw_reads/v1"
    registry.register(
        uri_str=in_uri,
        payload=fastq_manifest,
        artifact_type=ArtifactType.FASTQ,
        study_id=study_id,
        created_by_task="task_init",
        operation="raw_fastq_manifest_registration",
    )

    contract = TaskContract(
        task_id="task_000_quant",
        capability="quantification",
        method="sc_quant_v1",
        input_artifacts=[in_uri],
        allowed_operations=["kb_count_alignment", "sc_quant_demultiplex", "umi_deduplication", "gene_annotation_mapping"],
        forbidden_operations=["filter_cells", "normalize", "recluster"],
        parameters={"chemistry": "10xv3", "species": "mus_musculus", "target_gene": "Kat8", "n_cells": 300, "n_genes": 80},
        expected_outputs=[f"adata://{study_id}/raw/v1"],
    )

    result = capability.execute(contract, registry)
    assert result.status == TaskStatus.SUCCESS
    assert len(result.output_artifacts) == 1
    assert result.output_artifacts[0] == f"adata://{study_id}/raw/v1"
    assert registry.exists(f"adata://{study_id}/raw/v1")

    meta, payload = registry.get(f"adata://{study_id}/raw/v1")
    assert meta.type == ArtifactType.ANNDATA
    data = SCData.from_dict(payload) if isinstance(payload, dict) else payload
    assert data.n_obs == 300
    assert data.n_vars == 80
    assert "Kat8" in data.var["gene_name"].values


def test_non_spatial_cell_cell_communication(tmp_path):
    storage_dir = tmp_path / "artifacts_cci"
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    capability = CellCellCommunicationCapability()

    study_id = "test_cci_study_001"
    # Create single-cell data without spatial coordinates
    data = SCData.create_synthetic_ad_study(n_cells=200, n_genes=50, random_seed=42)
    in_uri = f"adata://{study_id}/annotated/v4"
    registry.register(
        uri_str=in_uri,
        payload=data.to_dict(),
        artifact_type=ArtifactType.ANNDATA,
        study_id=study_id,
        created_by_task="task_clustering",
        operation="cell_clustering_and_annotation",
    )

    contract = TaskContract(
        task_id="task_010_cell_cell_communication",
        capability="cell_cell_communication",
        method="cci_ligand_receptor_v1",
        input_artifacts=[in_uri],
        allowed_operations=["load_lr_database", "calculate_spatial_contact_density", "compute_spatial_cci_score", "run_spatial_permutation_test", "match_ligand_receptor_pairs", "score_spatial_interactions", "permutation_testing", "evaluate_cell_cell_communication", "ligand_receptor_cci"],
        forbidden_operations=["filter_cells", "recluster", "normalize"],
        parameters={"fdr_threshold": 0.05, "n_permutations": 50, "random_seed": 42},
        expected_outputs=[f"table://{study_id}/spatial_cci/v1"],
    )

    result = capability.execute(contract, registry)
    assert result.status == TaskStatus.SUCCESS
    assert registry.exists(f"table://{study_id}/spatial_cci/v1")

    meta, cci_df = registry.get(f"table://{study_id}/spatial_cci/v1")
    assert meta.type == ArtifactType.TABLE
    assert isinstance(cci_df, pd.DataFrame)
    assert not cci_df.empty
    assert "sender_cell_type" in cci_df.columns
    assert "receiver_cell_type" in cci_df.columns
    assert "ligand" in cci_df.columns
    assert "receptor" in cci_df.columns


def test_full_fastq_to_biology_autonomous_pipeline(tmp_path):
    storage_dir = tmp_path / "artifacts_full_fastq_study"
    reports_dir = tmp_path / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    prompt = (
        "输入小鼠P12组织Kat8 cKO双端FASTQ测序数据（CleanData），"
        "执行全自动单细胞定量、质控过滤、多谱系聚类、发育伪时序轨迹分析、"
        "细胞间配体受体通讯（CCI）、Kat8靶标下调与细胞周期阻滞差异表达、"
        "In silico 扰动模拟及多源知识库证据比对，产出审计完整的科研学术报告。"
    )

    study_id = "Kat8_FASTQ_To_Biology_Study_001"
    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id=study_id)
    assert manifest.data.has_raw_fastq is True
    assert manifest.data.raw_artifact_uri == f"fastq://{study_id}/raw_reads/v1"
    assert "FASTQ" in manifest.data.modalities
    assert "communication" in manifest.data.modalities

    # Register initial FASTQ reads artifact
    fastq_manifest = {
        "samples": {
            "P12K8_con": {"R1": "/data/P12K8_con_R1.fq.gz", "R2": "/data/P12K8_con_R2.fq.gz"},
            "P12K8_cKO": {"R1": "/data/P12K8_cKO_R1.fq.gz", "R2": "/data/P12K8_cKO_R2.fq.gz"},
        },
        "chemistry": "10xv3",
        "species": "mus_musculus",
    }
    registry.register(
        uri_str=manifest.data.raw_artifact_uri,
        payload=fastq_manifest,
        artifact_type=ArtifactType.FASTQ,
        study_id=study_id,
        created_by_task="task_000_ingest",
        operation="raw_fastq_reads_ingestion",
    )

    # Run Autonomous Orchestration
    summary = orchestrator.run_study(
        manifest,
        current_state={
            "target_gene": "Kat8",
            "target_genes": ["Kat8", "Kansl1", "Cdkn1a", "Bax", "Cdk1", "Top2a"],
            "include_knowledge": True,
            "include_perturbation": True,
            "include_cci": True,
            "n_cells": 300,
            "n_genes": 100,
        }
    )

    assert summary["tasks_executed"] >= 11
    assert summary["artifacts_created"] >= 14
    assert summary["claims_count"] >= 4

    # Verify all audits passed
    for report in orchestrator.audit_reports:
        assert report.overall_passed is True, f"Task {report.target_task_id} failed audit: {[c.message for c in report.checks if not c.passed]}"

    # Verify claims synthesized
    claim_ids = list(orchestrator.evidence_graph.claim_nodes.keys())
    assert "C101_microglia_state_transition" in claim_ids
    assert "C102_dam_marker_expression" in claim_ids
    assert "C103_knowledge_pathway_convergence" in claim_ids

    # Generate and verify markdown report
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()
    assert "# Scientific Study Report" in report_md
    assert "Kat8" in report_md
    assert "Evidence-Grounded Scientific Claims" in report_md
