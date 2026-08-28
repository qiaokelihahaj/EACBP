'''
EACBP Autonomous Biological Study Runner: Kat8 cKO Single-Cell Analysis.
Executes 6-plane OS: Compute, Auditor, Knowledge, Model/Simulation, Evidence, and Provenance.
Outputs comprehensive publication-grade markdown scientific study report.
'''

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Ensure EACBP root in PYTHONPATH
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

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
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.report.markdown_report import ScientificReportGenerator


def run_kat8_cKO_autonomous_study(output_dir: Path, data_path: Path = None, n_cells: int = 1500):
    output_dir.mkdir(parents=True, exist_ok=True)
    storage_dir = output_dir / "artifacts_kat8_cKO"
    reports_dir = output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("EACBP Autonomous Scientific OS: Initiating Kat8 cKO Single-Cell Study")
    print("=" * 80)

    # 1. Initialize Artifact Registry and Scientific Orchestrator
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    study_id = "Kat8_P12_cKO_SingleCell_Study_001"

    # 2. Intent Parsing & Study Manifest Construction
    prompt = (
        "分析小鼠P12组织Kat8条件性敲除(cKO)单细胞转录组数据，"
        "探索Kat8催化的H4K16ac表观调控缺失引起的细胞亚群重塑、核心靶基因差异表达、"
        "细胞周期与凋亡应激激活，以及发育分化轨迹阻滞机制。"
    )
    print(f"\n[Plane 1: Orchestration] Parsing Research Intent...")
    print(f"User Query: {prompt}")

    manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id=study_id)
    manifest.title = "Single-Cell Epigenetic Transcriptomic Dissection of Kat8 (Mof) cKO at Postnatal Stage P12"

    print(f"Study ID: {manifest.study_id}")
    print(f"Species: {manifest.biological_design.species}")
    print(f"Disease/Intervention: {manifest.biological_design.disease}")
    print(f"Conditions: {manifest.biological_design.conditions}")
    print(f"Prior-guided Hypotheses: {manifest.hypotheses.user_provided}")

    # 3. Data Ingestion (from real h5ad if available or high-res structured profile)
    print(f"\n[Plane 2: Compute] Ingesting Single-Cell Data into Content-Addressed Storage...")
    if data_path and data_path.exists() and data_path.suffix == ".h5ad":
        print(f"Loading data from file: {data_path}")
        raw_data = SCData.from_h5ad(str(data_path), max_cells=n_cells)
    else:
        print(f"Synthesizing high-fidelity P12 Kat8 cKO single-cell cohort ({n_cells} cells, 4 cKO vs 4 con mice)...")
        raw_data = SCData.create_synthetic_kat8_study(
            n_cells=n_cells,
            n_genes=400,
            n_cko_mice=4,
            n_con_mice=4,
            random_seed=42,
        )

    raw_uri = f"adata://{study_id}/raw/v1"
    manifest.data.raw_artifact_uri = raw_uri

    if not registry.exists(raw_uri):
        reg_meta = registry.register(
            uri_str=raw_uri,
            payload=raw_data.to_dict(),
            artifact_type=ArtifactType.ANNDATA,
            study_id=study_id,
            created_by_task="task_000_ingest",
            operation="raw_kat8_data_ingest",
        )
    else:
        reg_meta = registry.get_metadata(raw_uri)
    print(f"Raw Data Registered: {raw_uri} (SHA256: {reg_meta.sha256_hash[:16]}...)")

    # 4. Execute Full Multi-Plane Autonomous Study Workflow
    print(f"\n[Orchestrator] Executing 6-Plane Scientific Task DAG...")
    study_summary = orchestrator.run_study(
        manifest,
        current_state={
            "target_gene": "Kat8",
            "target_genes": ["Kat8", "Kansl1", "Msl1", "Cdkn1a", "Bax", "Cdk1", "Top2a", "Trp53", "H2ax"],
            "include_knowledge": True,
            "include_perturbation": True,
            "run_compound_perturbation": False,
        }
    )

    print(f"\n[Summary] Autonomous Workflow Execution Complete:")
    print(f"  - Tasks Executed: {study_summary['tasks_executed']}")
    print(f"  - Artifacts Created: {study_summary['artifacts_created']}")
    print(f"  - Evidence Nodes Synthesized: {study_summary['evidence_nodes_count']}")
    print(f"  - Multimodal Claims Formulated: {study_summary['claims_count']}")

    # 5. Review Independent Scientific Audits
    print(f"\n[Plane 3: Auditor] Scientific Audit Verification:")
    all_passed = True
    for report in orchestrator.audit_reports:
        status_str = "[PASSED]" if report.overall_passed else "[FAILED]"
        print(f"  {status_str} Task: {report.target_task_id:<30} Target: {report.target_artifact_uri}")
        if not report.overall_passed:
            all_passed = False
            for c in report.checks:
                if not c.passed:
                    print(f"      -> Warning: {c.message}")

    print(f"Overall Audit Status: {'ALL CHECKS PASSED' if all_passed else 'AUDIT WARNINGS FLAGGED'}")

    # 6. Print Evidence DAG and Calibrated Claims
    print(f"\n[Plane 4: Evidence & Claim] Synthesized Scientific Claims (4-Tier Protocol):")
    for claim_id, claim in orchestrator.evidence_graph.claim_nodes.items():
        print(f"\n  Claim ID: {claim_id}")
        print(f"  Statement: {claim.statement}")
        print(f"  Language Tier: {claim.language_tier.value}")
        print(f"  Causal Status: {claim.causal_status}")
        print(f"  Confidence: Overall={claim.confidence.overall:.3f}, Association={claim.confidence.association:.3f}, Mechanistic={claim.confidence.mechanistic:.3f}, Causal={claim.confidence.causal:.3f}")
        print(f"  Supporting Evidences: {claim.support_evidence_ids}")

    # 7. Generate Full Provenance-Tracked Markdown Scientific Report
    print(f"\n[Plane 5 & 6: Provenance & Report] Rendering Scientific Markdown Report...")
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()

    report_path = reports_dir / "Kat8_cKO_Study_Report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_md)

    print(f"Report Generated and Saved to: {report_path}")
    print("\n" + "=" * 80)
    print("EACBP Study Execution Finished Successfully!")
    print("=" * 80)
    return report_path, study_summary


if __name__ == "__main__":
    out_dir = PROJECT_ROOT / "outputs"
    run_kat8_cKO_autonomous_study(output_dir=out_dir)
