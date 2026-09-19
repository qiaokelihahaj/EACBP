'''
EACBP Autonomous FASTQ-to-Biology Scientific OS Pipeline Runner.
Executes end-to-end single-cell transcriptomics study directly from raw FASTQ reads:
  1. FASTQ Quantification (kb-python / STARsolo / native SC quantifier)
  2. Quality Control & Low-quality Filtering (QCCapability)
  3. Log1p Normalization & HVG Selection (NormalizationCapability)
  4. Batch Effect Correction & Integration (BatchIntegrationCapability via Harmony)
  5. Multi-lineage Unsupervised Clustering & Cell Type Annotation (ClusteringCapability)
  6. Subpopulation Subsetting (SubsetCapability)
  7. Cell Subpopulation Differential Abundance (DifferentialAbundanceCapability)
  8. Donor-level Pseudobulk Differential Expression (DifferentialExpressionCapability)
  9. Developmental Trajectory & Pseudotime Inference (TrajectoryCapability via PAGA/DPT)
 10. Ligand-Receptor Cell-Cell Communication (CellCellCommunicationCapability)
 11. Multi-source Prior Knowledge Engine & Literature Grounding (KnowledgeEngine)
 12. In Silico Genetic Perturbation Simulation (GeneticPerturbationCapability)
 13. Independent Scientific Auditing across 6 Planes (ScientificAuditor)
 14. Multimodal Evidence DAG Synthesis & Calibrated Scientific Claims (ClaimEngine)
 15. Publication-grade Markdown Scientific Study Report (ScientificReportGenerator)
'''

import sys
import os
import argparse
import json
import re
import uuid
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from eacbp.schemas.artifact import ArtifactType
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.report.markdown_report import ScientificReportGenerator


_SAFE_COMPONENT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def _validate_output_component(value: str, label: str) -> str:
    value = str(value or "")
    if value in {".", ".."} or not _SAFE_COMPONENT.fullmatch(value):
        raise ValueError(f"{label} must be a single safe path component, got {value!r}")
    return value


def run_fastq_to_biology(
    fastq_dir: Path = None,
    output_dir: Path = None,
    study_id: str = "Kat8_P12_cKO_FASTQ_To_Biology_001",
    target_gene: str = "Kat8",
    species: str = "mus_musculus",
    chemistry: str = "10xv3",
    n_cells: int = 1500,
    n_genes: int = 400,
    mode: str = "real",
    sample_manifest: dict = None,
    run_id: str = None,
    index_path: Path = None,
    t2g_path: Path = None,
    threads: int = 16,
    quant_tool: str = "kb_python_v1",
    star_bin: str = "STAR",
    genome_dir: Path = None,
    whitelist_path: Path = None,
    gtf_path: Path = None,
    method_profile: str = "baseline",
    root_cell_id: str = None,
):
    mode = str(mode).strip().lower()
    if mode not in {"real", "demo"}:
        raise ValueError("mode must be 'real' or 'demo'")
    study_id = _validate_output_component(study_id, "study_id")
    output_dir = output_dir or PROJECT_ROOT / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    run_id = str(run_id or (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:10]))
    run_id = _validate_output_component(run_id, "run_id")
    # Each invocation has an isolated output directory.  Existing studies and
    # failed runs remain available for inspection and recovery.
    run_dir = output_dir / "runs" / study_id / run_id
    storage_dir = run_dir / "artifacts_fastq_to_biology"
    storage_dir.mkdir(parents=True, exist_ok=False)
    reports_dir = run_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 85)
    print(" EACBP Autonomous Scientific OS: Unified FASTQ-to-Biology Pipeline")
    print("=" * 85)

    # 1. Initialize Registry & Orchestrator
    registry = ArtifactRegistry(storage_dir=str(storage_dir))
    orchestrator = ScientificOrchestrator(artifact_registry=registry)

    # 2. Parse Intent & Construct Study Manifest
    prompt = (
        f"输入{species}组织{target_gene} cKO双端FASTQ测序数据（CleanData），"
        "执行全自动单细胞定量、质控过滤、多谱系聚类、发育伪时序轨迹分析、"
        "细胞间配体受体通讯（CCI）、核心靶基因差异表达、"
        "In silico 扰动模拟及多源知识库证据比对，产出严格审计的科研学术报告。"
    )
    print(f"\n[Plane 1: Orchestration] Parsing Scientific Intent from Query...")
    manifest = IntentParser.parse_prompt_to_manifest(
        prompt,
        study_id=study_id,
        sample_manifest=sample_manifest,
    )
    manifest.title = f"End-to-End Single-Cell Epigenetic Dissection of {target_gene} cKO from Raw Sequencing Reads"

    print(f"  Study ID: {manifest.study_id}")
    print(f"  Modalities: {manifest.data.modalities}")
    print(f"  Raw Artifact URI: {manifest.data.raw_artifact_uri}")
    print(f"  Prior Hypotheses: {manifest.hypotheses.user_provided}")

    # 3. Register Raw FASTQ Artifact
    samples = {}
    explicit_samples = sample_manifest.get("samples", sample_manifest) if isinstance(sample_manifest, dict) else None
    if isinstance(explicit_samples, dict) and explicit_samples:
        # A supplied manifest is authoritative.  In particular, preserve all
        # lane paths and metadata instead of guessing from file names.
        samples = {str(name): dict(spec) for name, spec in explicit_samples.items()}
    elif fastq_dir and fastq_dir.exists():
        print(f"\n[Plane 2: Compute] Scanning raw FASTQ files recursively in: {fastq_dir}")
        all_fqs = sorted(set(list(fastq_dir.rglob("*.fastq.gz")) + list(fastq_dir.rglob("*.fq.gz"))))
        print(f"  Discovered {len(all_fqs)} total FASTQ sequencing files.")
        
        # Group by sample (e.g. P12K8_con, P12K8_cKO)
        sample_names = set()
        for fq in all_fqs:
            for part in fq.parts:
                if any(k in part.lower() for k in ["con", "cko", "k8", "p12", "sample"]):
                    sample_names.add(part)
        
        for s_name in sorted(sample_names):
            s_fqs = [f for f in all_fqs if s_name in f.parts and "cdna" in str(f).lower()]
            if not s_fqs:
                s_fqs = [f for f in all_fqs if s_name in f.parts]
            r1 = sorted(str(f) for f in s_fqs if "_R1" in f.name or "_1." in f.name)
            r2 = sorted(str(f) for f in s_fqs if "_R2" in f.name or "_2." in f.name)
            if r1 and r2:
                samples[s_name] = {"R1": r1, "R2": r2}
                print(f"  Matched Sample '{s_name}': {len(r1)} lane pair(s)")

    if not samples and mode == "demo":
        print(f"\n[Plane 2: Compute] Initializing FASTQ manifest for {study_id}...")
        samples = {
            "P12K8_con": {
                "R1": "demo/P12K8_con_R1.fastq.gz",
                "R2": "demo/P12K8_con_R2.fastq.gz",
                "metadata": {"donor": "demo_con_01", "condition": "con", "batch": "demo_batch_1"},
            },
            "P12K8_cKO": {
                "R1": "demo/P12K8_cKO_R1.fastq.gz",
                "R2": "demo/P12K8_cKO_R2.fastq.gz",
                "metadata": {"donor": "demo_cko_01", "condition": "cKO", "batch": "demo_batch_1"},
            },
        }

    fastq_payload = {
        "samples": samples,
        "chemistry": chemistry,
        "species": species,
        "target_gene": target_gene,
        "mode": mode,
    }

    raw_uri = manifest.data.raw_artifact_uri
    if not registry.exists(raw_uri):
        reg_meta = registry.register(
            uri_str=raw_uri,
            payload=fastq_payload,
            artifact_type=ArtifactType.FASTQ,
            study_id=study_id,
            created_by_task="task_000_ingest",
            operation="raw_fastq_reads_ingestion",
        )
        print(f"  Raw FASTQ Registered: {raw_uri} (SHA-256: {reg_meta.sha256_hash[:16]}...)")
    else:
        reg_meta = registry.get_metadata(raw_uri)
        print(f"  Raw FASTQ Retrieved: {raw_uri}")

    # 4. Execute Full 6-Plane OS Workflow
    print(f"\n[Orchestrator] Executing Full 6-Plane Scientific Task DAG...")
    summary = orchestrator.run_study(
        manifest,
        current_state={
            "target_gene": target_gene,
            "target_genes": [target_gene, "Kansl1", "Msl1", "Cdkn1a", "Bax", "Cdk1", "Top2a", "Trp53", "H2ax"],
            "include_knowledge": True,
            "include_perturbation": True,
            "include_cci": True,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "mode": mode,
            "run_id": run_id,
            "threads": int(threads),
            "quant_tool": quant_tool,
            "star_bin": star_bin,
            "genome_dir": str(genome_dir) if genome_dir else None,
            "whitelist_path": str(whitelist_path) if whitelist_path else None,
            "gtf_path": str(gtf_path) if gtf_path else None,
            "method_profile": method_profile,
            "capability_parameters": {"trajectory_inference": {"root_cell_id": root_cell_id}},
            "work_dir": str(run_dir / "kb_quant"),
            "index_path": str(index_path) if index_path else None,
            "t2g_path": str(t2g_path) if t2g_path else None,
        }
    )

    print(f"\n[Summary] Autonomous Workflow Execution Complete:")
    print(f"  - Tasks Executed: {summary['tasks_executed']}")
    print(f"  - Artifacts Created: {summary['artifacts_created']}")
    print(f"  - Evidence Nodes Synthesized: {summary['evidence_nodes_count']}")
    print(f"  - Multimodal Claims Formulated: {summary['claims_count']}")

    # 5. Review Independent Scientific Audits
    print(f"\n[Plane 3: Auditor] Scientific Audit Verification:")
    all_passed = summary.get("status") == "success"
    for report in orchestrator.audit_reports:
        status_str = "[PASSED]" if report.overall_passed else "[FAILED]"
        print(f"  {status_str} Task: {report.target_task_id:<32} Target: {report.target_artifact_uri}")
        if not report.overall_passed:
            all_passed = False
            for c in report.checks:
                if not c.passed:
                    print(f"      -> Warning: {c.message}")

    print(f"Overall Audit Status: {'ALL CHECKS PASSED' if all_passed else 'INCOMPLETE / FAILURES PRESENT'}")
    if summary.get("failures"):
        for failure in summary["failures"]:
            print(f"  Failure: {failure.get('task_id')}: {failure.get('error')}")

    # 6. Print Evidence DAG and Calibrated Claims
    print(f"\n[Plane 4: Evidence & Claim] Synthesized Scientific Claims (4-Tier Protocol):")
    for claim_id, claim in orchestrator.evidence_graph.claim_nodes.items():
        print(f"\n  Claim ID: {claim_id}")
        print(f"  Statement: {claim.statement}")
        print(f"  Language Tier: {claim.language_tier.value}")
        print(f"  Confidence: Overall={claim.confidence.overall:.3f}, Association={claim.confidence.association:.3f}, Mechanistic={claim.confidence.mechanistic:.3f}, Causal={claim.confidence.causal:.3f}")
        print(f"  Supporting Evidences ({len(claim.support_evidence_ids)}): {claim.support_evidence_ids[:5]}...")

    # 7. Render Provenance-Tracked Markdown Scientific Study Report
    print(f"\n[Plane 5 & 6: Provenance & Report] Rendering Scientific Markdown Report...")
    report_gen = ScientificReportGenerator(
        manifest=manifest,
        evidence_graph=orchestrator.evidence_graph,
        artifact_registry=registry,
        task_history=orchestrator.task_history,
        audit_reports=orchestrator.audit_reports,
    )
    report_md = report_gen.generate_markdown()

    report_path = reports_dir / f"{study_id}_Report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_md)

    print(f"Report Generated and Saved to: {report_path}")
    print("\n" + "=" * 85)
    if summary.get("status") == "success":
        print(" EACBP FASTQ-to-Biology Pipeline Finished Successfully!")
    else:
        print(" EACBP FASTQ-to-Biology Pipeline finished INCOMPLETE; see the report and failure summary.")
    print("=" * 85)
    return report_path, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="EACBP FASTQ-to-Biology Autonomous Pipeline")
    parser.add_argument("--fastq-dir", type=str, default=None, help="Path to raw FASTQ directory")
    parser.add_argument("--output-dir", type=str, default=str(PROJECT_ROOT / "outputs"), help="Output directory")
    parser.add_argument("--target-gene", type=str, default="Kat8", help="Target gene of interest")
    parser.add_argument("--mode", choices=("real", "demo"), default="real", help="Quantification mode (default: real)")
    parser.add_argument("--sample-manifest", type=str, default=None, help="JSON file with FASTQ paths and donor/condition/batch metadata")
    parser.add_argument("--run-id", type=str, default=None, help="Optional unique run identifier")
    parser.add_argument("--index-path", type=str, default=None, help="Explicit kb-python index path for mode=real")
    parser.add_argument("--t2g-path", type=str, default=None, help="Explicit transcript-to-gene path for mode=real")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--quant-tool", choices=("kb_python_v1", "starsolo_v1"), default="kb_python_v1")
    parser.add_argument("--star-bin", default="STAR")
    parser.add_argument("--genome-dir", type=Path)
    parser.add_argument("--whitelist-path", type=Path)
    parser.add_argument("--gtf-path", type=Path)
    parser.add_argument("--method-profile", choices=("baseline", "standard"), default="standard")
    parser.add_argument("--root-cell-id")
    args = parser.parse_args()

    fq_path = Path(args.fastq_dir) if args.fastq_dir else None
    manifest_data = None
    if args.sample_manifest:
        with open(args.sample_manifest, "r", encoding="utf-8") as handle:
            manifest_data = json.load(handle)
    _report_path, run_summary = run_fastq_to_biology(
        fastq_dir=fq_path,
        output_dir=Path(args.output_dir),
        target_gene=args.target_gene,
        mode=args.mode,
        sample_manifest=manifest_data,
        run_id=args.run_id,
        index_path=Path(args.index_path) if args.index_path else None,
        t2g_path=Path(args.t2g_path) if args.t2g_path else None,
        threads=args.threads,
        quant_tool=args.quant_tool,
        star_bin=args.star_bin,
        genome_dir=args.genome_dir,
        whitelist_path=args.whitelist_path,
        gtf_path=args.gtf_path,
        method_profile="baseline" if args.mode == "demo" else args.method_profile,
        root_cell_id=args.root_cell_id,
    )
    raise SystemExit(0 if run_summary.get("status") == "success" else 1)
