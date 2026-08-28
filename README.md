# EACBP: Evidence-aware Agentic Computational Biology Platform

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://www.python.org/)
[![Tests](https://img.shields.io/badge/Tests-92%20passed%20(100%25)-brightgreen.svg)]()

> **Core Philosophy**: EACBP is not a monolithic "single-cell AutoGPT", but an **Evidence-aware Scientific Workflow Operating System**. It grounds computational biology analysis in auditable task DAGs, immutable data artifacts, independent statistical/biological auditing, plug-and-play agent adapters with strict contract guardrails, multi-source knowledge verification, in silico perturbation simulation, and a multi-hop traceable scientific claim engine.

$$
\boxed{
\text{Scientific Question}
\rightarrow
\text{Auditable Task DAG}
\rightarrow
\text{Capability Invocation}
\rightarrow
\text{Data Artifacts}
\rightarrow
\text{Statistical/Biological Evidence}
\rightarrow
\text{Traceable Claims}
}
$$

---

## 🏛️ System Architecture: 6 Planes, 2 DAGs, 1 Protocol

```text
                         User / API / UI
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                 Scientific Orchestrator                     │
│                                                             │
│ Intent Parser → Study Planner → Dynamic DAG Manager → Router│
│                       │                │                    │
│                       ▼                ▼                    │
│                Scientific Policy    Capability Registry     │
└───────────────────────┬─────────────────────────────────────┘
                        │
        ┌───────────────┼─────────────────────┐
        ▼               ▼                     ▼
 Knowledge Plane   Compute Plane       Model/Simulation Plane
        │               │                     │
 Literature       scRNA / scATAC        Genetic Perturbation (CRISPR)
 GO / Reactome    Spatial Transcript.   Compound Response (CMAP)
 Dual-Mode KG     Trajectory / GRN      Agent Adapters (SpaCell, ChatCell, GeneAgent)
        │               │                     │
        └───────────────┼─────────────────────┘
                        ▼
                Scientific Auditor (Author != Reviewer)
         ┌──────────────┼──────────────┐
         ▼              ▼              ▼
   Computational   Statistical      Biological
   Validator       Validator        Validator
         └──────────────┼──────────────┘
                        ▼
                  Claim Engine (4-Tier Language Protocol)
                        ▼
            Traceable Report / Notebook / Manuscript

=============================================================
             Data / State / Provenance Plane
=============================================================
 AnnData | SpatialData | Tables | Figures | Lineage | Artifacts
 Immutable Versioning (v1->v2->v3) | Strict Task Contracts
=============================================================
```

---

## 🔒 The Four Architectural Invariants

1. **Raw data immutable**: Input datasets are read-only and strictly immutable.
2. **Every computation creates an Artifact**: Every operation produces a versioned, hashed, reproducible artifact with exact lineage (`adata://<study>/<name>/v<N>`). In-place overwriting is prohibited.
3. **Every scientific claim requires Evidence**: Assertions cannot be generated from raw LLM hallucinations. Claims are backed by linked evidence nodes with explicit confidence and 4-tier language calibration.
4. **Executor cannot certify itself**: The Scientific Auditor operates independently from execution to eliminate self-confirmation bias (*Author $\neq$ Reviewer*).

---

## 🔬 Calibrated 4-Tier Scientific Language Protocol

EACBP strictly enforces epistemic calibration in all generated conclusions:

| Tier | Language Level | Example Statement | Evidence Requirement |
| :--- | :--- | :--- | :--- |
| **Level 1** | **Observation** | *"APOE expression increased along pseudotime."* | Descriptive metrics |
| **Level 2** | **Statistical inference** | *"APOE showed significant association with pseudotime ($p < 10^{-4}$)."* | Validated statistical test (FDR controlled) |
| **Level 3** | **Supported interpretation** | *"This pattern is consistent with an AD-associated microglial state."* | Statistical evidence + Literature/Biological prior |
| **Level 4** | **Hypothesis** | *"APOE may participate in driving this microglial transition."* | Prospective hypothesis (observational / in silico perturbation) |

---

## 📁 Repository Structure

```text
EACBP/
├── pyproject.toml                     # Project configuration & test paths
├── README.md                          # Platform documentation
├── eacbp/
│   ├── schemas/                       # Unified Schema Protocol
│   │   ├── study.py                   # StudyManifest, BiologicalDesign, ExperimentalDesign
│   │   ├── task.py                    # TaskContract, TaskResult, TaskStatus
│   │   ├── artifact.py                # ArtifactMetadata, ArtifactRef, LineageNode
│   │   └── evidence.py                # EvidenceNode, ClaimNode, LanguageTier, ConfidenceScore
│   ├── artifact/                      # Data, State & Provenance Plane
│   │   ├── uri.py                     # Canonical URI parser (adata://, table://, fig://)
│   │   ├── storage.py                 # File-system storage backend & SHA-256 validation
│   │   ├── lineage.py                 # Lineage Graph & branch comparison (v4a vs v4b)
│   │   └── registry.py                # ArtifactRegistry interface
│   ├── capabilities/                  # Capability Registry & Compute Plane
│   │   ├── base.py                    # BaseCapability, ImplementationType
│   │   ├── registry.py                # Global CapabilityRegistry
│   │   ├── side_effect.py             # SideEffectValidator (deterministic state hashing)
│   │   ├── sc_data.py                 # SCData container (AnnData & Spatial-compatible)
│   │   ├── qc.py                      # Dataset audit & QC filtering
│   │   ├── normalization.py           # Library size scaling, log1p & HVG selection
│   │   ├── integration.py             # Harmony & No-correction batch integration
│   │   ├── clustering.py              # PCA, KNN, Leiden clustering & marker annotation
│   │   ├── subset.py                  # Subpopulation extraction (e.g. Microglia)
│   │   ├── deg.py                     # Pseudobulk Wald/t-test vs Cell-level test
│   │   ├── trajectory.py              # PAGA/DPT pseudotime & stability analysis
│   │   ├── spatial/                   # Spatial Single-Cell Analytics
│   │   │   ├── domain.py              # Spatial graph smoothing & domain clustering
│   │   │   ├── deg.py                 # Moran's I & Geary's C spatial autocorrelation
│   │   │   └── cci.py                 # Contact-density & proximity cell-cell interaction
│   │   └── perturbation/              # In Silico Perturbation Simulation
│   │       ├── genetic.py             # CRISPR knockout / overexpression via GRN propagation
│   │       └── compound.py            # Drug response & CMAP discordance state shifting
│   ├── adapters/                      # External Agent Adapter Plane
│   │   ├── base.py                    # BaseAgentAdapter with deterministic contract guardrails
│   │   ├── spacell_adapter.py         # SpaCellAgent adapter
│   │   ├── chatcell_adapter.py        # ChatCell dialogue & state query adapter
│   │   └── gene_adapter.py            # GeneAgent pathway & function adapter
│   ├── knowledge/                     # Multi-Source Knowledge Engine Plane
│   │   ├── literature.py              # PubMed & bioRxiv retriever
│   │   ├── databases.py               # GO ORA, Reactome & NCBI Gene retriever
│   │   └── engine.py                  # Discovery Mode vs Prior-Guided Mode controller
│   ├── auditor/                       # Independent Scientific Auditor Plane
│   │   ├── base.py                    # BaseAuditor, ValidationReport, ValidationCheck
│   │   ├── computational.py           # Matrix dimensions, non-finite values (NaN/Inf)
│   │   ├── statistical.py             # Pseudoreplication, FDR, Moran's I, Perturbation bounds
│   │   └── biological.py              # Marker coherence & disease consistency
│   ├── evidence/                      # Evidence & Claim Plane
│   │   ├── graph.py                   # EvidenceDAG (5-Pillar Evidence Graph)
│   │   ├── confidence.py              # 3D Confidence scoring with contradiction penalty
│   │   ├── language.py                # 4-Tier Language protocol enforcer
│   │   └── claim.py                   # ClaimEngine synthesis
│   ├── orchestrator/                  # Scientific Orchestrator Plane
│   │   ├── intent.py                  # Natural language intent parser
│   │   ├── policy.py                  # Scientific policies & Stop Rules
│   │   ├── router.py                  # 3-Tier Router (Hard constraints -> Policy -> LLM)
│   │   ├── dag.py                     # Dynamic Computational DAG planner
│   │   └── loop.py                    # Full execution orchestration loop
│   └── report/                        # Traceability & Reporting Plane
│       ├── provenance.py              # Multi-hop sentence-to-data provenance resolver
│       └── markdown_report.py         # 4-tier structured scientific report generator
└── tests/                             # 13 test suites (92 tests, 100% pass rate)
    ├── test_schemas.py                # Schema serialization tests
    ├── test_artifacts.py              # Immutability, lineage & branch tests
    ├── test_capabilities.py           # Capability execution tests
    ├── test_spatial.py                # Spatial domain, Moran's I, DEG & CCI tests
    ├── test_adapters.py               # External agent adapter execution tests
    ├── test_adversarial_guardrails.py # Interception of rogue mutations & tampering
    ├── test_adversarial_stress.py     # Extreme mathematical stress & boundary audits
    ├── test_knowledge.py              # PubMed, bioRxiv, GO/Reactome ORA & Dual-mode tests
    ├── test_perturbation.py           # CRISPR KO, GRN propagation & drug response tests
    ├── test_auditors.py               # Computational & statistical audits
    ├── test_evidence_graph.py         # Confidence calculator & language enforcer tests
    ├── test_orchestrator.py           # Dynamic DAG planning & routing tests
    └── test_end_to_end_study.py       # Full scRNA + Spatial + Perturbation simulated studies
```

---

## 🚀 Quickstart & Example

```python
from eacbp.orchestrator.intent import IntentParser
from eacbp.orchestrator.loop import ScientificOrchestrator
from eacbp.artifact.registry import ArtifactRegistry
from eacbp.capabilities.sc_data import SCData
from eacbp.schemas.artifact import ArtifactType
from eacbp.report.markdown_report import ScientificReportGenerator

# 1. Initialize Registry & Orchestrator
registry = ArtifactRegistry(storage_dir=".artifacts")
orchestrator = ScientificOrchestrator(artifact_registry=registry)

# 2. Parse User Intent into a Structured Study Manifest
prompt = "Investigate spatial DAM localization and in silico Trem2 knockout in Alzheimer's mouse brain."
manifest = IntentParser.parse_prompt_to_manifest(prompt, study_id="AD_spatial_001")

# 3. Ingest Raw Dataset (Single-cell + Spatial coordinates)
raw_data = SCData.create_synthetic_ad_study(n_cells=1200, n_genes=500, n_ad_mice=6, n_ctrl_mice=6, has_spatial=True)
raw_uri = "adata://AD_spatial_001/raw/v1"
manifest.data.raw_artifact_uri = raw_uri
registry.register(
    uri_str=raw_uri,
    payload=raw_data.to_dict(),
    artifact_type=ArtifactType.SPATIAL_DATA,
    study_id="AD_spatial_001",
    created_by_task="task_000_ingest",
    operation="raw_data_ingest",
)

# 4. Execute Full Workflow Loop
study_results = orchestrator.run_study(manifest)
print(f"Executed {study_results['tasks_executed']} tasks, generated {study_results['claims_count']} verified claims.")

# 5. Generate Evidence-Grounded Scientific Report
report_gen = ScientificReportGenerator(
    manifest=manifest,
    evidence_graph=orchestrator.evidence_graph,
    artifact_registry=registry,
    task_history=orchestrator.task_history,
)
report_markdown = report_gen.generate_markdown()
print(report_markdown)
```

---

## 🧪 Running Tests

```bash
# Run pytest test suite across all 92 unit and integration tests
.venv/Scripts/pytest -v --basetemp=.pytest_temp
```
