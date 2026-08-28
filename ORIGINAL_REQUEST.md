# Original User Request

## 2026-08-28T06:48:28Z

<USER_REQUEST>
Extend the Evidence-aware Agentic Computational Biology Platform (EACBP) to its full Scientific Workflow OS form (Phases V3–V5), implementing Spatial single-cell analytics, External Agent Adapters with strict contract guardrails, a multi-source Knowledge Engine, and in silico Perturbation Simulation.

Working directory: `c:/Users/qiaok/Desktop/EACBP`
Integrity mode: `development`

## Reference Architecture
- Architectural Protocol: 6 Planes (Orchestration, Knowledge, Compute, Model/Simulation, Auditor, Provenance), 2 DAGs (Computational Task DAG, Evidence DAG), 1 Unified Schema Protocol (StudyManifest, TaskContract, ArtifactMetadata, Evidence/Claim).
- Four Invariants: Raw data immutable, every computation produces hashed/versioned Artifacts, every claim requires multi-tier calibrated Evidence, execution cannot self-certify.

## Requirements

### R1. Spatial Single-Cell Analytics Plane
Implement native spatial computational capabilities under `eacbp/capabilities/spatial/` or integrated modules:
- Spatial domain identification and neighborhood microenvironment analysis.
- Spatial autocorrelation (e.g., Moran's I / Geary's C) and spatial differential expression (Spatial DEG).
- Cell-cell communication (CCI) and ligand-receptor interaction analysis.
- Output versioned spatial artifacts (`adata://<study>/spatial/vN`) with full lineage tracking.

### R2. External Agent Adapter Plane with Contract Guardrails
Implement adapter layer under `eacbp/adapters/` enabling plug-and-play third-party biological agents:
- `BaseAgentAdapter` interface integrating with `TaskContract` and `SideEffectValidator`.
- Adapters for specialized biological agents:
  - `SpaCellAgentAdapter`: Spatial domain & cellular neighborhood analysis.
  - `ChatCellAdapter`: Cellular dialogue and state prediction.
  - `GeneAgentAdapter`: Biological gene function and pathway reasoning.
- Strict side-effect validation: Enforce contract constraints, automatically intercepting and failing any agent attempting unauthorized in-place mutations, uncontracted reclustering, or silent cell filtering.

### R3. Multi-Source Knowledge Engine & Verification Loop
Implement Knowledge Plane under `eacbp/knowledge/`:
- Literature retriever interface: Querying PubMed / bioRxiv metadata and citations.
- Biological Database & Knowledge Graph retriever: Querying GO terms, Reactome pathways, NCBI Gene, and entity relationships.
- Support dual operational modes:
  - **Discovery Mode**: Provides disease background and technical constraints without biasing feature selection.
  - **Prior-Guided Mode**: Evaluates user-specified hypotheses (e.g., DAM TREM2-APOE axis) with mandatory `prior-guided` epistemic tagging in reports.

### R4. In Silico Perturbation & Simulation Plane
Implement Model/Simulation Plane under `eacbp/capabilities/perturbation/`:
- In silico genetic perturbation simulation (e.g., CRISPR knockout / overexpression response).
- In silico drug response & counterfactual state transition modeling.
- Generate calibrated in silico evidence nodes (`EvidenceType.PERTURBATION`) updating causal confidence metrics in the Evidence DAG.

### R5. Orchestrator Integration & End-to-End Scientific OS Validation
Integrate all capabilities into `ScientificOrchestrator`:
- Update `ComputationalDAGPlanner` and `CapabilityRouter` to dynamically route spatial, adapter, knowledge, and perturbation tasks based on `StudyManifest` data specifications.
- Update `EvidenceGraph` and `ClaimEngine` to synthesize comprehensive multimodal claims integrating spatial localization ($++++$), pseudobulk DEG ($++++$), trajectory ($+++$), literature ($++++$), and perturbation ($+$).
- Update `ScientificReportGenerator` to render complete multi-evidence provenance reports.

## Acceptance Criteria

### Spatial Analysis & Capabilities
- [ ] Spatial domain, spatial DEG, and neighborhood analysis execute and produce immutable, hashed artifacts with verified lineage.
- [ ] Computational and statistical auditors evaluate spatial metrics without regressions.

### Agent Adapters & Guardrails
- [ ] `SpaCellAgentAdapter`, `ChatCellAdapter`, and `GeneAgentAdapter` run within `TaskContract` bounds.
- [ ] Any agent violating forbidden operations (e.g., unauthorized reclustering) is blocked with a `POLICY_VIOLATION` task status.

### Knowledge Engine & Dual Modes
- [ ] Literature and Biological DB retrieval modules provide structured evidence nodes.
- [ ] Prior-guided analysis is visibly flagged in claim metadata and report summaries.

### Perturbation Simulation
- [ ] In silico perturbation capability models gene knockouts/state shifts and creates versioned simulation artifacts.
- [ ] Evidence nodes reflect causal confidence updates without violating observational language rules.

### Test Suite & End-to-End Workflow
- [ ] Comprehensive unit and integration test suite in `tests/` passes with 100% success (`pytest -v --basetemp=.pytest_temp`).
- [ ] Full end-to-end simulated study executing scRNA + Spatial + Agent Adapter + Knowledge + Perturbation succeeds and produces a traceable markdown report.
</USER_REQUEST>
