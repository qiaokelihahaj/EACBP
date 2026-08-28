# EACBP E2E Test Infrastructure & Test Architecture Specification (Phases V3–V5)

**Document Version**: 1.0.0  
**Target System**: Evidence-aware Agentic Computational Biology Platform (EACBP)  
**Scope**: 4-Tier Test Architecture (Category-Partition, Boundary Value Analysis, Pairwise Combinatorial, Real-World Workload) for EACBP V3–V5 Scientific Workflow OS.

---

## 1. Test Philosophy & Engineering Principles

### 1.1 Opaque-Box & Requirement-Driven Testing
The EACBP test infrastructure adheres strictly to **opaque-box, requirement-driven verification**:
- **Interface Contract Decoupling**: Tests interact exclusively with public contracts (`StudyManifest`, `TaskContract`, `TaskResult`, `ArtifactRegistry`, `SCData`, `EvidenceGraph`, `ClaimNode`, `ScientificOrchestrator`) and immutable artifacts (`adata://`, `table://`, `fig://`). No tests depend on private implementation helpers or fragile internal data layout.
- **Specification-Grounding**: Every test is mapped directly to authoritative requirements defined in `ORIGINAL_REQUEST.md`, `PROJECT.md`, and scientific algorithms in `survey_explorer_2/analysis.md`.
- **Epistemic & Mathematical Grounding**: Test assertions enforce verified mathematical laws (e.g., Moran's $I$ analytical bounds, hypergeometric hypergeometric distribution CDFs, GRN matrix propagation convergence, FDR Benjamini-Hochberg monotonicity, 4-tier language calibration).

### 1.2 The Four Invariants as Test Oracles
All test cases across all tiers continuously validate the Four Core Architectural Invariants:
1. **Raw Data Immutability**: `raw/v1` payload buffers and disk files cannot be altered in-place or deleted.
2. **Hashed & Versioned Artifacts**: Every task yields a content-addressed SHA-256 verified artifact (`adata://<study>/<stage>/vN`) registered in `ArtifactRegistry`. In-place overwrites raise `ArtifactAlreadyExistsError`.
3. **Multi-Tier Calibrated Evidence**: Claims strictly reflect quantitative 3D confidence scores (Association, Mechanistic, Causal). Observational studies strictly ban causal verbs (`causes`, `drives`, `proves`).
4. **Independent Scientific Auditing**: Execution cannot self-certify (*Author $\neq$ Reviewer*). Every task result is evaluated by independent `ComputationalValidator`, `StatisticalValidator`, and `BiologicalValidator`.

### 1.3 Test Independence, Isolation & Reproducibility
- **Filesystem Isolation**: Tests execute in ephemeral fixture directories via `pytest -v --basetemp=.pytest_temp`.
- **State Independence**: Every test instantiates isolated `ArtifactRegistry`, `ArtifactStorageBackend`, `CapabilityRegistry`, and `EvidenceGraph` instances.
- **Determinism**: All stochastic operations (random subsampling, Leiden clustering, PCA initialization, permutation tests) use fixed random seeds (`random_seed=42`).

---

## 2. 24-Feature Inventory & Multi-Tier Test Mapping

The platform capabilities are inventoried across 24 distinct functional features spanning all 6 architectural planes:

| # | Feature ID | Feature Name | Plane | Milestone | Tier 1 (Coverage) | Tier 2 (Boundary) | Tier 3 (Pairwise) | Tier 4 (E2E) |
|---|---|---|---|---|---|---|---|---|
| 1 | `F01` | Spatial Data Model & Coordinates | Compute | M1 | `T1-F01-01..05` | `T2-F01-01..05` | `T3-PAIR-01,02` | `T4-E2E-STEP-01` |
| 2 | `F02` | Spatial Domain Identification | Compute | M1 | `T1-F02-01..05` | `T2-F02-01..05` | `T3-PAIR-01` | `T4-E2E-STEP-07` |
| 3 | `F03` | Spatial Autocorrelation (Moran's I / Geary's C) | Compute | M1 | `T1-F03-01..05` | `T2-F03-01..05` | `T3-PAIR-02` | `T4-E2E-STEP-08` |
| 4 | `F04` | Spatial DEG (Spatially Variable Genes) | Compute | M1 | `T1-F04-01..05` | `T2-F04-01..05` | `T3-PAIR-02` | `T4-E2E-STEP-08` |
| 5 | `F05` | Cell-Cell Communication (Spatial CCI) | Compute | M1 | `T1-F05-01..05` | `T2-F05-01..05` | `T3-PAIR-03` | `T4-E2E-STEP-09` |
| 6 | `F06` | Spatial Artifact Lineage & Storage | Provenance | M1 | `T1-F06-01..05` | `T2-F06-01..05` | `T3-PAIR-07` | `T4-E2E-STEP-07` |
| 7 | `F07` | BaseAgentAdapter Protocol | Model/Sim | M2 | `T1-F07-01..05` | `T2-F07-01..05` | `T3-PAIR-01,04` | `T4-E2E-STEP-07` |
| 8 | `F08` | SpaCellAgentAdapter | Model/Sim | M2 | `T1-F08-01..05` | `T2-F08-01..05` | `T3-PAIR-01` | `T4-E2E-STEP-07` |
| 9 | `F09` | ChatCellAdapter | Model/Sim | M2 | `T1-F09-01..05` | `T2-F09-01..05` | `T3-PAIR-04` | `T4-E2E-STEP-15` |
| 10 | `F10` | GeneAgentAdapter | Model/Sim | M2 | `T1-F10-01..05` | `T2-F10-01..05` | `T3-PAIR-05` | `T4-E2E-STEP-13` |
| 11 | `F11` | SideEffectValidator Guardrail | Auditor | M2 | `T1-F11-01..05` | `T2-F11-01..05` | `T3-PAIR-01,07` | `T4-E2E-STEP-07` |
| 12 | `F12` | Contract-Bound Artifact Wrapping | Provenance | M2 | `T1-F12-01..05` | `T2-F12-01..05` | `T3-PAIR-07` | `T4-E2E-STEP-07` |
| 13 | `F13` | Literature Retriever (PubMed/bioRxiv) | Knowledge | M3 | `T1-F13-01..05` | `T2-F13-01..05` | `T3-PAIR-02,08` | `T4-E2E-STEP-11` |
| 14 | `F14` | Biological DB Retriever (GO/Reactome/NCBI) | Knowledge | M3 | `T1-F14-01..05` | `T2-F14-01..05` | `T3-PAIR-02,05` | `T4-E2E-STEP-12` |
| 15 | `F15` | Discovery Mode Workflow | Knowledge | M3 | `T1-F15-01..05` | `T2-F15-01..05` | `T3-PAIR-06` | `T4-E2E-STEP-01` |
| 16 | `F16` | Prior-Guided Mode & Epistemic Tagging | Knowledge | M3 | `T1-F16-01..05` | `T2-F16-01..05` | `T3-PAIR-06` | `T4-E2E-STEP-11` |
| 17 | `F17` | Knowledge Evidence Extraction | Evidence | M3 | `T1-F17-01..05` | `T2-F17-01..05` | `T3-PAIR-08` | `T4-E2E-STEP-16` |
| 18 | `F18` | Genetic Perturbation Simulation (CRISPR KO/OE)| Model/Sim | M4 | `T1-F18-01..05` | `T2-F18-01..05` | `T3-PAIR-03,05` | `T4-E2E-STEP-14` |
| 19 | `F19` | Drug Response & Counterfactual Simulation | Model/Sim | M4 | `T1-F19-01..05` | `T2-F19-01..05` | `T3-PAIR-04` | `T4-E2E-STEP-15` |
| 20 | `F20` | Perturbation Evidence Calibration | Evidence | M4 | `T1-F20-01..05` | `T2-F20-01..05` | `T3-PAIR-08` | `T4-E2E-STEP-16` |
| 21 | `F21` | Perturbation Artifact Lineage | Provenance | M4 | `T1-F21-01..05` | `T2-F21-01..05` | `T3-PAIR-07` | `T4-E2E-STEP-14` |
| 22 | `F22` | Dynamic DAG Planning & Routing | Orchestration| M5 | `T1-F22-01..05` | `T2-F22-01..05` | `T3-PAIR-07` | `T4-E2E-STEP-01` |
| 23 | `F23` | Multimodal Claim Synthesis (5-Pillar Calibration)| Evidence | M5 | `T1-F23-01..05` | `T2-F23-01..05` | `T3-PAIR-06,08` | `T4-E2E-STEP-17` |
| 24 | `F24` | Multimodal Provenance Reporting | Provenance | M5 | `T1-F24-01..05` | `T2-F24-01..05` | `T3-PAIR-06` | `T4-E2E-STEP-18` |

---

## 3. Tier 1: Feature Coverage Specifications (Isolation & Happy Path)

Tier 1 contains **120 test cases** ($\ge 5$ per feature) executing isolated happy-path and nominal parameter variations for each feature in the inventory.

### 3.1 Spatial Analytics Plane (F01 – F06)
- **F01: Spatial Data Model & Coordinates**
  - `T1-F01-01`: Verify 2D coordinate parsing into `.obsm['spatial']` ($(N, 2)$ array) with positive physical distances.
  - `T1-F01-02`: Verify 3D coordinate parsing ($(N, 3)$ array) for volumetric spatial datasets.
  - `T1-F01-03`: Verify $k$-NN spatial connectivity matrix construction in `.obsp['spatial_connectivities']`.
  - `T1-F01-04`: Verify spatial pairwise Euclidean distance matrix calculation in `.obsp['spatial_distances']`.
  - `T1-F01-05`: Verify `SCData.create_synthetic_spatial_ad_study()` generates valid spatial coordinates and plaque distance annotations.
- **F02: Spatial Domain Identification**
  - `T1-F02-01`: Verify `SpatialDomainCapability` clusters spots into plaque-adjacent vs distal microenvironments.
  - `T1-F02-02`: Verify spatially regularized PCA embedding $\mathbf{\tilde{Z}}_i = (1 - \lambda)\mathbf{Z}_i + \lambda \sum \tilde{w}_{ij}\mathbf{Z}_j$ with $\lambda = 0.3$.
  - `T1-F02-03`: Verify neighborhood cell-type composition vector generation ($\mathbf{c}_i \in \mathbb{R}^{|\mathcal{T}|}$).
  - `T1-F02-04`: Verify spatial domain cluster assignment is written to `.obs['spatial_domain']`.
  - `T1-F02-05`: Verify silhouette score calculation for spatial domains ($\ge 0.15$).
- **F03: Spatial Autocorrelation (Moran's I / Geary's C)**
  - `T1-F03-01`: Verify Moran's $I$ calculation for localized synthetic marker (*Apoe* $I > 0.40$).
  - `T1-F03-02`: Verify Moran's $I$ calculation for randomly shuffled expression yields $I \approx -\frac{1}{N-1}$.
  - `T1-F03-03`: Verify analytical variance computation $\text{Var}(I)$ under randomization null.
  - `T1-F03-04`: Verify Geary's $C$ metric calculation ($C < 0.60$ for spatially clustered genes).
  - `T1-F03-05`: Verify Benjamini-Hochberg FDR calculation across $G$ tested genes.
- **F04: Spatial DEG (Spatially Variable Genes)**
  - `T1-F04-01`: Verify GLM detection of plaque-distance-dependent differential expression.
  - `T1-F04-02`: Verify spatial domain DEG table output contains `gene_symbol`, `log2_fc`, `morans_i`, `fdr_q_val`.
  - `T1-F04-03`: Verify filtering of SVGs using dual thresholds: $I_g \ge 0.20$ and $q_g < 0.05$.
  - `T1-F04-04`: Verify spatial DEG execution generates output artifact `table://<study>/spatial_deg/v1`.
  - `T1-F04-05`: Verify statistical auditor sign-off on spatial DEG table format.
- **F05: Cell-Cell Communication (Spatial CCI)**
  - `T1-F05-01`: Verify ligand-receptor pair matching from curated database (*Apoe-Trem2*, *App-Cd74*).
  - `T1-F05-02`: Verify spatial contact density calculation $W_{spatial}(A, B)$ between cell types $A$ and $B$.
  - `T1-F05-03`: Verify proximity-weighted score $S_{spatial} = S_{raw} \times W_{spatial}$.
  - `T1-F05-04`: Verify coordinate permutation testing ($B = 1000$ shuffles) empirical $p$-value computation.
  - `T1-F05-05`: Verify output artifact `table://<study>/spatial_cci/v1` containing interaction scores and FDR values.
- **F06: Spatial Artifact Lineage & Storage**
  - `T1-F06-01`: Verify spatial AnnData persisted under `adata://<study>/spatial_domains/vN` with SHA-256 hash.
  - `T1-F06-02`: Verify spatial tables persisted under `table://<study>/spatial_deg/vN`.
  - `T1-F06-03`: Verify `LineageGraph` connects spatial artifacts to parent normalized/integrated AnnData.
  - `T1-F06-04`: Verify retrieval of spatial coordinates and domain annotations from persisted `.h5ad` payload.
  - `T1-F06-05`: Verify `ArtifactMetadata.type == ArtifactType.SPATIAL_DATA`.

### 3.2 External Agent Adapter Plane & Guardrails (F07 – F12)
- **F07: BaseAgentAdapter Protocol**
  - `T1-F07-01`: Verify `BaseAgentAdapter` subclass instantiation and registration in `CapabilityRegistry`.
  - `T1-F07-02`: Verify contract parameter passing and execution within `execute(contract, registry)`.
  - `T1-F07-03`: Verify adapter returns valid `TaskResult` with `TaskStatus.SUCCESS`.
  - `T1-F07-04`: Verify adapter records `executed_operations` accurately.
  - `T1-F07-05`: Verify adapter execution time and resource metadata tracking.
- **F08: SpaCellAgentAdapter**
  - `T1-F08-01`: Verify `SpaCellAgentAdapter` identifies microenvironment niches from spatial AnnData.
  - `T1-F08-02`: Verify `SpaCellAgentAdapter` creates output spatial niche table.
  - `T1-F08-03`: Verify execution with `allowed_operations=['compute_spatial_neighbors', 'identify_domains']`.
  - `T1-F08-04`: Verify adapter output preserves original cell barcodes and metadata keys.
  - `T1-F08-05`: Verify adapter links generated table to source AnnData in `ArtifactRegistry`.
- **F09: ChatCellAdapter**
  - `T1-F09-01`: Verify `ChatCellAdapter` answers cellular state transition queries for DAM microglia.
  - `T1-F09-02`: Verify `ChatCellAdapter` generates state prediction probabilities between Homeostatic and DAM states.
  - `T1-F09-03`: Verify dialogue summary output format as structured `Report` or `Table` artifact.
  - `T1-F09-04`: Verify execution under contract with `allowed_operations=['query_cell_dialogue', 'predict_state_transition']`.
  - `T1-F09-05`: Verify determinism of dialogue predictions under fixed temperature / mock backend.
- **F10: GeneAgentAdapter**
  - `T1-F10-01`: Verify `GeneAgentAdapter` performs pathway reasoning on DEG gene lists (*Trem2*, *Apoe*, *Tyrobp*).
  - `T1-F10-02`: Verify `GeneAgentAdapter` identifies functional lipid metabolism and phagocytosis annotations.
  - `T1-F10-03`: Verify generation of structured gene reasoning summary table.
  - `T1-F10-04`: Verify contract compliance with `allowed_operations=['query_gene_ontology', 'map_reactome_pathways']`.
  - `T1-F10-05`: Verify ortholog resolution from mouse to human gene identifiers.
- **F11: SideEffectValidator Guardrail**
  - `T1-F11-01`: Verify `SideEffectValidator` passes compliant tasks when `executed_operations` match `allowed_operations`.
  - `T1-F11-02`: Verify `SideEffectValidator` passes tasks when output cell count equals input cell count.
  - `T1-F11-03`: Verify `SideEffectValidator` passes tasks when cluster column `obs['leiden']` is unaltered.
  - `T1-F11-04`: Verify `SideEffectValidator` validates input artifact SHA-256 integrity post-execution.
  - `T1-F11-05`: Verify `SideEffectValidator` returns clean `(True, "")` validation status on compliant execution.
- **F12: Contract-Bound Artifact Wrapping**
  - `T1-F12-01`: Verify adapter payload is wrapped with `ArtifactMetadata` including creator task ID.
  - `T1-F12-02`: Verify content hash is computed over adapter output payload.
  - `T1-F12-03`: Verify storage backend persists adapter output to designated disk path.
  - `T1-F12-04`: Verify registration of wrapped artifact in `ArtifactRegistry`.
  - `T1-F12-05`: Verify lineage edge creation from input artifacts to adapter wrapped artifacts.

### 3.3 Multi-Source Knowledge Engine & Dual Modes (F13 – F17)
- **F13: Literature Retriever (PubMed/bioRxiv)**
  - `T1-F13-01`: Verify PubMed query formulation for Alzheimer's Disease and microglia markers (*Trem2*, *Apoe*).
  - `T1-F13-02`: Verify extraction of PMID, DOI, title, year, journal, and abstract text.
  - `T1-F13-03`: Verify retrieval of bioRxiv preprint metadata for emerging targets.
  - `T1-F13-04`: Verify deterministic cached retriever returns expected PMIDs in offline test environment.
  - `T1-F13-05`: Verify co-occurrence score computation between target genes and disease terms.
- **F14: Biological DB Retriever (GO/Reactome/NCBI)**
  - `T1-F14-01`: Verify Gene Ontology (GO) Biological Process hypergeometric enrichment test (ORA).
  - `T1-F14-02`: Verify Reactome pathway enrichment on DAM DEG list identifies "Microglial Pathogen Phagocytosis".
  - `T1-F14-03`: Verify NCBI Gene ortholog resolution between human *TREM2* (GeneID: 54209) and mouse *Trem2* (GeneID: 83433).
  - `T1-F14-04`: Verify ORA multiple testing correction via Benjamini-Hochberg ($q < 0.05$).
  - `T1-F14-05`: Verify retrieval of canonical marker sets from reference biological database.
- **F15: Discovery Mode Workflow**
  - `T1-F15-01`: Verify execution of `KnowledgeEngine.execute_discovery()` without prior hypothesis gene weighting.
  - `T1-F15-02`: Verify HVG selection and Leiden clustering proceed purely from data variance.
  - `T1-F15-03`: Verify knowledge retrieval executes post-DEG to contextualize empirical clusters.
  - `T1-F15-04`: Verify generated claims contain `epistemic_tag="unbiased_discovery"`.
  - `T1-F15-05`: Verify report summary omits confirmatory hypothesis testing callouts in discovery mode.
- **F16: Prior-Guided Mode & Epistemic Tagging**
  - `T1-F16-01`: Verify execution of `KnowledgeEngine.execute_prior_guided()` with hypothesis "DAM TREM2-APOE axis".
  - `T1-F16-02`: Verify mandatory injection of `epistemic_tag="prior-guided"` on all hypothesis-derived `ClaimNode`s.
  - `T1-F16-03`: Verify `EvidenceNode` instances record `biological_context={"mode": "prior_guided"}`.
  - `T1-F16-04`: Verify separate reporting of confirmatory vs exploratory $p$-values.
  - `T1-F16-05`: Verify Markdown report generator renders explicit warning callout: `> [!WARNING] Prior-Guided Hypothesis Testing`.
- **F17: Knowledge Evidence Extraction**
  - `T1-F17-01`: Verify generation of `EvidenceNode` with `EvidenceType.LITERATURE_SUPPORT`.
  - `T1-F17-02`: Verify generation of `EvidenceNode` with `EvidenceType.PATHWAY_ENRICHMENT`.
  - `T1-F17-03`: Verify linkage of `EvidenceNode` to source table artifacts (`table://<study>/literature_evidence/v1`).
  - `T1-F17-04`: Verify quantitative metric recording (`pmid_count`, `pathway_fdr_q_value`).
  - `T1-F17-05`: Verify insertion of knowledge evidence nodes into `EvidenceGraph`.

### 3.4 In Silico Perturbation Simulation Plane (F18 – F21)
- **F18: Genetic Perturbation Simulation (CRISPR KO/OE)**
  - `T1-F18-01`: Verify simulated CRISPR knockout of *Trem2* sets target gene expression $\approx 0$ ($\delta = 0.95$).
  - `T1-F18-02`: Verify linear GRN propagation $\mathbf{\Delta x} = (\mathbf{I} - \alpha \mathbf{A})^{-1}\mathbf{v}$ attenuates downstream DAM signature (*Apoe*, *Clec7a*, *Tyrobp*).
  - `T1-F18-03`: Verify simulated overexpression of *Trem2* increases target expression and boosts DAM signature.
  - `T1-F18-04`: Verify state reversion metric calculation: $\ge 45\%$ shift towards homeostatic baseline.
  - `T1-F18-05`: Verify output artifact `adata://<study>/perturbation_ko_trem2/v1` contains non-negative expression values.
- **F19: Drug Response & Counterfactual Simulation**
  - `T1-F19-01`: Verify CMAP drug signature loading and cosine discordance calculation $-\cos(\mathbf{s}_{disease}, \mathbf{s}_{drug})$.
  - `T1-F19-02`: Verify identification of therapeutic candidate compound with positive reversal score ($> +0.50$).
  - `T1-F19-03`: Verify counterfactual state transition modeling predicts reduction in DAM cell frequency.
  - `T1-F19-04`: Verify drug exacerbation detection for negative reversal scores ($< -0.50$).
  - `T1-F19-05`: Verify output artifact `table://<study>/drug_reversal_scores/v1`.
- **F20: Perturbation Evidence Calibration**
  - `T1-F20-01`: Verify creation of `EvidenceNode` with `EvidenceType.PERTURBATION`.
  - `T1-F20-02`: Verify `ConfidenceCalculator` updates `Causal Confidence` ($C_{causal} \in [0.20, 0.50]$ capped at $0.50$).
  - `T1-F20-03`: Verify causal status label is set to `'in_silico_perturbed'`.
  - `T1-F20-04`: Verify composite confidence formula applies 30% causal weighting when perturbation evidence is present.
  - `T1-F20-05`: Verify `LanguageEnforcer` permits prospective Level 4 hypothesis claims for simulated interventions.
- **F21: Perturbation Artifact Lineage**
  - `T1-F21-01`: Verify simulation AnnData persisted under `adata://<study>/perturbation/vN` with SHA-256 checksum.
  - `T1-F21-02`: Verify `LineageGraph` connects perturbation artifact to source microglia subset AnnData.
  - `T1-F21-03`: Verify lineage metadata records target gene, perturbation method, and network attenuation parameter.
  - `T1-F21-04`: Verify immutability of parent AnnData post-simulation.
  - `T1-F21-05`: Verify branch versioning for multiple perturbation experiments (`v6_ko_trem2`, `v6_ko_apoe`).

### 3.5 Orchestrator Integration & Multimodal Synthesis (F22 – F24)
- **F22: Dynamic DAG Planning & Routing**
  - `T1-F22-01`: Verify `ComputationalDAGPlanner` inserts spatial tasks when `manifest.data.has_spatial_coordinates=True`.
  - `T1-F22-02`: Verify `ComputationalDAGPlanner` inserts knowledge and perturbation tasks when requested in manifest.
  - `T1-F22-03`: Verify `CapabilityRouter` resolves methods for spatial, adapter, knowledge, and perturbation capabilities.
  - `T1-F22-04`: Verify topological sorting of planned DAG executes preprocessing prior to simulation.
  - `T1-F22-05`: Verify task contract input/output chaining across all planned tasks.
- **F23: Multimodal Claim Synthesis (5-Pillar Calibration)**
  - `T1-F23-01`: Verify synthesis of Claim C101 (DAM State Transition) backed by pseudobulk DEG + Trajectory.
  - `T1-F23-02`: Verify synthesis of Claim C102 (Plaque Niche Localization) backed by Spatial Moran's I ($++++$).
  - `T1-F23-03`: Verify synthesis of Claim C103 (Knowledge Convergence) backed by Literature + Reactome ($++++$).
  - `T1-F23-04`: Verify synthesis of Claim C104 (In Silico KO Reversal) backed by Perturbation evidence ($+$).
  - `T1-F23-05`: Verify 4-tier language assignment: Level 2 for DEG, Level 3 for Spatial+DEG, Level 4 for In Silico KO.
- **F24: Multimodal Provenance Reporting**
  - `T1-F24-01`: Verify `ScientificReportGenerator` renders complete Markdown report structure.
  - `T1-F24-02`: Verify sentence-to-artifact clickable provenance card links.
  - `T1-F24-03`: Verify Mermaid diagram rendering for Artifact Lineage DAG.
  - `T1-F24-04`: Verify Mermaid diagram rendering for Evidence-Claim DAG.
  - `T1-F24-05`: Verify scientific auditor sign-off summary table included in final report.

---

## 4. Tier 2: Boundary Value Analysis & Adversarial Corner Cases

Tier 2 contains **120 test cases** ($\ge 5$ per feature) evaluating extreme boundary values, fault tolerance, data corruption, unauthorized operations, and adversarial attacks.

### 4.1 Boundary & Corner Cases Specification Matrix

| Feature ID | Test ID | Adversarial / Boundary Condition | Expected System Behavior & Invariant Assertion |
|---|---|---|---|
| **F01** | `T2-F01-01` | Spatial AnnData missing `.obsm['spatial']` slot entirely | `SpatialCapability` raises `ValueError("Spatial coordinates missing in .obsm['spatial']")`; task fails with `CODE_ERROR`. |
| **F01** | `T2-F01-02` | Spatial coordinates contain `NaN` or `Inf` values | `ComputationalValidator` check `embedding_spatial_finite` fails with `ValidationSeverity.ERROR`; task halted. |
| **F01** | `T2-F01-03` | All spatial coordinates identical (zero spatial variance) | `SpatialNeighborhoodGraph` detects singular geometry, raises `ValueError("Degenerate spatial coordinates")`. |
| **F01** | `T2-F01-04` | Single spot dataset ($N = 1$) provided to spatial neighbor builder | Gracefully handled, returns empty connectivity matrix without crashing or dividing by zero. |
| **F01** | `T2-F01-05` | Mismatch between `n_obs` in expression $X$ and rows in `.obsm['spatial']` | `ComputationalValidator` check `matrix_dimension_match` fails with `ValidationSeverity.ERROR`. |
| **F02** | `T2-F02-01` | Spatial smoothing parameter $\lambda = 0.0$ (no spatial regularization) | Evaluates standard PCA clustering; outputs valid domain labels without NaN or error. |
| **F02** | `T2-F02-02` | Spatial smoothing parameter $\lambda = 1.0$ (pure spatial neighborhood) | Evaluates pure topological smoothing; ensures matrix values remain bounded. |
| **F02** | `T2-F02-03` | Disconnected spatial graph with multiple isolated singleton spots | Spatial clustering handles disconnected components without infinite loops or unhandled graph exceptions. |
| **F02** | `T2-F02-04` | Number of requested spatial domains $K \ge N_{cells}$ | Capability raises `ValueError("Cluster count cannot exceed sample count")` with clean error reporting. |
| **F02** | `T2-F02-05` | Clustering yields single monolithic domain (silhouette $< 0.0$) | `StatisticalValidator` check `clustering_silhouette` raises `ValidationSeverity.WARNING`. |
| **F03** | `T2-F03-01` | Constant gene expression across all spots (zero expression variance) | Moran's $I$ denominator $= 0$; capability returns $I = 0.0, z = 0.0, p = 1.0$ without `ZeroDivisionError`. |
| **F03** | `T2-F03-02` | Extreme spatial checkerboard pattern (maximal negative autocorrelation) | Moran's $I$ yields $I \approx -1.0$, Geary's $C > 1.5$; FDR correctly assigned. |
| **F03** | `T2-F03-03` | Spatial weight matrix $W$ with all zero weights ($S_0 = 0$) | Autocorrelation calculator catches $S_0 == 0$, returns neutral score with warning log. |
| **F03** | `T2-F03-04` | All $G$ genes tested return $p = 1.0$ (no spatial signal) | Benjamini-Hochberg yields $q = 1.0$ for all genes; SVG list is cleanly empty (0 items). |
| **F03** | `T2-F03-05` | Floating-point underflow in analytical $p$-value calculation ($z > 40$) | Two-tailed normal CDF uses logarithmic precision; returns $p = 0.0$ or minimum float without overflow exception. |
| **F04** | `T2-F04-01` | Plaque distances contain negative numbers ($\exists d_i < 0$) | Spatial DEG pre-check flags invalid physical distances, raises `ValueError`. |
| **F04** | `T2-F04-02` | Plaque distance vector length $\neq N_{cells}$ | Dimension validator fails task with `ValidationSeverity.ERROR`. |
| **F04** | `T2-F04-03` | Single donor replicate provided for spatial DEG | `StatisticalValidator` raises pseudoreplication warning; flags spatial DEG as exploratory only. |
| **F04** | `T2-F04-04` | DEG table generated with missing FDR column (`fdr_q_val` omitted) | `StatisticalValidator` check `multiple_testing_correction` fails with `ValidationSeverity.ERROR`. |
| **F04** | `T2-F04-05` | Zero genes pass FDR threshold ($q < 0.05$) | Capability returns empty table artifact with valid headers; downstream pipeline handles empty SVG set. |
| **F05** | `T2-F05-01` | Curated ligand-receptor database contains zero overlapping genes | CCI capability returns empty interaction table with informative warning; no crash. |
| **F05** | `T2-F05-02` | Target cell type has 0 cells in tissue slice | CCI calculator skips empty cell type pair gracefully, avoids zero division in $W_{spatial}$. |
| **F05** | `T2-F05-03` | Permutation count $B = 0$ specified in parameters | Parameter validator rejects $B < 10$, raises `ValueError("n_permutations must be >= 10")`. |
| **F05** | `T2-F05-04` | Extremely high expression outlier in single spot ($10^6 \times$ normal) | Expression normalization scales outlier; permutation test remains robust. |
| **F05** | `T2-F05-05` | Spatial contact matrix $W_{spatial}$ contains negative values | Capability validator checks non-negativity of weights, rejects corrupted proximity matrix. |
| **F06** | `T2-F06-01` | Attempt to overwrite existing spatial artifact `adata://AD/spatial/v1` | `ArtifactStorageBackend.save()` raises `ArtifactAlreadyExistsError`; enforces Invariant 2. |
| **F06** | `T2-F06-02` | Corrupted `.h5ad` file on disk with invalid HDF5 magic header | `PayloadSerializer.deserialize()` catches `IOError`, raises `CorruptedArtifactError`. |
| **F06** | `T2-F06-03` | Lineage graph contains cyclic reference (`v1 -> v2 -> v1`) | `LineageGraph.add_edge()` detects cycle, raises `DAGCycleError`. |
| **F06** | `T2-F06-04` | Storage payload SHA-256 does not match metadata registry hash | `ForensicValidator` flags data tampering, halts execution with `STOP_RULE`. |
| **F06** | `T2-F06-05` | Parent URI referenced in metadata does not exist in registry | `LineageGraph` check flags missing root, raises `OrphanArtifactError`. |
| **F07** | `T2-F07-01` | Third-party adapter raises unhandled Python exception during `execute()` | `CapabilityRegistry` catches exception, sets `TaskResult.status = TaskStatus.EXECUTION_FAILURE` with stack trace. |
| **F07** | `T2-F07-02` | Adapter returns `None` instead of `TaskResult` | `CapabilityRegistry` validates return type, flags `INVALID_ADAPTER_OUTPUT`. |
| **F07** | `T2-F07-03` | Adapter times out after exceeding contract `max_execution_seconds` | Execution harness terminates worker, sets `TaskStatus.TIMEOUT`. |
| **F07** | `T2-F07-04` | Adapter requests unsupported artifact modality (`IMAGE_TIF` in scRNA) | Pre-execution contract validation rejects incompatible modality. |
| **F07** | `T2-F07-05` | Adapter attempts to execute shell commands or write to unauthorized paths | Sandboxed execution / filesystem guardrail blocks access; flags security violation. |
| **F08** | `T2-F08-01` | `SpaCellAgentAdapter` stealthily drops 50 outlier cells without contract | `SideEffectValidator` detects $n_{obs}$ changed ($1200 \rightarrow 1150$); sets `TaskStatus.POLICY_VIOLATION`. |
| **F08** | `T2-F08-02` | `SpaCellAgentAdapter` attempts unauthorized global reclustering | `SideEffectValidator` detects mutated `obs['leiden']`; sets `TaskStatus.POLICY_VIOLATION`. |
| **F08** | `T2-F08-03` | `SpaCellAgentAdapter` executes operation not in `allowed_operations` | `SideEffectValidator` detects unwhitelisted op; sets `TaskStatus.POLICY_VIOLATION`. |
| **F08** | `T2-F08-04` | `SpaCellAgentAdapter` receives empty spatial AnnData ($0$ cells) | Pre-execution check fails with `ComputationalValidator.matrix_non_empty`. |
| **F08** | `T2-F08-05` | `SpaCellAgentAdapter` attempts in-place mutation of parent artifact | SHA-256 hash mismatch of input artifact detected; sets `POLICY_VIOLATION`. |
| **F09** | `T2-F09-01` | `ChatCellAdapter` queried for non-existent cell type (`"AlienNeuron"`) | Adapter returns clean fallback message without hallucinating synthetic gene profiles. |
| **F09** | `T2-F09-02` | `ChatCellAdapter` attempts to mutate raw gene expression values | `SideEffectValidator` detects raw expression tampering; sets `POLICY_VIOLATION`. |
| **F09** | `T2-F09-03` | `ChatCellAdapter` prompt injection attack in user query string | Sanitizer strips prompt override instructions; adapter operates within biological bounds. |
| **F09** | `T2-F09-04` | `ChatCellAdapter` produces non-deterministic output across identical seeds | Validation test asserts identical token/probability outputs when `seed=42`. |
| **F09** | `T2-F09-05` | `ChatCellAdapter` returns malformed JSON prediction payload | JSON schema validator rejects payload, marks task as `EXECUTION_FAILURE`. |
| **F10** | `T2-F10-01` | `GeneAgentAdapter` receives gene list with invalid symbols (`"XYZ123"`, `""`) | Adapter filters unmapped symbols, reports unmapped rate in check metrics. |
| **F10** | `T2-F10-02` | `GeneAgentAdapter` attempts to recalculate DEG $p$-values | `SideEffectValidator` detects `recalculate_deg` in forbidden ops; sets `POLICY_VIOLATION`. |
| **F10** | `T2-F10-03` | `GeneAgentAdapter` receives empty gene list (0 genes) | Adapter returns empty pathway table with warning; avoids crash. |
| **F10** | `T2-F10-04` | `GeneAgentAdapter` cross-species lookup fails (no human ortholog) | Adapter flags gene as species-specific, preserves native mouse symbol. |
| **F10** | `T2-F10-05` | `GeneAgentAdapter` produces circular pathway graph | Knowledge DAG builder rejects cyclic pathway dependencies. |
| **F11** | `T2-F11-01` | Third-party agent attempts subtle single-cell barcode deletion | `SideEffectValidator` detects $N_{out} = N_{in} - 1$; triggers `POLICY_VIOLATION`. |
| **F11** | `T2-F11-02` | Third-party agent swaps cluster labels of two cells | `SideEffectValidator` compares cluster array element-wise; triggers `POLICY_VIOLATION`. |
| **F11** | `T2-F11-03` | Third-party agent injects `recluster` into `executed_operations` | `SideEffectValidator` intercepts denylisted operation; triggers `POLICY_VIOLATION`. |
| **F11** | `T2-F11-04` | Third-party agent modifies `.var` index names | `SideEffectValidator` detects feature index alteration; triggers `POLICY_VIOLATION`. |
| **F11** | `T2-F11-05` | Third-party agent mutates AnnData `.uns` dictionary metadata | `SideEffectValidator` checks metadata contract preservation; triggers `POLICY_VIOLATION`. |
| **F12** | `T2-F12-01` | Adapter output payload is non-serializable object (e.g. open socket/file handle)| `PayloadSerializer` raises `TypeError`; artifact creation safely aborted. |
| **F12** | `T2-F12-02` | Target disk partition is full / read-only | `ArtifactStorageBackend` catches `IOError`, returns clean `STORAGE_FAILURE` result. |
| **F12** | `T2-F12-03` | Duplicate registration of identical artifact URI | Storage backend raises `ArtifactAlreadyExistsError`. |
| **F12** | `T2-F12-04` | Empty payload provided to artifact wrapper | Wrapper rejects empty payload, raises `ValueError("Payload cannot be empty")`. |
| **F12** | `T2-F12-05` | Metadata creation missing required creator task ID | Pydantic schema validation raises `ValidationError`. |
| **F13** | `T2-F13-01` | PubMed API returns HTTP 500 / network timeout | `LiteratureRetriever` falls back to deterministic local knowledge cache without crashing. |
| **F13** | `T2-F13-02` | Query string contains SQL/Lucene injection special characters (`AND OR NOT ""`) | Query builder sanitizes search terms safely. |
| **F13** | `T2-F13-03` | Literature search returns 0 matching publications | Retriever returns empty evidence list; downstream synthesizer marks literature support as neutral. |
| **F13** | `T2-F13-04` | Duplicate PMIDs returned across multiple query terms | Retriever deduplicates citations by unique PMID/DOI. |
| **F13** | `T2-F13-05` | Publication abstract contains corrupted non-UTF-8 bytes | Text parser decodes with replacement characters; avoids decoding crash. |
| **F14** | `T2-F14-01` | GO enrichment query provided with single gene | Hypergeometric test executes with $k=1, n=1$; returns valid $p$-value without crash. |
| **F14** | `T2-F14-02` | DEG gene set larger than total genome universe ($n > M$) | Validator detects invalid universe parameters, raises `ValueError`. |
| **F14** | `T2-F14-03` | Reactome pathway database file missing or corrupted | Retriever raises informative configuration error; initiates fallback cache. |
| **F14** | `T2-F14-04` | All genes in DEG set unmapped to any known pathway | Hypergeometric test returns empty enrichment table; downstream pipeline continues safely. |
| **F14** | `T2-F14-05` | Division by zero in hypergeometric test when universe $M = 0$ | Pre-check asserts $M > 0$, prevents math error. |
| **F15** | `T2-F15-01` | In Discovery Mode, user supplies biased gene weighting list | `DiscoveryModePolicy` strips biased gene weights, logs policy enforcement warning. |
| **F15** | `T2-F15-02` | Discovery Mode attempts to force clustering resolution based on prior literature| `DiscoveryModePolicy` rejects literature-directed clustering; relies on silhouette optimization. |
| **F15** | `T2-F15-03` | Discovery Mode claim improperly tagged as `prior-guided` | Epistemic auditor detects tag contradiction; resets tag to `unbiased_discovery`. |
| **F15** | `T2-F15-04` | Discovery Mode DEG list contains 0 significant genes | Downstream claim engine notes absence of differential signal; avoids false claims. |
| **F15** | `T2-F15-05` | Discovery Mode pipeline attempts to skip independent auditing | Orchestrator enforces Invariant 4; mandatory auditing cannot be bypassed. |
| **F16** | `T2-F16-01` | Prior-guided analysis executed but `epistemic_tag` omitted from claim | `EpistemicAuditor` halts execution, raises `EpistemicTagMissingError`. |
| **F16** | `T2-F16-02` | User hypothesis references non-existent gene in dataset | Knowledge engine flags missing target gene; reports hypothesis as untestable on dataset. |
| **F16** | `T2-F16-03` | Prior-guided claim uses exploratory $p$-value as confirmatory | Statistical auditor rejects confirmatory claim; forces exploratory labeling. |
| **F16** | `T2-F16-04` | Report generator omits prior-guided warning callout in prior-guided study | Report validation test fails if warning badge is missing from document text. |
| **F16** | `T2-F16-05` | Contradictory evidence found against user prior hypothesis | Evidence graph links contradictory evidence node with penalty; reduces overall confidence. |
| **F17** | `T2-F17-01` | Knowledge evidence node created with negative confidence score | Pydantic schema validation rejects `score < 0.0` or `score > 1.0`. |
| **F17** | `T2-F17-02` | Evidence node missing source artifact URI | Schema validation rejects missing `source_artifact_uris`. |
| **F17** | `T2-F17-03` | Evidence graph cycle introduced via knowledge node | DAG validator detects cycle, raises `DAGCycleError`. |
| **F17** | `T2-F17-04` | Evidence strength set to `INSUFFICIENT` | Confidence calculator assigns 0.0 weight to insufficient evidence node. |
| **F17** | `T2-F17-05` | Duplicate evidence nodes with identical ID inserted into graph | Graph rejects duplicate node ID, raises `DuplicateNodeError`. |
| **F18** | `T2-F18-01` | Genetic KO simulated for target gene not in dataset `var_names` | Capability raises `KeyError(f"Target gene '{gene}' not found in expression matrix")`. |
| **F18** | `T2-F18-02` | GRN adjacency matrix spectral radius $\rho(\mathbf{A}) \ge 1.0$ (divergent series) | Capability normalizes matrix eigenvalues to guarantee $(\mathbf{I} - \alpha \mathbf{A})^{-1}$ convergence. |
| **F18** | `T2-F18-03` | Perturbation attenuation parameter $\alpha < 0.0$ or $\alpha \ge 1.0$ | Parameter validator raises `ValueError("network_attenuation must be in [0.0, 0.99]")`. |
| **F18** | `T2-F18-04` | Simulated expression yields negative counts ($X_{perturbed} < 0$) | Capability applies $\max(0, \mathbf{X} + \mathbf{\Delta X})$ non-negativity constraint. |
| **F18** | `T2-F18-05` | Perturbation simulation yields `NaN` or `Inf` in output matrix | `ComputationalValidator` catches NaNs; task fails with `ValidationSeverity.ERROR`. |
| **F19** | `T2-F19-01` | Drug signature vector contains zero variance (all values $= 0$) | Cosine discordance returns $0.0$; capability logs neutral drug response. |
| **F19** | `T2-F19-02` | Drug signature gene set has 0 overlap with disease DEG set | Capability raises `ValueError("Zero overlapping genes between drug and disease signatures")`. |
| **F19** | `T2-F19-03` | Drug response simulation yields negative cell state transition probabilities | Softmax / normalization enforces $\sum p_k = 1.0$ and $p_k \ge 0.0$. |
| **F19** | `T2-F19-04` | Compound signature contains `NaN` values | Pre-check detects NaNs, raises `ValueError`. |
| **F19** | `T2-F19-05` | Reversal score exactly $+1.0$ or $-1.0$ (collinear signatures) | Handled smoothly without precision overflow; valid cosine bounds $[-1.0, 1.0]$ preserved. |
| **F20** | `T2-F20-01` | Perturbation evidence node attempts to set `causal_status = "experimental_perturbed"` | Auditor checks provenance: if wet-lab assay missing, forces status to `'in_silico_perturbed'`. |
| **F20** | `T2-F20-02` | In silico perturbation causal confidence score computed $> 0.50$ | `ConfidenceCalculator` strictly caps in silico causal confidence at $0.50$. |
| **F20** | `T2-F20-03` | Observational claim upgraded to Level 1 / Level 2 with causal verbs post-simulation| `LanguageEnforcer` catches causal verb (e.g. *"Apoe drives DAM"*); throws `EPISTEMIC VIOLATION`. |
| **F20** | `T2-F20-04` | Contradicting perturbation evidence (KO exacerbates disease state) | Confidence calculator applies contradiction penalty, lowering overall claim confidence. |
| **F20** | `T2-F20-05` | Perturbation simulation stability score $< 0.50$ under replicate bootstrap | Auditor flags unstable simulation; downgrades evidence strength to `WEAK`. |
| **F21** | `T2-F21-01` | Attempt to overwrite existing simulation artifact `adata://AD/perturbation/v1` | `ArtifactStorageBackend` raises `ArtifactAlreadyExistsError`. |
| **F21** | `T2-F21-02` | Perturbation artifact lineage missing parent microglia subset link | Lineage integrity auditor flags missing parent edge. |
| **F21** | `T2-F21-03` | Corrupted simulation AnnData payload on disk | Deserializer catches `IOError`, halts loading. |
| **F21** | `T2-F21-04` | Multiple KO branches created with identical branch name | Registry enforces unique branch naming (`v1_ko_trem2`, `v1_ko_apoe`). |
| **F21** | `T2-F21-05` | Simulation metadata missing random seed | Schema validation enforces `random_seed` recording for reproducibility. |
| **F22** | `T2-F22-01` | `ComputationalDAGPlanner` receives cyclic task dependency specification | Planner detects cycle, raises `DAGCycleError`. |
| **F22** | `T2-F22-02` | `CapabilityRouter` receives unregistered capability request (`"magic_cure_v1"`) | Router raises `KeyError("Unregistered capability")` with available method suggestions. |
| **F22** | `T2-F22-03` | Task dependency references non-existent upstream task ID | DAG validator catches broken dependency link before execution begins. |
| **F22** | `T2-F22-04` | Orchestrator receives `STOP_RULE` from statistical auditor | Orchestrator halts execution of downstream dependent tasks immediately. |
| **F22** | `T2-F22-05` | Maximum retry count exceeded on failed task | Orchestrator transitions study status to `FAILED_AUDIT`, generates partial report. |
| **F23** | `T2-F23-01` | Level 1 observational claim contains forbidden verb *"causes"* | `LanguageEnforcer` throws `EPISTEMIC VIOLATION: Statement contains strong causal verb`. |
| **F23** | `T2-F23-02` | Level 2 statistical claim contains forbidden verb *"proves"* | `LanguageEnforcer` throws `EPISTEMIC VIOLATION: Statement contains uncalibrated proof claim`. |
| **F23** | `T2-F23-03` | Claim confidence computed as negative ($C_{overall} < 0.0$) due to penalties | `ConfidenceCalculator` clamps overall confidence to minimum $0.0$. |
| **F23** | `T2-F23-04` | Claim linked to 0 supporting evidence nodes | `ClaimEngine` rejects ungrounded claim, raises `UngroundedClaimError`. |
| **F23** | `T2-F23-05` | Claim confidence score exceeds $1.0$ due to weight accumulation | `ConfidenceCalculator` clamps overall confidence to maximum $1.0$. |
| **F24** | `T2-F24-01` | Report generator references deleted or missing artifact URI | Report generator catches `KeyError`, renders missing artifact placeholder with warning. |
| **F24** | `T2-F24-02` | Empty Evidence DAG passed to report generator | Generator renders valid report structure with "No Claims Synthesized" section. |
| **F24** | `T2-F24-03` | Special characters in study title break Markdown rendering (`<script>`, `#`, `|`) | Report generator escapes HTML/Markdown special characters. |
| **F24** | `T2-F24-04` | Mermaid diagram syntax error generated by malformed node IDs | Generator sanitizes node IDs (replacing spaces/slashes with underscores). |
| **F24** | `T2-F24-05` | Sentence provenance tracker encounters unmapped sentence | Tracker logs unmapped statement, assigns generic study header provenance. |

---

## 5. Tier 3: Cross-Feature Pairwise Combinatorial Interactions

Tier 3 contains **8 deep cross-plane interaction tests** verifying contract guardrails, data flow consistency, and epistemic synchronization across subsystems.

```
┌────────────────────────────────────────────────────────────────────────────────────────┐
│                        TIER 3 PAIRWISE COMBINATORIAL INTERACTIONS                      │
├────────────────────────────────────────────────────────────────────────────────────────┤
│ T3-PAIR-01: Spatial Analytics (F01-F06) × SpaCellAgentAdapter (F08) × Guardrails (F11) │
│ T3-PAIR-02: Spatial DEG (F04) × Knowledge Engine Literature/BioDB (F13, F14)          │
│ T3-PAIR-03: Spatial CCI (F05) × In Silico Genetic Perturbation Simulation (F18)        │
│ T3-PAIR-04: ChatCellAdapter (F09) × In Silico Compound Response Simulation (F19)       │
│ T3-PAIR-05: GeneAgentAdapter (F10) × In Silico CRISPR KO Target Selection (F18)        │
│ T3-PAIR-06: Prior-Guided Mode (F16) × Multimodal Claim Synthesis (F23) × Report (F24)  │
│ T3-PAIR-07: Dynamic DAG Planner (F22) × SideEffectValidator (F11) × Lineage Graph (F06)│
│ T3-PAIR-08: 5-Pillar Evidence (F23) × Contradiction Penalty (F20) × Language (F23)     │
└────────────────────────────────────────────────────────────────────────────────────────┘
```

### 5.1 Pairwise Interaction Test Cases

#### `T3-PAIR-01`: Spatial Analytics $\times$ `SpaCellAgentAdapter` $\times$ `SideEffectValidator`
- **Objective**: Verify that `SpaCellAgentAdapter` executes spatial domain clustering on `adata://.../spatial/vN` under strict contract bounds, and that rogue attempts to recluster or drop cells are blocked with `POLICY_VIOLATION`.
- **Inputs**: Spatial AnnData ($N = 600$ cells, 2D coordinates in `.obsm['spatial']`).
- **Contract**: `allowed_operations=['compute_spatial_neighbors', 'identify_domains']`, `forbidden_operations=['filter_cells', 'recluster_global']`.
- **Assertions**:
  1. Compliant execution produces `adata://.../spatial_domains/vN` and `table://.../spatial_niches/v1`.
  2. Rogue adapter attempting `filter_cells` is blocked; `TaskResult.status == TaskStatus.POLICY_VIOLATION` with `ExecutionFailureType.UNAUTHORIZED_SIDE_EFFECT`.
  3. Raw data `raw/v1` remains unaltered with identical SHA-256 hash.

#### `T3-PAIR-02`: Spatial DEG $\times$ Knowledge Engine (`LiteratureRetriever` + `BiologicalDBRetriever`)
- **Objective**: Verify that Spatially Variable Genes (SVGs) detected via Moran's $I$ are automatically piped into GO/Reactome ORA and PubMed literature retrieval.
- **Inputs**: Spatial DEG table containing top SVGs (*Apoe*, *Trem2*, *Cst7*, *Tyrobp*, *C3*).
- **Assertions**:
  1. `BiologicalDBRetriever` yields significant GO term "microglial cell activation" ($q < 0.01$).
  2. `LiteratureRetriever` returns PMIDs matching spatial amyloid plaque colocalization.
  3. Structured `EvidenceNode` instances for both `SPATIAL_LOCALIZATION` and `PATHWAY_ENRICHMENT` are inserted into `EvidenceGraph`.

#### `T3-PAIR-03`: Spatial CCI $\times$ In Silico Genetic Perturbation Simulation
- **Objective**: Verify that key ligand-receptor interactions identified by spatial cell-cell communication (*Apoe-Trem2*) are perturbed in silico, evaluating disruption of spatial signaling.
- **Inputs**: Spatial AnnData with active *Apoe-Trem2* interaction score $S_{spatial} = 0.82$.
- **Assertions**:
  1. In silico KO of *Trem2* reduces downstream cell-cell communication interaction score to $S_{spatial} \approx 0.05$.
  2. Simulation yields versioned artifact `adata://<study>/perturbation_ko_trem2/vN`.
  3. `EvidenceType.PERTURBATION` evidence node created and linked to spatial CCI claim.

#### `T3-PAIR-04`: `ChatCellAdapter` $\times$ In Silico Compound Response Simulation
- **Objective**: Verify that `ChatCellAdapter` can query post-perturbation counterfactual states and generate conversational summaries of drug reversal efficacy.
- **Inputs**: Simulated compound response AnnData with positive reversal score ($+0.68$).
- **Assertions**:
  1. `ChatCellAdapter` queries state transition probabilities from DAM to Homeostatic microglia.
  2. Adapter generates dialogue summary indicating high likelihood of phenotypic reversion.
  3. Output wrapped as `table://<study>/chatcell_reversal_summary/v1` with verified SHA-256.

#### `T3-PAIR-05`: `GeneAgentAdapter` $\times$ In Silico CRISPR KO Target Selection
- **Objective**: Verify that `GeneAgentAdapter` pathway reasoning on donor pseudobulk DEGs identifies *Trem2* as a key hub gene, which dynamically triggers CRISPR KO simulation.
- **Inputs**: Pseudobulk DEG table containing 50 upregulated genes in AD microglia.
- **Assertions**:
  1. `GeneAgentAdapter` highlights *Trem2* as central regulator in phagocytosis subnetwork.
  2. Orchestrator routes target *Trem2* to `GeneticPerturbationCapability`.
  3. GRN network propagation attenuates 80% of connected downstream target genes.

#### `T3-PAIR-06`: Prior-Guided Mode $\times$ Multimodal Claim Synthesis $\times$ Report Generator
- **Objective**: Verify end-to-end propagation of the `[PRIOR-GUIDED]` epistemic badge from `StudyManifest` through `KnowledgeEngine`, `ClaimEngine`, to `ScientificReportGenerator`.
- **Inputs**: `StudyManifest` with `analysis_policy.prior_guided_analysis = True` and hypothesis `"DAM TREM2-APOE axis"`.
- **Assertions**:
  1. `ClaimEngine` creates `ClaimNode` with `epistemic_tag="prior-guided"`.
  2. `ScientificReportGenerator` includes warning callout banner in rendered Markdown.
  3. Level 3 Supported Interpretation claim explicitly references prior hypothesis origin.

#### `T3-PAIR-07`: Dynamic DAG Planner $\times$ `SideEffectValidator` $\times$ Lineage Graph
- **Objective**: Verify that if an intermediate capability in a multi-plane DAG triggers a contract violation, the orchestrator halts downstream dependent tasks while preserving complete upstream lineage.
- **Inputs**: Planned 10-task DAG where Task 5 (`SpaCellAgentAdapter`) is injected with a rogue cell filtering behavior.
- **Assertions**:
  1. Task 5 fails with `TaskStatus.POLICY_VIOLATION`.
  2. Tasks 6–10 are not executed (status `CANCELLED` / `BLOCKED`).
  3. Upstream artifacts (Tasks 1–4) remain valid, immutable, and fully queryable in `LineageGraph`.

#### `T3-PAIR-08`: 5-Pillar Multimodal Evidence $\times$ Contradiction Penalties $\times$ Language Enforcer
- **Objective**: Verify that when empirical evidence strongly supports a claim ($++++$ spatial, $++++$ DEG, $++++$ literature) but an in silico perturbation produces contradictory data, the confidence calculator applies contradiction penalties and `LanguageEnforcer` adjusts the language tier.
- **Inputs**: Supporting spatial and DEG evidence ($s=0.90$) + Contradicting perturbation evidence ($s=0.85$, polarity=`CONTRADICTS`).
- **Assertions**:
  1. Contradiction penalty $P_{contra} = 0.85 \times 0.35 \times 0.25 \approx 0.074$ is subtracted from overall confidence.
  2. Confidence reduced from $0.88 \rightarrow 0.80$.
  3. `LanguageEnforcer` formats claim at Level 3 (Supported Interpretation) noting the simulation discordance.

---

## 6. Tier 4: Real-World Application Scenario (Full E2E Scientific Study)

Tier 4 specifies a comprehensive end-to-end scientific study: **Alzheimer's Disease Microglial Disease-Associated Microglia (DAM) Activation Axis**.

### 6.1 Study Specification
- **Study ID**: `AD_microglia_spatial_perturb_study_001`
- **Species**: *Mus musculus*
- **Disease**: Alzheimer's Disease (5xFAD transgenic mouse model vs Wild-Type Control)
- **Experimental Design**: 12 biological donor mice (6 5xFAD AD, 6 WT Ctrl, balanced across 2 technical batches).
- **Modalities**: Dual-modality scRNA-seq ($N = 1,200$ cells, $G = 500$ genes) + Spatial Transcriptomics (2D coordinates in $\mu\text{m}$, plaque proximity annotations).
- **Analysis Policy**: Dual-mode evaluation (Discovery Mode baseline + Prior-Guided DAM TREM2-APOE hypothesis testing).

### 6.2 18-Step End-to-End Workflow Execution Trace

```text
Step 01: User Prompt Ingestion & Intent Parsing
         Prompt: "Investigate Alzheimer's microglial spatial niches and test the DAM Trem2-Apoe axis."
         ──► StudyManifest created (prior_guided=True, has_spatial_coordinates=True).
         
Step 02: Raw Data Ingestion & Immutability Lock
         ──► Ingests 1,200 cells with .obsm['spatial'] to adata://AD_001/raw/v1 (SHA-256 locked).

Step 03: Dataset Audit Task
         ──► DatasetAuditCapability verifies 12 donor replicates, 2 batches -> table://AD_001/audit/v1.

Step 04: Quality Control & Cell Filtering
         ──► QCCapability filters low-quality cells -> adata://AD_001/qc/v1.

Step 05: Log1p Normalization & HVG Selection
         ──► NormalizationCapability (10,000 CPM + log1p) -> adata://AD_001/norm/v2.

Step 06: Harmony Batch Integration
         ──► IntegrationCapability corrects technical batch variance -> adata://AD_001/integrated/v3.

Step 07: Leiden Clustering & Microglia Subsetting
         ──► ClusteringCapability annotates Microglia, Astrocytes, Neurons, Oligodendrocytes.
         ──► SubsetCapability extracts Microglial population -> adata://AD_001/microglia/v5.

Step 08: Spatial Domain Identification via SpaCellAgentAdapter
         ──► SpaCellAgentAdapter clusters Plaque-Adjacent vs Distal microenvironments
         ──► adata://AD_001/spatial_domains/v6 (Audited by SpatialValidator).

Step 09: Spatial Autocorrelation & Spatial DEG
         ──► SpatialDEGCapability computes Moran's I for all genes.
         ──► Apoe (I=0.52, q<1e-6), Trem2 (I=0.48, q<1e-5), Clec7a (I=0.41, q<1e-4).
         ──► table://AD_001/spatial_moran/v1.

Step 10: Spatial Cell-Cell Communication (CCI)
         ──► CellCellCommunicationCapability evaluates Apoe-Trem2 and App-Cd74 LR pairs.
         ──► table://AD_001/spatial_cci/v1.

Step 11: Donor-Level Pseudobulk DEG
         ──► DifferentialExpressionCapability aggregates counts per donor (6 AD vs 6 Ctrl).
         ──► StatisticalValidator verifies pseudobulk design (confirmatory=True).
         ──► table://AD_001/deg_pseudobulk/v1.

Step 12: Trajectory Inference & Stability Testing
         ──► TrajectoryCapability computes PAGA pseudotime from Homeostatic -> DAM microglia.
         ──► 5-fold bootstrap stability test yields score = 0.78 (passes >= 0.60 stop rule).
         ──► table://AD_001/trajectory/v1.

Step 13: Multi-Source Knowledge Retrieval (Prior-Guided Mode)
         ──► KnowledgeEngine queries PubMed (retrieves PMIDs on DAM Trem2 axis).
         ──► BiologicalDBRetriever runs Reactome ORA ("Microglial Pathogen Phagocytosis", q<1e-4).
         ──► Mandatory tag applied: [PRIOR-GUIDED HYPOTHESIS TESTING: DAM TREM2-APOE AXIS].
         ──► table://AD_001/literature/v1, table://AD_001/pathways/v1.

Step 14: GeneAgentAdapter Functional Reasoning
         ──► GeneAgentAdapter maps Trem2-Tyrobp-Apoe signaling cascade.
         ──► table://AD_001/gene_reasoning/v1.

Step 15: In Silico Genetic Perturbation Simulation (CRISPR KO of Trem2)
         ──► GeneticPerturbationCapability simulates Trem2 KO via GRN linear propagation.
         ──► Attenuates DAM state signature by 58% towards homeostatic baseline.
         ──► adata://AD_001/perturbation_ko_trem2/v7.

Step 16: In Silico Pharmacological Counterfactual Simulation & ChatCell Dialogue
         ──► CompoundPerturbationCapability simulates neuroinflammation-reversing compound.
         ──► ChatCellAdapter predicts 72% reversion probability.
         ──► table://AD_001/drug_reversal/v1, table://AD_001/chatcell_summary/v1.

Step 17: Independent Scientific Auditing (All 4 Validators)
         ──► ComputationalValidator: Zero NaNs/Infs in X and .obsm['spatial'].
         ──► StatisticalValidator: Pseudobulk confirmatory=True, FDR controlled, stability=0.78.
         ──► BiologicalValidator: Canonical microglial markers (Tmem119, Hexb, Cx3cr1) coherent.
         ──► ForensicValidator: All 14 artifact SHA-256 hashes match registry records.

Step 18: 5-Pillar Evidence DAG Assembly & Multimodal Claim Synthesis
         ──► 5 Evidence Pillars synthesized:
             - E_spatial (Moran's I = 0.52, ++++)
             - E_deg (Pseudobulk log2FC = 2.4, q < 1e-6, ++++)
             - E_traj (Stability = 0.78, +++)
             - E_lit (PubMed PMIDs >= 3, Reactome q < 1e-4, ++++)
             - E_perturb (In silico Trem2 KO reversion = 58%, +)
         ──► 4-Tier Traceable Claims generated with strict causal language enforcement.

Step 19: Provenance-Tracked Scientific Report Generation
         ──► ScientificReportGenerator outputs complete Markdown report with clickable sentence
             provenance cards, prior-guided warning callouts, and Mermaid lineage/evidence DAGs.
```

### 6.3 E2E Success & Verification Assertions
1. **Zero Policy Violations**: All 14 computational tasks complete with `TaskStatus.SUCCESS`.
2. **Four Invariants Fully Satisfied**:
   - `adata://AD_001/raw/v1` remains 100% byte-identical to original ingest.
   - All 14 artifacts have unique SHA-256 checksums and unbroken parent lineage in `LineageGraph`.
   - All claims are backed by $\ge 1$ supporting `EvidenceNode`; observational claims contain zero banned causal verbs.
   - Independent scientific auditors pass without `STOP_RULE` or unhandled `ERROR`.
3. **Multimodal Evidence Convergence**: Calculated overall confidence for the DAM activation claim $C_{overall} \ge 0.85$.
4. **Epistemic Traceability**: Clickable sentence provenance links in Markdown report navigate directly to verified artifact SHA-256 checksums.

---

## 7. Test Execution Protocol & Pass Criteria

### 7.1 Test Execution Command
The test suite is executed using the standard pytest runner in the active virtual environment:

```powershell
.venv\Scripts\pytest -v --basetemp=.pytest_temp
```

### 7.2 Pass Criteria & Quality Gates
- **Pass Rate**: **100% pass** (0 failed, 0 errors) across all test files.
- **Feature Coverage**: 100% coverage of all 24 features across Tiers 1–4.
- **Zero Regressions**: All baseline tests in `tests/test_*.py` pass synchronously.
- **Execution Time**: Complete test suite finishes within $< 60\text{ seconds}$ on standard runner hardware.
- **Audit Compliance**: All test artifacts pass independent computational, statistical, and biological validation.
