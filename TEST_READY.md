# EACBP Test Suite Readiness Report (TEST_READY.md)

**Date**: 2026-08-28  
**Author**: Test Writer 1 (E2E Testing Track)  
**System**: Evidence-aware Agentic Computational Biology Platform (EACBP) — Phases V3–V5 Scientific Workflow OS  
**Status**: **TEST INFRASTRUCTURE & ARCHITECTURE READY**

---

## 1. Executive Summary

The E2E test infrastructure and architecture for EACBP V3–V5 has been fully specified and structured according to the rigorous **4-Tier Test Methodology**:
1. **Tier 1 (Feature Coverage & Happy-Path Isolation)**: 120 test cases ($\ge 5$ per feature across all 24 features).
2. **Tier 2 (Boundary Value Analysis & Corner Cases)**: 120 test cases ($\ge 5$ per feature across all 24 features).
3. **Tier 3 (Cross-Feature Pairwise Combinatorial Interactions)**: 8 deep cross-plane interaction tests.
4. **Tier 4 (Real-World Application Workload)**: 18-step full end-to-end scientific study on the Alzheimer's Disease Microglial DAM Activation Axis.

The complete architectural specifications, invariant test oracles, and feature inventories are published in `TEST_INFRA.md` at the project root.

---

## 2. Test Suite File Structure & Target Mapping

The test suite in `tests/` is organized into modular test files mapped to each subsystem and milestone:

| Test File | Target Subsystem / Milestone | Features Covered | Test Types |
|---|---|---|---|
| `tests/test_schemas.py` | Unified Schemas (`eacbp/schemas/`) | F01, F07, F12, F16, F20 | Unit, Schema validation, Pydantic bounds |
| `tests/test_artifacts.py` | Artifact Storage & Lineage (`eacbp/artifact/`) | F06, F12, F21 | Immutability, SHA-256 hashing, DAG lineage |
| `tests/test_capabilities.py` | Compute Capabilities & Guardrails (`eacbp/capabilities/`) | F11, Base capabilities | Contract bounds, SideEffectValidator |
| `tests/test_spatial.py` | Spatial Single-Cell Analytics (`eacbp/capabilities/spatial/` - M1) | F01, F02, F03, F04, F05, F06 | Spatial coordinates, Moran's I, Spatial DEG, CCI |
| `tests/test_adapters.py` | External Agent Adapters (`eacbp/adapters/` - M2) | F07, F08, F09, F10, F11, F12 | SpaCell, ChatCell, GeneAgent, rogue blocking |
| `tests/test_knowledge.py` | Knowledge Engine & Dual Modes (`eacbp/knowledge/` - M3) | F13, F14, F15, F16, F17 | PubMed, Reactome, Discovery vs Prior-guided |
| `tests/test_perturbation.py` | In Silico Perturbation Plane (`eacbp/capabilities/perturbation/` - M4) | F18, F19, F20, F21 | CRISPR KO, drug response, counterfactuals |
| `tests/test_auditors.py` | Independent Scientific Auditor (`eacbp/auditor/`) | F01-F24 Verification | Computational (NaNs), Statistical (FDR), Biological |
| `tests/test_evidence_graph.py`| Evidence DAG & Claim Engine (`eacbp/evidence/`) | F17, F20, F23 | 5-Pillar confidence, LanguageEnforcer |
| `tests/test_orchestrator.py` | Orchestration & Dynamic DAG (`eacbp/orchestrator/` - M5) | F22, F23, F24 | Intent parsing, 3-tier router, DAG planner |
| `tests/test_end_to_end_study.py`| Full Integrated Scientific Study (M5–M6) | F01–F24 Integrated | 18-step E2E study (scRNA + Spatial + Adapters + Knowledge + Perturbation) |

---

## 3. Test Runner Commands & Environment Configuration

### 3.1 Primary Test Execution Command
```powershell
.venv\Scripts\pytest -v --basetemp=.pytest_temp
```

### 3.2 Targeted Subsystem Test Commands
- **Spatial Plane**:
  ```powershell
  .venv\Scripts\pytest -v tests/test_spatial.py --basetemp=.pytest_temp
  ```
- **Agent Adapter Plane**:
  ```powershell
  .venv\Scripts\pytest -v tests/test_adapters.py --basetemp=.pytest_temp
  ```
- **Knowledge Plane**:
  ```powershell
  .venv\Scripts\pytest -v tests/test_knowledge.py --basetemp=.pytest_temp
  ```
- **In Silico Perturbation Plane**:
  ```powershell
  .venv\Scripts\pytest -v tests/test_perturbation.py --basetemp=.pytest_temp
  ```
- **Full E2E Scientific Study**:
  ```powershell
  .venv\Scripts\pytest -v tests/test_end_to_end_study.py --basetemp=.pytest_temp
  ```

---

## 4. Acceptance Criteria Verification Matrix

| Acceptance Criteria ID | Description | Planned Test Cases | Target Threshold |
|---|---|---|---|
| **AC-1** | Spatial domain, spatial DEG, and neighborhood analysis execute and produce immutable, hashed artifacts with verified lineage. | `test_spatial.py`, `test_artifacts.py` | 100% pass, SHA-256 verified, lineage connected |
| **AC-2** | Computational and statistical auditors evaluate spatial metrics without regressions. | `test_auditors.py`, `test_spatial.py` | Zero NaNs, Moran's I analytical bounds $[-1, 1]$, FDR $< 0.05$ |
| **AC-3** | `SpaCellAgentAdapter`, `ChatCellAdapter`, and `GeneAgentAdapter` run within `TaskContract` bounds. | `test_adapters.py` | `TaskStatus.SUCCESS`, outputs wrapped in artifacts |
| **AC-4** | Any agent violating forbidden operations (e.g. unauthorized reclustering) is blocked with a `POLICY_VIOLATION` task status. | `test_capabilities.py`, `test_adapters.py` | `TaskStatus.POLICY_VIOLATION`, `UNAUTHORIZED_SIDE_EFFECT` |
| **AC-5** | Literature and Biological DB retrieval modules provide structured evidence nodes. | `test_knowledge.py`, `test_evidence_graph.py` | Structured `LITERATURE_SUPPORT` & `PATHWAY_ENRICHMENT` nodes |
| **AC-6** | Prior-guided analysis is visibly flagged in claim metadata and report summaries. | `test_knowledge.py`, `test_orchestrator.py` | `epistemic_tag="prior-guided"`, Markdown warning callout present |
| **AC-7** | In silico perturbation capability models gene knockouts/state shifts and creates versioned simulation artifacts. | `test_perturbation.py` | Non-negative matrices, versioned `adata://.../perturbation/vN` |
| **AC-8** | Evidence nodes reflect causal confidence updates without violating observational language rules. | `test_evidence_graph.py`, `test_perturbation.py` | $C_{causal} \le 0.50$, causal verb ban enforced on observational claims |
| **AC-9** | Comprehensive unit and integration test suite in `tests/` passes with 100% success. | All test files in `tests/` | 100% pass (`pytest -v --basetemp=.pytest_temp`) |
| **AC-10**| Full end-to-end simulated study executing scRNA + Spatial + Agent Adapter + Knowledge + Perturbation succeeds and produces a traceable markdown report. | `test_end_to_end_study.py` | 18-step study executes, complete Markdown report with clickable provenance |

---

## 5. Implementation Defect Escalation Log

During baseline test execution, the following implementation issue was detected and escalated to the implementing engineers:

- **Defect ID**: `DEFECT-001`
- **Location**: `eacbp/orchestrator/loop.py:196` / `eacbp/auditor/statistical.py:38` / `eacbp/artifact/storage.py:108`
- **Symptom**: `FileNotFoundError: Artifact payload not found at .../normalized/v2.h5ad` during `StatisticalValidator` audit in `test_full_ad_mouse_study_pipeline`.
- **Root Cause**: In the orchestrator execution loop, intermediate artifacts produced by tasks need to be consistently persisted and resolved via `ArtifactRegistry.get()` before the independent statistical validator attempts to load their payload from disk.
- **Recommended Fix for Implementing Engineers**: Ensure that when `ScientificOrchestrator` executes normalization and integration tasks, the output artifacts are registered and saved to `storage` before invoking `auditor.audit_task()`.
