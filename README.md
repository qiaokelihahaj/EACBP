# EACBP

Evidence-aware computational biology workflows in Python.

EACBP combines explicit task contracts, versioned local artifacts, independent audits, and evidence-linked Markdown reports. Real-data orchestration uses the standard library-backed methods; explicit demos use lightweight baselines. The knowledge sources remain local curated fixtures. This is a research workflow, not a validated clinical pipeline.

## Architecture

Current code-based architecture map, review findings, and improvement priorities (Chinese): [ARCHITECTURE_REVIEW.md](ARCHITECTURE_REVIEW.md).

Runtime configuration separation, task execution helpers, and durable artifact audit access are described in [ARCHITECTURE_OPTIMIZATION.md](ARCHITECTURE_OPTIMIZATION.md). The earlier architecture review is historical; its resolved findings are marked at the top.

Capability declarations, multi-target branches, execution events, and conservative transaction cleanup are described in [ARCHITECTURE_P2.md](ARCHITECTURE_P2.md).

The installed `eacbp` CLI (`plan`, `run`, `resume`, `inspect`, `cleanup`, `report`) and verified evidence snapshots are described in [ARCHITECTURE_P3.md](ARCHITECTURE_P3.md). `python -m eacbp` provides the same commands from a source environment.

`StudyManifest → dependency-ordered TaskContracts → CapabilityRegistry → Artifacts → ScientificAuditor → EvidenceGraph → ClaimEngine → report`

- `schemas/`: study, task, artifact and evidence contracts.
- `orchestrator/`: conservative intent parsing, explicit dependencies, routing, audit gates, task journals and resume.
- `capabilities/`: quantification, QC, normalization, clustering, donor/cell-level statistics, root-distance ordering, spatial analysis and simulations.
- `artifact/`: immutable version files, persistent metadata/lineage and read-time SHA-256 validation.
- `adapters/`: local SpaCell/ChatCell/GeneAgent-style computations with contract checks; no external agent service is connected.
- `knowledge/`: local curated literature and biological dictionaries. Records are not independently verified online.
- `auditor/`, `evidence/`, `report/`: checks, admitted evidence, constrained claims and provenance.

## Installation and tests

```sh
python -m pip install ".[dev,standard,fate]"
python -m pytest -q
```

AnnData/Scanpy are required for full single-cell processing. CI is configured for tests and wheel builds on Windows/Linux with Python 3.11/3.12/3.13. Local test results do not establish that remote CI or real FASTQ alignment has run.

## Real data and demo data

FASTQ quantification defaults to `mode="real"`. Real mode requires valid sample metadata, every FASTQ pair, the selected executable (kb-python or STAR), and explicit references. Missing dependencies, invalid inputs and alignment failures return a structured failure. STARsolo preserves all lanes and reads its own run's filtered Matrix Market output; raw output must be explicitly requested with `use_raw_counts=True`. Missing filtered output never silently falls back to raw or a previous run.

Only explicit `mode="demo"` may generate synthetic data. Synthetic origin propagates through artifact lineage, evidence and report labels. Demo statistical outputs are demonstrations, not biological measurements.

Use the script help for the current command-line interface:

```sh
python scripts/run_fastq_to_biology_pipeline.py --help
python scripts/run_kat8_study.py --help
python scripts/run_standard_study.py --help
```

Script output is isolated by study/run directory. Starting another run never deletes an earlier artifact tree. Real sequencing inputs must preserve every lane and carry explicit condition/donor metadata; a filename alone is not an experimental design.

`scripts/run_standard_study.py` is a compatibility adapter around
`eacbp.cli.run_study`. It keeps the historical flags and
`outputs/runs/<study-id>/<uuid>/` layout, while `run_config.json`, the artifact
registry, evidence snapshots, reports, and resume behavior come from the same
package lifecycle as `eacbp run`. Tests that need isolated scratch output
should use a short, scoped `.pytest_temp_<topic>` directory (for example
`.pytest_temp_rf`) and `-p no:cacheprovider`; these directories
are ignored without hiding source or research `temp*` directories.

Descriptor callbacks use the documented signatures `validator(contract, result,
registry)` and `evidence_extractor(contract, result, report, registry)`; shorter
or keyword-only callbacks need an explicit wrapper. Parameter-model compatibility
with Pydantic, dataclasses, and callables remains available. `ScientificPolicy`
is deprecated and retained for existing imports; runtime policy lives in the
router and independent auditors.

The generic h5ad entry point has no default disease or cell type:

```sh
python scripts/run_standard_study.py --data study.h5ad --species homo_sapiens --tissue kidney --root-cell-id CELL_ID --paga
```

Omit `--root-cell-id` when no justified trajectory root is available; that branch is then omitted. Use `--condition-a/--condition-b` for an explicit contrast, `--marker-reference markers.json` for annotation, and `--terminal-states terminals.json` for CellRank. FASTQ CLI supports `--quant-tool starsolo_v1 --star-bin STAR --genome-dir INDEX --whitelist-path WHITELIST` (with `--sample-manifest`); `--gtf-path` is optional when the index already contains annotation.

## Implemented methods

| Capability | Actual local implementation |
|---|---|
| Normalization | Library-size normalization, log1p; preserves counts layer |
| Integration | PCA batch-mean centering or no correction |
| Clustering | Lloyd K-means and marker-score annotation |
| Display embedding | First two embedding coordinates; not UMAP |
| Differential expression | Donor count aggregation, library normalization and Welch tests; explicitly exploratory cell-level fallback |
| Trajectory-like ordering | Root-relative Euclidean embedding distance, reproducible subsampling and FDR-adjusted Spearman associations |
| Spatial analysis | Nearest-neighbor spatial graph, spatial autocorrelation and proximity-weighted ligand/receptor scores |
| Perturbation | Bounded in-silico network/compound models; no experimental causal confirmation |

The table above describes the explicit `baseline` profile. The `standard` profile calls actual third-party libraries:

| Method ID | Implementation |
|---|---|
| `harmonypy_v1` | Harmony soft clustering and batch correction using harmonypy 0.0.10 |
| `scanpy_leiden_umap_v1` | Scanpy neighbor graph, Leiden via leidenalg, and UMAP |
| `scanpy_dpt_v1` | Scanpy diffusion map and diffusion pseudotime; optional PAGA graph |
| `cellrank_fate_v1` | CellRank kernel and GPCCA estimator fate probabilities with supplied terminal cells |

For a single batch, the planner selects `no_correction_v1`. Standard clustering retains existing `cell_type` annotations or uses an explicit `marker_reference` dictionary; it never substitutes ground truth annotations. Without a reference it returns named clusters. DPT requires `root_cell_id`; without it, automatic trajectory analysis is omitted and reported. PAGA additionally requires an observed grouping. DPT ordering alone is not a validated trajectory stability result.

```python
summary = orchestrator.run_study(manifest, {
    "method_profile": "standard",
    "capability_parameters": {
        "clustering": {"marker_reference": {"RequestedType": ["GENE1", "GENE2"]}},
        "trajectory_inference": {"root_cell_id": "actual_cell_id", "run_paga": True},
    },
    # Optional: IDs must occur in the analyzed population.
    "cellrank_terminal_states": {"FateA": ["terminal_cell_a"], "FateB": ["terminal_cell_b"]},
})
```

CellRank uses DPT pseudotime in this pipeline. Its standalone capability can also accept an explicit row-stochastic `obsp['transition_matrix']` with `kernel='precomputed'`. It does not infer RNA velocity from raw spliced/unspliced counts. Fate probabilities are conditional on the supplied transition model and terminal states.

Historical identifiers such as `harmony`, `leiden_knn_v1`, `paga_dpt` and `cellrank` may be recognized as migration aliases. Returned method identifiers describe the actual algorithm. They do not connect those named third-party algorithms. Legacy Leiden/UMAP output fields are opt-in and labelled as aliases.

## Evidence and failure policy

- Failed execution/audits do not contribute evidence. Dependent tasks become blocked; independent branches can continue.
- Empty or missing support IDs cannot produce admitted claims. Statistical inference requires audited adjusted-significance evidence. Level 3 additionally requires verified knowledge support.
- Nonsignificant results and empty dynamic-gene tables are valid scientific outcomes, not a reason to fabricate findings.
- Claims restate specific results; disease/target keywords do not automatically generate marker, mechanism or plaque-location conclusions.
- Curated knowledge has no invented p values and supplies explicitly unverified contextual hypotheses.
- Confidence values are heuristic summaries, not calibrated probabilities.

## Resume and provenance

Reuse the same registry directory and pass `current_state={..., "resume": True}` to `run_study`, with the same manifest and settings. A successful saved task is reused only when its contract and input hashes match and its output files still pass integrity checks. Strict reproducibility also fingerprints package source and installed dependencies. Restored outputs are audited again. Each computation attempt writes into a private transaction directory. A successful attempt publishes all artifact metadata and its recovery receipt in one atomic index replacement. A crash before publication leaves no public outputs and the attempt can run again; a crash after publication reuses the receipt even if the study journal was not saved. Evidence is admitted only after the audited checkpoint is saved. A new invocation resets in-memory state and evidence.

Transaction directories contain the actual payloads of committed artifacts and must not be deleted wholesale. Uncommitted attempts remain on disk for diagnosis but are invisible through the public registry; automatic garbage collection is not implemented. Atomic publication applies to orchestrated task artifacts, not external programs' side effects or direct individual `register()` calls. Older runs that published partial outputs without transaction receipts still require a separate run directory.

Metadata audit adapts the remaining plan: missing or ambiguous contrasts omit automatic differential statistics and their consumers; insufficient donors omit abundance inference. Explicit invalid contrasts fail. With no requested cell type the workflow analyzes all cells; an absent requested cell type fails rather than silently selecting another population. Omitted branches and reasons appear in the report.

The task journal uses an OS process lock and atomic replacement. Artifact hashes detect payload changes relative to the local metadata index; this is not a tamper-proof external signature of the entire filesystem. An older artifact directory without an index is not silently reconstructed as trusted data.

## Scope and migration

### Advanced biological analysis

`--advanced` on `scripts/run_standard_study.py` selects actual PyDESeq2 donor-level DEG and donor leave-one-out sensitivity. Install `.[advanced-statistics]`. Donor counts, pairing and covariates must describe biological samples. Inadequate leave-one-out replication produces an explicit skip, not a robustness claim. Statistical null estimates remain missing, with a reason.

`--analysis-config config.json` accepts `capability_parameters`, `method_overrides`, `analysis_extensions`, and `advanced_analysis`. For example:

```json
{
  "advanced_analysis": true,
  "capability_parameters": {
    "deg": {"condition_a": "treated", "condition_b": "control", "donor_col": "donor", "paired": true}
  },
  "analysis_extensions": {
    "functional_activity": {
      "network_path": "C:/data/network.csv",
      "network_source": "Explicit curated resource",
      "network_version": "1",
      "tmin": 5
    }
  }
}
```

The local network needs `source`, `target`, and `weight` columns and must match the study species. Its contents enter the resume signature. Functional activity is inferred through decoupler and compared at the donor level, not a direct protein-activity measurement. Scrublet and local CellTypist have independent QC audits. LIANA supports donor-specific inference, audited condition comparisons, reports and resume. CellBender is available with explicit external inputs and an installed runtime. See [BIOINFORMATICS_IMPLEMENTATION.md](BIOINFORMATICS_IMPLEMENTATION.md) for verification evidence and limits.

Install `.[advanced-qc]` for Scrublet and CellTypist. In `analysis_extensions`, enable `doublet_detection` with a library `batch_key` (use `null` only for an explicitly single-library input), and enable `cell_annotation` with a local `model_path`. CellTypist keeps reference labels in `cell_type_celltypist`, records unknown/conflicting labels, and preserves existing labels by default. These QC observations do not add mechanism or causal support. For downstream communication using the reference labels, select the corresponding `cell_type_col` explicitly. See [ADVANCED_QC.md](ADVANCED_QC.md) for the label-conflict policy and parameters.

CellBender runs in the environment where the orchestrator runs; a Windows process cannot directly execute a Linux binary path. With EACBP and CellBender installed in a Linux/WSL runtime, an example configuration is:

```json
{
  "analysis_extensions": {
    "background_removal": {
      "unfiltered_input_path": "/data/raw_feature_bc_matrix.h5",
      "executable": "/opt/cellbender/bin/cellbender",
      "output_path": "/data/run/corrected.h5",
      "run_cwd": "/data/run",
      "extra_args": ["--expected-cells", "500", "--total-droplets-included", "2000", "--cpu-threads", "4"]
    }
  }
}
```

Choose the droplet/cell settings for the actual dataset, create the run directory first, and use a fresh output path. The explicit input must include empty droplets; the imported target artifact can contain a subset of the same cell/gene IDs. The planner runs background removal before QC and normalizes `corrected_counts` while retaining raw `counts`. Input, executable and any explicitly supplied existing `--checkpoint` file are hashed; the output destination is not hashed as an input. The independent audit re-reads and maps the actual external output. Report/log warnings remain visible, and CLI success does not certify biological correction quality. A real public-example run passed the adapter, full workflow and resume checks, while retaining ELBO convergence warnings.

To enable LIANA, include an explicit local resource in the analysis configuration:

```json
{
  "analysis_extensions": {
    "liana_communication": {
      "lr_resource_path": "/absolute/path/ligand_receptor.csv",
      "lr_resource_version": "your-resource-version",
      "lr_resource_source": "your-resource-citation",
      "cell_type_col": "reference_cell_type",
      "donor_col": "donor_id",
      "condition_col": "condition",
      "condition_a": "treated",
      "condition_b": "control",
      "min_cells": 5
    }
  }
}
```

The resource must contain `ligand` and `receptor` columns appropriate for the manifest species; no resource is downloaded automatically. Each donor-condition group needs at least two cell types meeting `min_cells`. Set `paired: true` for repeated donors. Condition comparisons use donor ranks with BH adjustment separately for each score; LIANA ranks themselves are not FDR. Interactions that cannot be evaluated because genes are absent are recorded with the missing genes in artifact metadata and task metrics, and checked independently. They are not evidence that communication is absent. These outputs describe inferred communication and cannot establish signaling or causality.

There is no bundled web UI, multi-user execution service or live PubMed/NCBI connector. Current simulations, marker references and statistical models require domain-specific validation on real data. Very large unsupported workloads fail explicitly instead of allocating unbounded dense matrices.

`PROJECT.md`, `TEST_INFRA.md`, and `TEST_READY.md` contain historical planning material; this README and executable tests describe current behavior. See `REPAIR_NOTES.md` for the corrective changes and remaining validation limits.
