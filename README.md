# EACBP

Evidence-aware computational biology workflows in Python.

EACBP combines explicit task contracts, versioned local artifacts, independent audits, and evidence-linked Markdown reports. Real-data orchestration uses the standard library-backed methods; explicit demos use lightweight baselines. The knowledge sources remain local curated fixtures. This is a research workflow, not a validated clinical pipeline.

## Architecture

Unified evidence admission, execution resource policy, shared application services,
capability contracts, and portable research bundles are described in
[ARCHITECTURE_P4.md](ARCHITECTURE_P4.md).

Current code-based architecture map, review findings, and improvement priorities (Chinese): [ARCHITECTURE_REVIEW.md](ARCHITECTURE_REVIEW.md).

Runtime configuration separation, task execution helpers, and durable artifact audit access are described in [ARCHITECTURE_OPTIMIZATION.md](ARCHITECTURE_OPTIMIZATION.md). The earlier architecture review is historical; its resolved findings are marked at the top.

Capability declarations, multi-target branches, execution events, and conservative transaction cleanup are described in [ARCHITECTURE_P2.md](ARCHITECTURE_P2.md).

The installed `eacbp` CLI (`plan`, `run`, `resume`, `inspect`, `cleanup`, `report`) and verified evidence snapshots are described in [ARCHITECTURE_P3.md](ARCHITECTURE_P3.md). `python -m eacbp` provides the same commands from a source environment. The `verify-evidence` command verifies a SHA-bound P4 acceptance index and returns engineering, scientific and evidence-version results without rerunning analysis; see [P5 numerical storage and evidence verification](docs/ARCHITECTURE_P5_NUMERIC_EVIDENCE.md) for its supported envelope and historical-source archive option.

`eacbp review-descriptive` reconciles pinned external single-cell summary exports
when biological replication is unverified or absent. It recomputes selected-gene
UMI/detection/CPM from a small cell table, checks frozen module summaries, and
deduplicates QC masks by actual cell membership. It generates a report with the
study-design boundary and never admits these exports as native inferential
evidence. See [descriptive review and Endometrium validation](docs/DESCRIPTIVE_REVIEW.md).

`StudyManifest → dependency-ordered TaskContracts → CapabilityRegistry → Artifacts → ScientificAuditor → EvidenceGraph → ClaimEngine → report`

- `schemas/`: study, task, artifact and evidence contracts.
- `orchestrator/`: conservative intent parsing, explicit dependencies, routing, audit gates, task journals and resume.
- `capabilities/`: quantification, QC, normalization, clustering, donor/cell-level statistics, root-distance ordering, spatial analysis and simulations.
- `artifact/`: immutable version files, persistent metadata/lineage and read-time SHA-256 validation.
- `adapters/`: local SpaCell/ChatCell/GeneAgent-style computations with contract checks; no external agent service is connected.
- `knowledge/`: local curated literature and biological dictionaries. Records are not independently verified online.
- `auditor/`, `evidence/`, `report/`: checks, admitted evidence, constrained claims and provenance.

## Installation and tests

Use Python 3.12 for the full optional test install below; the `fate`/CellRank extra requires Python 3.12 or newer.

```sh
python -m pip install ".[dev,standard,advanced-statistics,advanced-qc,communication,fate]"
python -m pytest -q
```

The core CI matrix installs only `.[dev]` and runs the explicitly listed base suite on Windows and Linux with Python 3.11/3.12/3.13. The shared suite manifest is `.github/ci-test-suites.json`; every test module must be assigned to core, a capability lane, or the full integration lane. `environment.yml` is the legacy Conda core/bio setup for Python 3.10/3.11 and does not install the advanced capability extras.

Separate Linux/Python 3.12 jobs install and run the `standard`, `advanced-statistics`, `advanced-qc`, `communication`, and `fate` capability suites. Each selected capability test must run without skips. A final Linux/Python 3.12 integration job installs the exact scientific environment from `docker/requirements-linux-py312-scientific.lock` with hashes, runs the complete test tree, checks installed dependencies, and records `pip freeze` and JUnit output. It also builds with `docker/requirements-linux-py312-build.lock` and checks the installed wheel's CLI, web UI and acceptance assets, and Slurm batch template in a fresh virtual environment. The wheel smoke output is saved with the JUnit report and dependency snapshot. Its only accepted skips are the Windows-only process-supervision case on Linux and the test that needs an external CellBender smoke report. These checks exercise package installation and tests, not real FASTQ alignment.

## Browser workbench

Run `python -m eacbp webui` (or `eacbp webui` after installation) to open the
local browser interface. It supports h5ad file selection, study configuration,
plan previews, background execution, logs, run history, strict resume, and
saved-report viewing/revalidation. No extra frontend build or web-framework
dependency is required. Use `--workspace PATH`, `--runs-dir PATH`,
`--port 8765`, and `--no-browser` as needed. See [WEBUI.md](WEBUI.md) for usage,
SSH access and the single-user/single-job scope. Existing scientific methods
and audit gates are shared with the CLI.

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
`eacbp.application.study_service.run_study`. It keeps the historical flags and
`outputs/runs/<study-id>/<uuid>/` layout, while `run_config.json`, the artifact
registry, evidence snapshots, reports, and resume behavior come from the same
package lifecycle as `eacbp run`. Pytest uses its managed system temporary
directory by default. For isolated runs, pass `--basetemp` with a dedicated
scratch directory outside the repository; pytest clears that directory before
running, so it must never contain research inputs or saved results. Test source
and CI suite definitions stay in the repository; generated environments,
bytecode, build directories and test scratch data are disposable. Keep research
outputs and validation summaries separately from scratch data.

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

Container deployment design and implementation gates (Chinese): [CONTAINERIZATION_PLAN.md](CONTAINERIZATION_PLAN.md). The [P1 deployment guide](P1_DEPLOYMENT.md) provides a Python 3.12 CPU candidate Dockerfile and a one-study Slurm/Apptainer adapter using the existing CLI. Linux image builds and real cluster execution remain unvalidated.

P1 also adds versioned inference contracts, explicit assumption and scientific-result states, and neutral summaries of audited inconclusive results. See the [P1 acceptance record](docs/P1_ACCEPTANCE.md) and [research acceptance workflow](docs/P1_RESEARCH_ACCEPTANCE.md) for direct PyDESeq2 comparisons and the distinction between software regression and real-study validation.

### Advanced biological analysis

`--advanced` on `scripts/run_standard_study.py` selects actual PyDESeq2 donor-level DEG and donor leave-one-out sensitivity. Install `.[advanced-statistics]`. Donor counts, pairing and covariates must describe biological samples. Inadequate leave-one-out replication produces an explicit skip, not a robustness claim. Statistical null estimates remain missing, with a reason.

`--analysis-config config.json` accepts `capability_parameters`, `method_overrides`, `analysis_extensions`, and `advanced_analysis`. Unknown top-level options, unregistered capability names, and unknown target-override names are rejected. Plugin options belong under the registered capability's `capability_parameters` or `analysis_extensions` entry; arbitrary top-level plugin options must be migrated there. The plan preview lists otherwise valid capability configuration that is not applied to its provisional tasks in `unused_configuration` and its notices. For example:

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

A local browser interface is bundled through `eacbp webui`; multi-user execution and live PubMed/NCBI connectors are not implemented. Current simulations, marker references and statistical models require domain-specific validation on real data. Very large unsupported workloads fail explicitly instead of allocating unbounded dense matrices.

`PROJECT.md`, `TEST_INFRA.md`, and `TEST_READY.md` contain historical planning material; this README and executable tests describe current behavior. See `REPAIR_NOTES.md` for the corrective changes and remaining validation limits.

Public PBMC3k and Kang IFN-response reproduction records, fixed inputs on F:, and replay commands are in [PUBLIC_DATA_VALIDATION.md](docs/PUBLIC_DATA_VALIDATION.md). Paired-donor sensitivity analyses and direct-reference comparisons are in [RESEARCH_VALUE_VALIDATION.md](docs/RESEARCH_VALUE_VALIDATION.md); numerical agreement and known-response benchmarks do not establish new biological discoveries.

`eacbp analyze` accepts a scientific question, an H5AD file, species/tissue metadata and resource limits. It selects a bounded deterministic policy from observed counts and experimental design, records its decisions, quarantines prior analysis outputs, and runs the shared audited study lifecycle. It does not accept marker answers, a workflow configuration, manually selected populations or a trajectory root. Existing population labels are used only when the question explicitly requests observed populations. Every eligible population receives a donor-based contrast and donor sensitivity analysis; ambiguous or inadequate designs remain explicit gaps.

```powershell
eacbp analyze --question "Explore the cell populations and their characteristic genes" --data F:/data/input.h5ad --species human --tissue PBMC --run-dir F:/results/exploration
```

Install `.[standard,advanced-statistics]` for count-based exploration and qualified donor comparisons; add `velocity` for real stochastic RNA velocity on validated spliced/unspliced counts. scVelo 0.3.4 currently requires NumPy below 2.4: a real integration test exposed an upstream scalar-assignment incompatibility with NumPy 2.5. The velocity extra and runtime check enforce this bound. The supplementary Linux/CPython 3.12 lock in `docker/requirements-linux-py312-velocity.lock` overrides the base scientific lock's NumPy; install it after the build lock with `--no-build-isolation`. The base Docker image does not include this optional velocity environment.

Native `analysis_request_status` records whether all requested analyses actually completed and passed audit, separately from engineering `status`. This status is preserved in snapshots, reports and portable exports. It does not certify external biological benchmark agreement. `scripts/run_classic_autonomy_curriculum.py` replays three question-only exercises using the installed image; `scripts/evaluate_classic_autonomy.py` independently checks frozen reference expectations after execution. The evaluator never supplies analysis settings or expected answers to EACBP. Data, dependencies, caches and outputs for this local curriculum live under `F:/EACBP-autonomous-practice`.

New autonomous runs freeze a question/source-bound requirement contract before choosing qualified targets. Reports list each requested dimension, population, contrast and admitted task; unsupported or non-estimable requirements remain partial, while an audited estimate without significant features can answer a question. Coverage validates the original computation receipt and audit identity, including result inputs and scientific status. Historical runs without this contract retain their original completion policy. A completed, fully admitted compute journal may publish a partial scientific question through presentation recovery without executing new analyses.

Biological donor roles are shared across native analyses and independent audits. Sample, library and generic replicate identifiers are not automatic donor aliases. Explicit role declarations record their source and limitations, and cannot split an already observed individual into supposed biological replicates. Spatial analyses require measured coordinates; ligand–receptor analysis also requires real resource matches. Neighborhood graphs and communication permutations respect observed section/donor partitions. Analytical Moran/Geary inference currently supports one measured coordinate system; unsupported multi-section inference fails explicitly. These identity checks do not independently establish random sampling or causal identification.

If presentation saving fails after every scientific task has completed and passed audit, `eacbp finalize-completed-journal` can rebuild the snapshot and report from the original journal. It requires the original compute image digest and a hash-pinned installed-source proof, verifies the original task signatures, full plan, durable admissions and payload hashes, and records the current presentation version separately. It does not execute models or issue new scientific audits; the failure summary and original compute records remain available. Ordinary resume retains its strict source/environment constraints.

The question-only curriculum, independent outcomes, failed runs and resulting architecture improvements are recorded in [AUTONOMOUS_CLASSIC_REPRODUCTION.md](docs/AUTONOMOUS_CLASSIC_REPRODUCTION.md). Native QC records a disposition for every input cell inside its registered artifact and independently audits the complete filtering mask. `analyze` now freezes the primary 3 scaled MAD policy and two QC perturbations (2 and 4 scaled MAD) before execution, then reruns eligible branches through the native study lifecycle. A separate `qc_sensitivity` directory records full cell dispositions, shared-cell clustering agreement, differential-result sensitivity and native run hashes. The primary snapshot and completion status remain unchanged; failed, unavailable and duplicate branches remain explicit. Preserve the entire sensitivity directory with the primary run, since the primary study export does not include the separate branch registries. Low RNA remains a threshold flag, not independently proven cell damage, and no reference-scored branch replaces the default QC policy.

New planned PyDESeq2 tasks register per-fit JSON diagnostics beside their scientific tables, covering primary fits, LOO baselines and every attempted donor fold. Independent auditing checks the counts, design, contrast and fold coverage. Actual dispersion trends and captured Python warnings remain diagnostics; fit completion does not certify convergence. Native DEG and functional-activity evidence now binds a unique source row, verified at extraction, snapshot and report boundaries. Supporting and contradicting evidence receive the same trace and metric checks; legacy and summary evidence explicitly lack row verification. See [ARCHITECTURE_P1_EXPERIMENTS.md](docs/ARCHITECTURE_P1_EXPERIMENTS.md) for experiment scope, retained failures and verification records.

New `analyze` runs also preserve excluded cells through the shared `qc_review` branch. It inherits the actual QC thresholds and independently audits four outputs: a summary, the complete raw-count candidate matrix, all cell dispositions, and raw feature detection descriptions. Candidates remain separate from primary inference and generate no scientific claims. Reports expose these files and hashes only after durable admission. Configured studies can enable `analysis_extensions.qc_review: true`; legacy studies keep their existing opt-in behavior. See [ARCHITECTURE_P2_QC_CANDIDATES.md](docs/ARCHITECTURE_P2_QC_CANDIDATES.md) for PBMC regression, new Baron human pancreas experiments and retained scientific limits.

Native Leiden studies now record effective representation, neighbor and partition settings; independent auditing rebuilds the graph and partition, and reports read these settings after durable admission. The standalone `eacbp.application.clustering_sensitivity_execution.run_clustering_sensitivity` API freezes resolution, seed, raw-cell subsampling and eligible donor leave-one-out perturbations before native execution. It retains fixed primary QC thresholds and reports shared-cell ARI/NMI, per-cluster coverage and matched Jaccard without selecting a winner. See [ARCHITECTURE_P3_CLUSTERING_STABILITY.md](docs/ARCHITECTURE_P3_CLUSTERING_STABILITY.md) for the frozen protocol, execution evidence and scientific limits.

Prepared clustering experiments now persist separate invocation, native execution/admission, comparison and presentation attempts. Call `execute_clustering_sensitivity(prepared, resume=True)` to re-admit and reuse completed native studies, continue unstarted branches and recover failed presentation. Incomplete native studies remain preserved for explicit study-level recovery. Prior completed evidence hashes and the frozen source/environment must match. Each invocation retains durable process-tree samples and elapsed-stage records; sampled RSS sums are observations, not a physical-memory peak or resource limit. Resume writes a new comparison/report version and preserves prior attempts. See [ARCHITECTURE_P4_RECOVERY_CONFIRMATION.md](docs/ARCHITECTURE_P4_RECOVERY_CONFIRMATION.md).

Use `eacbp recovery-status --prepared F:/results/prepared --bindings-sha256 <independently-checked-sha256>` for a read-only recovery diagnosis, then `eacbp resume-experiment` with the same arguments for explicit continuation. The browser's run-history recovery card uses this shared gate and requires an existing prepared directory inside its selected results root. A ready diagnosis does not authorize reuse: execution must re-admit each existing native study. An incomplete native study still needs `eacbp resume --run-dir` separately. Experiment completion does not establish biological stability or external benchmark agreement. See [P6 recovery entry points](docs/ARCHITECTURE_P6_RECOVERY_ENTRYPOINTS.md).

New TABLE artifacts use a versioned `.table.npz` contract with explicit values, dtypes, missing-value representation and ordered axes. Legacy CSV artifacts retain their original decoding behavior. Durable audit hashes cover the entire contract and payload. Browser queries read typed tables, and CSV downloads produce a separate export; the evidence bundle retains the authoritative typed artifact. Unsupported table types fail explicitly. See [P5 numerical storage and evidence verification](docs/ARCHITECTURE_P5_NUMERIC_EVIDENCE.md).

The public experiment inputs and outputs from P0–P5 were deleted at the user's request on 2026-10-04. Code, reports and archived verification receipts remain. Their recorded checks are historical evidence; full payload replay requires obtaining the original data again. The cleanup receipt is `F:/EACBP-cleanup-2026-10-04/cleanup-verification-v2.json`. P6 uses a new small synthetic fixture on F: to validate recovery behavior.
