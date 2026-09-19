# Advanced QC and annotation interfaces

The optional capabilities in `eacbp.capabilities.advanced_qc` call the named
scientific implementation directly.  They do not replace a missing package or
external executable with a baseline calculation.  A dependency failure is
returned by `CapabilityRegistry.execute_contract` as an execution failure.

## Capability IDs and outputs

| Capability | Implementation ID | Required input | Primary output |
| --- | --- | --- | --- |
| `doublet_detection` | `scanpy_scrublet_v1` | one AnnData artifact with raw integer-valued counts in `counts_layer` (default `counts`) or `X`, and a complete library-batch column (default `batch`) | AnnData, default `adata://{study_id}/doublet_detection/v1` |
| `cell_annotation` | `celltypist_local_v1` | one AnnData artifact plus an explicit local `model_path` | AnnData, default `adata://{study_id}/cell_annotation/v1` |
| `background_removal` | `cellbender_cli_v1` | one AnnData artifact plus explicit unfiltered droplet input and executable paths | AnnData, default `adata://{study_id}/background_removal/v1` |

When `TaskContract.expected_outputs` is present, its first URI is used exactly
and must be an AnnData URI in the same study.  An optional second URI may be a
`table://` per-cell/audit table or a `json://` run report.  More than two
expected outputs is rejected rather than silently omitted.

## Scrublet doublet scoring

`ScrubletDoubletCapability` (`ScanpyScrubletCapability` is an alias) runs
`scanpy.pp.scrublet` on a working AnnData copy whose `X` is the validated raw
counts matrix.  If `counts_layer` (default `counts`) is present it is used;
otherwise raw `X` is used.  Values must be finite, non-negative and integer
valued.  The caller can set `batch_key` (or `batch_col`) to the library batch;
the default requires `obs['batch']` and runs Scrublet independently by batch.
Set `batch_key=None` only for an explicitly single-batch contract.

The default `filter_doublets=False` writes `obs['doublet_score']` and
`obs['predicted_doublet']` and retains every cell.  Set
`filter_doublets=True` to request filtering explicitly.  The output keeps the
original `X`, all existing layers, and a `counts` layer (created from `X` when
the input had no counts layer).  The output `uns['eacbp_scrublet']` and
artifact metrics record cell counts before and after, marked and filtered
doublets, batch counts, requested and effective principal components, the
requested threshold, and the threshold found for each batch.  The Scrublet
payload itself remains in `uns['scrublet']`.

Relevant parameters include `counts_layer`, `batch_key`, `filter_doublets`,
`threshold`, `expected_doublet_rate`, `sim_doublet_ratio`, `n_prin_comps`,
`n_neighbors`, and `random_seed`.  `threshold=None` delegates threshold
selection to Scanpy and therefore requires the Scanpy Scrublet optional
dependencies needed by that version.

The executed operations are `validate_raw_counts`,
`score_doublets_by_batch`, `record_doublet_threshold`, and
`mark_predicted_doublets`; an explicit filter adds `filter_cells` so a task
contract can forbid or allow that mutation.

## Local CellTypist annotation

`CellTypistAnnotationCapability` (`CellTypistCapability` is an alias) requires
`parameters.model_path` to be an existing local model file.  URL/model-name
resolution and automatic download are not used.  The file is SHA-256 hashed
before loading.  If a planner supplies
`parameters.external_resource_sha256`, the value for `model_path` must match
the observed hash; a mismatch fails the task.  The hash, resolved path, and
`model_metadata` are copied into `uns['celltypist']` and task metrics.

Model species and tissue can be supplied as `model_metadata={'species': ...,
'tissue': ...}` or via `model_species`/`model_tissue`.  The implementation also
uses metadata exposed by the local model where available.  Missing values are
recorded as `"unknown"` and `metadata_complete=False`; they are never guessed.
If the input has species/tissue metadata that conflicts with the model, the
conflict is recorded and explicit downstream replacement is suppressed.

CellTypist labels are written to `obs['cell_type_celltypist']` by default.
Existing `obs['cell_type']` values are copied to `cell_type_existing`, retained
unchanged, and never silently overwritten.  Unknown labels (`Unknown`,
`Unassigned`, empty values, and equivalent tokens) and row-level conflicts are
marked in `cell_type_celltypist_unknown` and
`cell_type_celltypist_conflict`.  Set `use_as_cell_type=True` to request use of
the model labels as the downstream `cell_type`; unknown labels, metadata
conflicts, and row-level conflicts still retain the existing value.  Use
`annotation_key` and `existing_cell_type_key` to change these column names.
`use_as_cell_type=True` is an explicit downstream annotation choice: it writes
the CellTypist result to the selected `existing_cell_type_key` (default
`cell_type`) subject to the unknown/conflict retention rule.  It leaves a
clustering column such as `cluster` untouched unless the caller explicitly
sets `existing_cell_type_key="cluster"`; even then conflicting rows remain
unchanged and the separate CellTypist column is available for downstream
selection.

The optional table output contains one row per input cell, the chosen
CellTypist label, unknown/conflict flags, and any one-dimensional prediction
columns returned by CellTypist.  The executed operations are
`validate_local_annotation_model`, `celltypist_annotation`, and
`record_annotation_conflicts`.

## CellBender external CLI

`CellBenderBackgroundRemovalCapability` (`CellBenderCapability` is an alias)
requires all of the following parameters:

```python
{
    "unfiltered_input_path": "/run/raw_feature_bc_matrix.h5",
    "executable": "/opt/cellbender/bin/cellbender",
    "output_path": "/run/cellbender_output.h5",
    "run_cwd": "/run",                 # optional; defaults to current cwd
    "extra_args": ["--epochs", "150"], # optional
}
```

`unfiltered_input_path` must exist and is never replaced with a filtered
artifact.  `output_path` is required so the external side effect is explicit
and must end in `.h5`, the CellRanger v3 output format.  An existing output is
rejected unless `allow_overwrite=True` is explicit.  The adapter invokes
`[executable, "remove-background", "--input", input,
"--output", output, *extra_args]` without a shell.  It records the resolved
executable, executable/input/output hashes, command, current run directory,
return code, and bounded stdout/stderr in `uns['cellbender']` and metrics.

The standard `output.h5` contains all original droplet barcodes; the adapter
allows it to be larger than the target AnnData artifact.  It maps the external
matrix by unique cell and gene IDs into the target cell/gene order.  Every
target ID must be present exactly once; duplicate or missing IDs fail, while
extra external droplets/features are recorded in the alignment report.
Dimensions-only matches are insufficient.  The original matrix remains in
`X`; an existing `layers['counts']` is preserved (or is created from `X`), and
the ID-mapped external matrix is stored only in
`layers['corrected_counts']`.  No synthetic corrected matrix is produced when
the executable, output, reader, hash, or alignment check fails.  The optional
second output is a JSON run report.  Executed operations are
`validate_unfiltered_droplet_input`, `run_cellbender_cli`,
`validate_cell_gene_alignment`, and `store_corrected_counts_layer`.

The adapter's provenance record identifies the matrix source as
`cellbender_cli_output` and records the resolved run directory, command,
external input/executable/output paths and SHA-256 hashes, mapping dimensions,
extra droplet/feature counts, and the preserved layer names.  `extra_args`
cannot override the adapter's explicit `--input` or `--output` arguments.

Standalone downstream normalization selects `input_layer='corrected_counts'`
explicitly. The advanced planner selects this layer when background removal
is requested. The background-removal capability itself preserves `X` and raw
`counts`; corrected expression does not replace the original measurements.
Donor count models still use their explicit `counts_layer` (default `counts`).
To fit them on corrected counts, set that parameter to `corrected_counts` in
the corresponding analysis contract; the integer-count and design checks
still apply.

## Audit requirements

Contracts should permit the operations listed above and should include
validation requirements for raw-count integrity, batch-specific thresholds,
model hash/metadata, or external CLI alignment as applicable.  The independent
`AdvancedQCValidator` re-hashes the CellBender input, executable and output,
checks contract paths and expected external hashes, reloads the external `.h5`
matrix, repeats the unique-ID mapping, and compares it with
`corrected_counts`.  A contract that forbids `filter_cells` will make an
explicit `filter_doublets=True` result fail the side-effect guard.  Failures do
not register a primary output, so downstream tasks cannot consume a synthetic
or partially written result.

CellBender 0.4.0 has been run in an isolated WSL Python 3.12 environment with
CPU PyTorch, using the official tiny mouse-heart example (37,760 droplets,
100 genes). The actual CLI produced a matrix, report, metrics and log. The
project adapter mapped a reordered 12-cell, 8-gene subset, preserved raw
counts, and passed independent structural auditing. Runtime records are in
`outputs/cellbender_smoke`.

The example report contains training/test ELBO warnings and says the output
could be suboptimal. These remain quality warnings; successful execution and
structural auditing do not certify correction quality. The adapter retains
visible report warnings and sidecar hashes. It ignores HTML styles/scripts,
and does not treat ordinary explanatory text about convergence as a warning.
`scripts/verify_cellbender_runtime.py` records runtime integrity separately
from scientific quality and returns a nonzero quality-caveat status when
warnings remain. This small example is integration evidence, not validation
of results for a user's dataset.
