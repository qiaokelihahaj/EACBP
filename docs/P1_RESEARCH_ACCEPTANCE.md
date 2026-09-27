# P1 research acceptance package

The acceptance runner validates a fixed h5ad input and study design, executes the existing EACBP PyDESeq2 pseudobulk capability, and compares it with a separate direct call to the installed PyDESeq2 API. Results are JSON with separate `engineering`, `method`, and `research_acceptance` statuses.

This workspace has no supplied real h5ad dataset, source/license record, or verified environment lockfile. The included manifest is a fill-in template. Synthetic test success demonstrates software behavior only and is never reported as biological evidence.

## Run

Copy `eacbp/acceptance/fixtures/research_manifest.template.json` and fill in the real data path, SHA256, accession/citation, license, species, counts layer, contrast, donor/pairing, batch/covariates, resource limits, exact installed package versions, and hashed lockfile. For example, calculate file digests in PowerShell with:

```powershell
Get-FileHash -Algorithm SHA256 .\dataset.h5ad
Get-FileHash -Algorithm SHA256 .\environment.lock
```

Run the package directly without changing the main EACBP CLI:

```powershell
.\.venv\Scripts\python.exe -m eacbp.acceptance .\acceptance_manifest.json --output .\acceptance_result.json
```

The process exits 0 when the input/design checks and the direct method comparison pass. A research candidate still reports `pending_domain_review`: numerical agreement cannot establish biological validity or justify interpretation. Synthetic manifests report `not_applicable_synthetic` and `biological_validity_claimed: false`.

## Input and design checks

The runner checks the h5ad's SHA256 before reading it; requires the requested raw counts layer, condition and donor columns, and declared batch column; validates counts through the existing PyDESeq2 preparation path; checks donor qualification and pairing; and asks PyDESeq2 to construct the design matrix before fitting. A rank-deficient or zero-residual-degree-of-freedom design is rejected. Optional resource limits cover cells, genes, nonzero counts, input bytes, and the existing pseudobulk working-memory estimate.

If `species_metadata_key` points to a matching `AnnData.uns` value, the species assertion is checked against the file. Otherwise species is explicitly recorded as a manifest attestation. Batch is recorded in the manifest and is included in the model when `batch_column` is set.

## Direct reference comparison

EACBP is run through `PyDESeq2PseudobulkCapability`. The reference independently groups raw cells into donor-condition count rows and directly calls `pydeseq2.dds.DeseqDataSet` and `pydeseq2.ds.DeseqStats` with the same design, contrast, alpha, Cook's filter, and independent-filter settings. The JSON records the runtime versions, fit configuration, reference comparison tolerances, finite comparisons, missing-value counts, and maximum absolute differences for effect, interval bounds, and FDR.

PyDESeq2 returns the log2 fold change and its standard error; the package compares Wald normal intervals computed from those values with the intervals recorded by EACBP. A sampled process RSS is recorded when `psutil` is available. The sampling interval and its peak-missing limitation are included in the output. The result also records cell, gene, nonzero-count, donor-condition group, and batch counts.

## Synthetic regression checks

Run only the P1 acceptance regression module with the repository's supported Windows interpreter and an explicit isolated pytest temp directory:

```powershell
.\.venv\Scripts\python.exe -m pytest tests/test_research_acceptance.py --basetemp=.pytest-basetemp-p1-research-20260927 -q
```

The regression module covers a non-significant paired design that must complete successfully, direct reference agreement, one-donor rejection, complete condition/batch confounding rejection, and recovery through the existing `ScientificOrchestrator` without recomputing or admitting duplicate evidence. These are software/method regressions; they are not real-study acceptance.

The workspace currently has PyDESeq2 0.5.4 in its Windows/Python 3.13.7 environment. That observation is emitted by a local run but is not a portable research lock. The template intentionally leaves the version lock unfilled; pin a clean target environment and hash its actual lockfile before claiming a fixed-version research run.
