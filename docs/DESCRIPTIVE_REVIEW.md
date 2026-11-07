# Descriptive evidence review

EACBP can now review an existing small-table scRNA export without importing an
h5ad or repeating the original upstream work. `review-descriptive` is for two
libraries whose independent biological replication is unverified or absent.
It makes a useful, auditable descriptive comparison possible without turning
cells or filtering rules into biological replicates.

## What changed and why

The Endometrium_MU research history contains finished counting, annotation and
molecular analysis, while its README still described an earlier waiting stage.
The usable completed results are the run log and hash-bound exports. Its main
epithelial population contains 464 reference and 51 query cells. Five existing
expression states have reference/query counts 163/16, 93/5, 76/11, 66/19 and 66/0.
None reaches the project's 20-cell display threshold in both libraries.

Native EACBP h5ad execution and native evidence verification answer different
questions: neither can honestly stand in for replaying externally frozen module
scores from these small exports. This addition is a separate application service
and CLI; it does not weaken native donor, statistical-audit or evidence-admission
gates. No remote jobs or original research files are changed.

The review explicitly separates:

- **Source integrity:** bundle pin, each source's SHA256, module definitions,
  original scoring manifest and fixed-background sensitivity declarations agree.
  Hashes establish consistency with supplied records, not an external signature.
- **Raw recomputation:** selected exported genes' UMI sums, detected cells,
  detection fractions and CPM are rebuilt for every supplied mask and expression
  state. CPM uses exported full-feature `total_counts`, not the selected panel's
  gene sum. The denominator's agreement with the full matrix is not rechecked.
- **External summary consistency:** group cell counts, finite score counts,
  comparison eligibility and query-minus-reference arithmetic are checked.
  Module scores, gene mapping and background controls are not recomputed.
- **Dependent sensitivity:** masks are deduplicated by cell IDs, not by matching
  cell counts. The subsets remain overlapping observations from the same input.
  A whole-population subtype view does not provide independent confirmation of
  the pooled view.

The report always shows sample sizes, the explicit replication status, state
coverage, mask equivalence, module ranges and selected-gene values. Ranges across
filters are not confidence intervals. Identity-inclusion sensitivity and gene
ablation require additional inputs and remain outside the current review.

## Commands

```powershell
python -m eacbp review-descriptive --bundle bundle.json --bundle-sha256 <SHA256> --output review.json --report report.md
```

Outputs must be distinct new files. Changed pins, malformed tables, missing
required rows, duplicate cell identities, incorrect counts/deltas and inconsistent
mask metadata fail with a nonzero exit code before a passing receipt is written.
The review loads each pinned file once and parses the same checked bytes. Each
input is limited to 64 MiB. No network access or optional scientific packages are
needed beyond EACBP's core dependencies.

`review.json` contains `native_evidence_admitted: false`,
`inference_status: "not_estimable"`, exact input pins, reviewer source hashes,
recomputed rows and the explicit verification scope. A passing receipt means
the stated checks passed; it does not certify animal sampling, cell identity,
upstream counting or a biological mechanism.

## Export contract

The supported format is `frozen_scoring_tables_v1`. The JSON bundle has
`schema_version: 1`, `study_id`, `population`, `reference_sample`, `query_sample`,
positive integer `min_cells`, and
`design: {"observation_unit": "library", "biological_replication": "unverified"}`.
`unreplicated` is also accepted. `sources` contains exactly these roles, each
with `path` and a lowercase `sha256`; relative paths resolve beside the bundle:

| Role | Contents |
| --- | --- |
| `cells` | Unique `cell_id`, explicit `sample_id`, `cluster`, full-feature `total_counts`, selected `GENE_raw_UMI` columns, paired `mask_NAME_pass` / `mask_NAME_metadata_known` booleans |
| `gene_summaries` | Pooled/state group, mask, sample and gene keys; retained cell count, library UMI, raw count sum, detected count, detection fraction, CPM |
| `module_summaries` | Group/mask/module keys; reference/query labels, retained and finite-score counts, medians, mask status and missing metadata counts, comparison status and delta |
| `scoring` | Completed descriptive scoring manifest, cell count, module definition hash/version and original normalization/score settings |
| `sensitivity` | Completed descriptive manifest binding the scoring hash, module hash, unchanged scores, contrast, minimum cells and score basis |
| `modules` | Module-set ID/version and unique module definitions |

This version supports the existing export naming convention:
`pooled / ALL_REVIEWED_EPITHELIAL`, `subcluster_state / subcluster_<cluster>`,
and an optional whole-population `subtype` matching the declared `population`.
It is not an arbitrary CSV importer. The package checker and tests define the
full required field and missing-value contract.

## Real Endometrium replay

```powershell
python scripts/verify_endometrium_descriptive.py --project F:/Endometrium_MU --run-dir outputs/endometrium_descriptive/new_run
```

The case adapter reads historical pins from
`results/epithelial_molecular_report_v1/validation.json`, verifies the delivered
tables, binds module definitions through the original scoring manifest, and
copies only six small input files into a new run. It saves the original index,
portable relative-path bundle, review, report and acceptance receipt. It does
not simply trust freshly computed hashes for changed research data.

The real exported-data replay verified 515 cells, 16 selected genes and 1,344
raw summary rows; it checked 19 modules and 931 module summary rows. Seven mask
names resolve to five distinct cell sets. The 11,760 summary rows outside the
selected panel are source-pinned but not independently reaggregated.

The observed p53/NF-κB/redox and pro-intrinsic-apoptosis score directions remain
descriptive library differences. The reviewer reports all modules, including
ones with changing signs, without selecting only favorable results. This does
not establish DNA breaks, activated signaling proteins, executed apoptosis,
mitochondrial dysfunction or a cell-autonomous knockout effect.

The source project has no completed biological sample table. Independent-animal
identities, pooling and experimental conditions still need experimental records
before inferential analysis can be supported. See the primary study on
[false discoveries from ignoring biological replication](https://www.nature.com/articles/s41467-021-25960-2).

## Regression checks

`tests/test_descriptive_review.py` covers small synthetic tables, incorrect
counts, full-library CPM denominators, exact member-set deduplication, missing
metadata, table omissions, duplicate identities and unsafe summary comparisons.
File/CLI tests cover pins, provenance binding and output preservation. The module
is in the core CI suite; the private real dataset is not needed by CI.
