# P1 implementation and acceptance record

P1 adds an explicit inference contract, reproducible method acceptance, and a
single-study container/Slurm deployment path. It does not establish biological
validity merely because software tests pass.

## Inference and evidence

Advanced donor-level analyses carry a versioned `InferenceContract` describing
the question, estimand, population, input artifacts, observation and independent
units, pairing, design, counts source, contrast, alpha, confidence level, FDR
family, and sensitivity plan. Its content determines its identity. The resolved
contract participates in task execution and resume signatures.

`ScientificResult` distinguishes supported estimates, inconclusive estimates,
not-estimable results, and failed assumptions. Assumption assessments distinguish
passed, failed, unknown, and not-assessable checks. A full-rank design does not
certify independence between biological donors.

Audited inconclusive results can produce an exact, neutral Level 1 result
description with zero discovery confidence. They cannot support a positive
statistical claim, absence-of-effect claim, or equivalence claim. Existing
significance and biological-interpretation gates remain in force. Failed or
unaudited results remain visible as limitations in the report.

Inference identities are carried through task results, artifact metadata,
evidence, claims, saved snapshots, and report provenance cards. Cross-reference
checks reject mismatched identities.

## Reproducible acceptance

See [P1_RESEARCH_ACCEPTANCE.md](P1_RESEARCH_ACCEPTANCE.md) for the manifest and
`python -m eacbp.acceptance` workflow. Input hashes, source/license, donor and
condition metadata, design, resource limits, and reference environment are
explicit. Method agreement compares effects, confidence intervals, and FDR
against a direct PyDESeq2 invocation. Engineering, method, and research status
are reported separately.

Synthetic fixtures are software regression evidence. A real-data run additionally
needs a specified dataset, its provenance and license, a fixed reference
environment, and scientific review of its design and interpretation.

## Deployment

See [P1_DEPLOYMENT.md](../P1_DEPLOYMENT.md) for the CPU image and Slurm/Apptainer
workflow. Each job invokes the existing `run`, `resume`, or `report` CLI for one
study. Scheduler controls do not introduce a second workflow engine.

The development host has no usable Linux Docker daemon, Slurm commands, or
Apptainer runtime. Image builds, container smoke execution, and actual scheduler
execution require a suitable host and remain deployment acceptance work.

## Validation record

Validation was performed in the existing Windows/Python 3.13 environment. A
synthetic CLI acceptance run completed with 24 tested features and no significant
features. Against direct PyDESeq2, maximum absolute differences were below
`1.2e-16` for effect, both Wald interval bounds, and FDR. The fixture emitted
dispersion-fit and small-residual-degree-of-freedom warnings; numerical agreement
does not remove those limitations or make it a real-study validation.

Both shell scripts passed Git Bash syntax checks. The wheel was built and
installed into a temporary environment, and the Slurm template resolved from its
installed data path when imported outside the repository. These checks do not
substitute for a Linux container build or an actual Slurm job.

Frozen integrated suite: **539 passed, 1 skipped** in 181.82 seconds. The skip is
the CellBender check that requires an external smoke report. The advanced
statistics CI runner completed **53 tests with zero skips** in 25.35 seconds.
`git diff --check` passed. The final wheel contains the acceptance manifest and
the exact canonical Slurm template.

```powershell
$env:CELLTYPIST_FOLDER = "$PWD\.cache\celltypist"
$env:MPLCONFIGDIR = "$PWD\.cache\matplotlib"
.\.venv\Scripts\python.exe -m pytest -q --basetemp=.pytest_temp_p1_full --disable-warnings
.\.venv\Scripts\python.exe scripts/verify_ci_capability_tests.py advanced-statistics
.\.venv\Scripts\python.exe -m pip wheel --no-deps --no-build-isolation --wheel-dir .cache/p1-wheel-final .
```

No real biological dataset or cluster configuration was supplied for this run.
Real-study acceptance, a clean Linux image build, container smoke execution, and
Slurm execution remain unperformed. No Git commit, push, image publication, or
cluster submission was made.
