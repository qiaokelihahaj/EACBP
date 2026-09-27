# P1 deployment: Python 3.12 CPU image and one-study Slurm jobs

This slice adds a candidate Linux container and a thin Slurm adapter. The container installs the existing `eacbp` wheel with its `standard`, `fate`, and `advanced-statistics` extras, including PyDESeq2 and decoupler. It then runs the existing `eacbp run`, `eacbp resume`, or `eacbp report` command. Each submitted job performs one of those operations for one run directory. The adapter talks directly to Slurm; it does not keep a queue or database. Advanced QC and communication extras remain outside this image target.

The canonical Slurm batch template stays at `slurm/run_eacbp_container.sbatch` in the source tree and is included in wheels as `share/eacbp/slurm/run_eacbp_container.sbatch`. The adapter locates that installed data file by default, with a source-tree fallback; `--script` selects a site-specific template. This keeps one template source for both repository and wheel installs.

## Build candidate image

Build from a clean Linux checkout with Docker BuildKit. `python:3.12-slim-bookworm` and the dependency ranges in `pyproject.toml` are candidate inputs, not a verified lock. The image build checks imports and `pip check`, then records the exact resolved Python distributions in `/opt/eacbp-build/pip-freeze.txt` inside that image. Save that file with the image after a successful Linux build, review the base image digest and dependency freeze, and only then treat the build as a deployment candidate.

```bash
docker build \
  --file docker/Dockerfile \
  --build-arg EACBP_SOURCE_REVISION="$(git rev-parse HEAD)" \
  --tag eacbp:py312-cpu .

docker run --rm eacbp:py312-cpu cat /opt/eacbp-build/pip-freeze.txt \
  > requirements-linux-py312-analysis.freeze
```

The image uses a non-root default user, CPU thread defaults of one, and writable `/outputs` and `/scratch` locations. Docker callers can override the user to match their host account and thread count. The image does not include FASTQ aligners, local reference files, input data, or private configuration. `.dockerignore` excludes the repository's local data and build/test outputs; the Dockerfile copies only package source and build metadata.

This candidate does not pin the base image by digest or contain a checked-in Linux lock. Do not publish or describe it as reproducible until a clean Linux build has supplied and reviewed both. Do not use the Windows dependency snapshot as the Linux lock.

## Prepare the cluster

Prepare a SIF from a reviewed OCI image digest on a supported build host, then stage the SIF on storage visible to the compute nodes. Compute jobs do not pull images. Record the SIF SHA-256 and pass it as `EACBP_IMAGE_SHA256`. The Slurm template uses Apptainer, does not request a GPU, disables the host home bind, and runs with `--cleanenv`. The cluster's Apptainer version and administrator bind policy still need validation on a compute node.

Create persistent directories before submission. The output root must support the locking, hard-link, and atomic-rename behavior required by EACBP; validate those filesystem properties on the target mount before production use. The output root and log directory must be visible to the relevant compute nodes. Choose a scratch path on node-local or shared storage according to the site's policy.

For a new run, export the following values in the shell from which you will submit. All host paths are absolute. `EACBP_DATA`, `EACBP_MANIFEST`, and optional `EACBP_CONFIG` are file paths relative to their corresponding input directory. Paths in the manifest/config that refer to inputs or reference resources must use the container paths `/data`, `/config`, or `/refs`.

```bash
export EACBP_ACTION=run
export EACBP_RUN_ID=study-2026-09-27
export EACBP_IMAGE=/shared/containers/eacbp-analysis.sif
export EACBP_IMAGE_SHA256=<64-hex-SIF-sha256>
export EACBP_OUTPUTS=/shared/eacbp/outputs
export EACBP_DATA_DIR=/shared/eacbp/data
export EACBP_DATA=study.h5ad
export EACBP_CONFIG_DIR=/shared/eacbp/config
export EACBP_MANIFEST=manifest.json
export EACBP_CONFIG=analysis.json  # optional; omit or unset to use CLI defaults
export EACBP_REFS_DIR=/shared/eacbp/refs  # optional
```

Submit one Slurm job for that operation. `--partition` and `--account` are site-specific and intentionally have no defaults. The helper creates the log directory before calling `sbatch` and sets Slurm output/error paths to `eacbp-%j.out` and `eacbp-%j.err` there.

```bash
python -m eacbp.deployment submit \
  --logs-dir /shared/eacbp/logs \
  --cpus-per-task 4 --mem 16G --time 04:00:00 \
  --partition <site-cpu-partition> --account <site-account>
```

Omit the site-specific flags if the cluster provides suitable defaults. `submit --dry-run` validates the image, paths, action, and run directory, prints the planned `sbatch` arguments, and does not create the log directory or contact Slurm. The command passes environment variables to Slurm; keep that environment limited to non-secret job settings and paths.

To resume an interrupted/failed run, wait until its earlier Slurm job is stopped, then set `EACBP_ACTION=resume` and retain the same `EACBP_RUN_ID`, SIF, and output root before submitting again. For a report rebuild, set `EACBP_ACTION=report`. Resume does not reopen the original h5ad, manifest, or config; those values are stored in the run directory. It does reuse absolute external resource paths stored from the original run, so re-bind any required resource roots at the same container locations: set `EACBP_DATA_DIR`, `EACBP_CONFIG_DIR`, or `EACBP_REFS_DIR` again if the saved configuration refers to `/data`, `/config`, or `/refs`. Resume uses EACBP's existing integrity checks and can reject a run if the environment or saved inputs no longer match. A killed job is not guaranteed to be resumable at every interruption point; keep the run directory and inspect it before retrying.

## Manage a submitted job

The helper uses only Slurm's `sbatch`, `squeue`, `sacct`, and `scancel` commands. It accepts numeric job IDs and array task IDs; it does not accept shell fragments.

```bash
python -m eacbp.deployment status 123456
python -m eacbp.deployment logs 123456 --logs-dir /shared/eacbp/logs --lines 200
python -m eacbp.deployment cancel 123456
```

`status` checks `squeue`, then falls back to accounting with `sacct` after the job leaves the active queue. `logs` prints the last lines of stdout and stderr. `cancel` requests `scancel`; check the final state and logs before submitting a resume job. Do not delete the complete run directory or its `_transactions` directory as scratch cleanup. Only remove a node-local scratch directory after the job is finished and its contents are no longer needed.

## Scope and verification status

The batch script has no project, account, partition, or absolute home-directory defaults. It checks the action and run ID, confines input file paths to their read-only roots, rejects a pre-existing run directory for `run`, and requires one for `resume`/`report`. The Python Slurm adapter constructs argv without a shell, validates scheduler identifiers and resource strings, supports a no-submit dry run, and wraps status, cancel, and log viewing.

The wheel was built and installed into a temporary virtual environment, then imported from outside the source tree; the default batch template resolved from the installed `share/eacbp/slurm` data path. The deployment tests check argument construction and local validation. Git Bash syntax checks passed for the entrypoint and Slurm template. These checks do not validate a SIF, Apptainer behavior, cluster storage, or scientific analysis in the container.

The current Windows session has a Docker client but no reachable Docker engine; the WSL status query was denied, and `sbatch`/Apptainer are unavailable locally. Therefore this change has **not** been built or run as a Linux image and has not been submitted to a Slurm cluster.
