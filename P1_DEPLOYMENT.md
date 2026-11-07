# P1 deployment: Python 3.12 CPU image and one-study Slurm jobs

This slice adds a Linux/amd64 CPU container and a thin Slurm adapter. The image installs the existing `eacbp` wheel with all five scientific runtime extras: `standard`, `fate`, `advanced-statistics`, `advanced-qc`, and `communication`. Its checked-in lock contains 97 exact third-party distributions with PyPI hashes for Python 3.12 on Linux x86_64; the separate build lock pins the wheel-builder tools. The `dev` extra is not part of the runtime image, and the redundant `bio` extra is covered by `standard`. The image includes no GPU or CUDA framework. It runs the existing `eacbp run`, `eacbp resume`, or `eacbp report` command. Each submitted job performs one of those operations for one run directory. The adapter talks directly to Slurm; it does not keep a queue or database.

The canonical Slurm batch template stays at `slurm/run_eacbp_container.sbatch` in the source tree and is included in wheels as `share/eacbp/slurm/run_eacbp_container.sbatch`. The adapter locates that installed data file by default, with a source-tree fallback; `--script` selects a site-specific template. This keeps one template source for both repository and wheel installs.

## Build candidate image

The Dockerfile pins the official Python 3.12.14 Bookworm image to OCI index digest `sha256:392307d22300de8b5986851a12d9176dfc0fc073e65bf6523ebd7dcbeb23564e` ([official Docker Hub tag](https://hub.docker.com/_/python/tags?name=3.12.14-slim-bookworm)). The dependency ranges in `pyproject.toml` are resolved into `docker/requirements-linux-py312-scientific.lock` and installed with hash checking; wheel-build tools are pinned in `docker/requirements-linux-py312-build.lock`. The base digest is a multi-platform index, while this lock and deployment target are Linux/amd64. Regenerate the lock files on Linux/amd64 with Python 3.12 so environment markers match the target. The image build checks the scientific imports and `pip check`, then records the installed distributions in `/opt/eacbp-build/pip-freeze.txt` inside the image.

Build from a clean Linux checkout with Docker BuildKit:

```bash
docker build \
  --platform linux/amd64 \
  --file docker/Dockerfile \
  --build-arg EACBP_SOURCE_REVISION="$(git rev-parse HEAD)" \
  --tag eacbp:py312-cpu .

docker run --rm eacbp:py312-cpu cat /opt/eacbp-build/pip-freeze.txt \
  > requirements-linux-py312-analysis.freeze

docker run --rm eacbp:py312-cpu
```

The locks were resolved with uv 0.11.8 against `https://pypi.org/simple` on 2026-09-27. Regenerate them from the repository root with:

```bash
uv pip compile pyproject.toml \
  --extra standard --extra fate --extra advanced-statistics \
  --extra advanced-qc --extra communication \
  --no-emit-package eacbp --python-version 3.12 --generate-hashes \
  --index-url https://pypi.org/simple \
  --output-file docker/requirements-linux-py312-scientific.lock

uv pip compile docker/requirements-linux-py312-build.in \
  --python-version 3.12 --generate-hashes \
  --index-url https://pypi.org/simple \
  --output-file docker/requirements-linux-py312-build.lock
```

Review the new pins and rebuild the image when updating either lock.

The image uses a non-root default user, CPU thread defaults of one, and writable `/outputs` and `/scratch` locations. Docker callers can override the user to match their host account and thread count. The image does not include FASTQ aligners, local reference files, input data, or private configuration. `.dockerignore` excludes the repository's local data and build/test outputs; the Dockerfile copies only package source and build metadata.

The checked-in default now pins the base image and the Linux scientific dependencies. The image still needs a successful build and runtime smoke check on a host that can fetch the pinned base before it is treated as a verified container release. Do not use the Windows dependency snapshot as the Linux lock.

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

The wheel was built and installed into a temporary virtual environment, then imported from outside the source tree; the default batch template resolved from the installed `share/eacbp/slurm` data path. The deployment tests check argument construction and local validation. Git Bash syntax checks passed for the entrypoint and Slurm template. After configuring the user-level Docker Desktop proxy, the pinned image built with its 97-package hashed Linux/Python 3.12 lock, and `pip check` passed. The rebuilt image uses the non-root `10001:10001` account; a network-disabled startup probe verified its cache and temporary directories are writable. The synthetic full-container lifecycle check passed run, resume, report, export, import, and restored report with networking disabled and a read-only root filesystem; see [DOCKER_ACCEPTANCE.md](docs/DOCKER_ACCEPTANCE.md) for the proxy recovery, image digest, and latest status. `sbatch` and Apptainer remain unavailable, so the SIF, cluster bind policy, compute-node storage, and Slurm submission remain unverified.
