#!/usr/bin/env bash
# One study operation per Slurm job. Resources and log paths are normally
# supplied by `python -m eacbp.deployment submit`; these are safe defaults for
# direct `sbatch` use. Do not add a site-specific account or partition here.
#SBATCH --job-name=eacbp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

fail() {
  echo "eacbp slurm: $*" >&2
  exit 2
}

require_absolute_path() {
  local name="$1"
  local value="$2"
  [[ "$value" == /* ]] || fail "$name must be an absolute path"
  [[ "$value" != *:* ]] || fail "$name cannot contain ':' because it is used in an Apptainer bind"
}

resolve_directory() {
  local name="$1"
  local value="$2"
  require_absolute_path "$name" "$value"
  [[ -d "$value" ]] || fail "$name is not a directory: $value"
  realpath -e -- "$value"
}

resolve_file() {
  local name="$1"
  local root="$2"
  local relative="$3"
  [[ -n "$relative" && "$relative" != /* ]] || fail "$name must be a non-empty relative file path"
  case "/$relative/" in
    */../*) fail "$name cannot escape its input directory" ;;
  esac
  local resolved
  resolved="$(realpath -e -- "$root/$relative")" || fail "$name does not exist under $root: $relative"
  [[ "$resolved" == "$root/"* && -f "$resolved" ]] || fail "$name must resolve to a file inside $root"
  printf '%s\n' "$resolved"
}

: "${EACBP_ACTION:?Set EACBP_ACTION to run, resume, or report}"
: "${EACBP_RUN_ID:?Set EACBP_RUN_ID to one run-directory name}"
: "${EACBP_IMAGE:?Set EACBP_IMAGE to a prepared .sif path}"
: "${EACBP_OUTPUTS:?Set EACBP_OUTPUTS to an existing writable output directory}"

case "$EACBP_ACTION" in
  run|resume|report) ;;
  *) fail "EACBP_ACTION must be run, resume, or report" ;;
esac
case "$EACBP_RUN_ID" in
  ''|.|..|*[!A-Za-z0-9._-]*) fail "EACBP_RUN_ID must be one safe path component" ;;
esac

[[ "$EACBP_IMAGE" == /* && "$EACBP_IMAGE" == *.sif ]] || fail "EACBP_IMAGE must be an absolute .sif path"
[[ -r "$EACBP_IMAGE" && -f "$EACBP_IMAGE" ]] || fail "Container image is not readable: $EACBP_IMAGE"
EACBP_IMAGE="$(realpath -e -- "$EACBP_IMAGE")"
OUTPUTS="$(resolve_directory EACBP_OUTPUTS "$EACBP_OUTPUTS")"
RUN_DIR="$OUTPUTS/$EACBP_RUN_ID"
if [[ "$EACBP_ACTION" == run ]]; then
  [[ ! -e "$RUN_DIR" ]] || fail "run directory already exists; choose a new run id: $RUN_DIR"
else
  [[ -d "$RUN_DIR" ]] || fail "$EACBP_ACTION requires an existing run directory: $RUN_DIR"
fi
[[ -w "$OUTPUTS" ]] || fail "Output directory is not writable: $OUTPUTS"

if [[ -n "${EACBP_IMAGE_SHA256:-}" ]]; then
  [[ "$EACBP_IMAGE_SHA256" =~ ^[[:xdigit:]]{64}$ ]] || fail "EACBP_IMAGE_SHA256 must contain 64 hexadecimal characters"
  actual_sha="$(sha256sum -- "$EACBP_IMAGE" | cut -d ' ' -f 1)"
  [[ "${actual_sha,,}" == "${EACBP_IMAGE_SHA256,,}" ]] || fail "Container SIF SHA-256 does not match EACBP_IMAGE_SHA256"
fi

action_args=()
bind_args=(--bind "$OUTPUTS:/outputs")
if [[ "$EACBP_ACTION" == run ]]; then
  : "${EACBP_DATA_DIR:?Set EACBP_DATA_DIR to the read-only input directory}"
  : "${EACBP_CONFIG_DIR:?Set EACBP_CONFIG_DIR to the read-only config directory}"
  : "${EACBP_DATA:?Set EACBP_DATA to a file relative to EACBP_DATA_DIR}"
  : "${EACBP_MANIFEST:?Set EACBP_MANIFEST to a file relative to EACBP_CONFIG_DIR}"
  DATA_ROOT="$(resolve_directory EACBP_DATA_DIR "$EACBP_DATA_DIR")"
  CONFIG_ROOT="$(resolve_directory EACBP_CONFIG_DIR "$EACBP_CONFIG_DIR")"
  DATA_FILE="$(resolve_file EACBP_DATA "$DATA_ROOT" "$EACBP_DATA")"
  MANIFEST_FILE="$(resolve_file EACBP_MANIFEST "$CONFIG_ROOT" "$EACBP_MANIFEST")"
  DATA_RELATIVE="${DATA_FILE#"$DATA_ROOT"/}"
  MANIFEST_RELATIVE="${MANIFEST_FILE#"$CONFIG_ROOT"/}"
  bind_args+=(--bind "$DATA_ROOT:/data:ro" --bind "$CONFIG_ROOT:/config:ro")
  action_args+=(--manifest "/config/$MANIFEST_RELATIVE" --data "/data/$DATA_RELATIVE" --run-dir "/outputs/$EACBP_RUN_ID")
  if [[ -n "${EACBP_CONFIG:-}" ]]; then
    CONFIG_FILE="$(resolve_file EACBP_CONFIG "$CONFIG_ROOT" "$EACBP_CONFIG")"
    CONFIG_RELATIVE="${CONFIG_FILE#"$CONFIG_ROOT"/}"
    action_args+=(--config "/config/$CONFIG_RELATIVE")
  fi
  if [[ -n "${EACBP_REFS_DIR:-}" ]]; then
    REFS_ROOT="$(resolve_directory EACBP_REFS_DIR "$EACBP_REFS_DIR")"
    bind_args+=(--bind "$REFS_ROOT:/refs:ro")
  fi
else
  # Resume/report use the paths persisted in run_config.json. Re-bind any
  # external resources at their original container paths when the saved
  # configuration refers to /data, /config, or /refs.
  if [[ -n "${EACBP_DATA_DIR:-}" ]]; then
    DATA_ROOT="$(resolve_directory EACBP_DATA_DIR "$EACBP_DATA_DIR")"
    bind_args+=(--bind "$DATA_ROOT:/data:ro")
  fi
  if [[ -n "${EACBP_CONFIG_DIR:-}" ]]; then
    CONFIG_ROOT="$(resolve_directory EACBP_CONFIG_DIR "$EACBP_CONFIG_DIR")"
    bind_args+=(--bind "$CONFIG_ROOT:/config:ro")
  fi
  if [[ -n "${EACBP_REFS_DIR:-}" ]]; then
    REFS_ROOT="$(resolve_directory EACBP_REFS_DIR "$EACBP_REFS_DIR")"
    bind_args+=(--bind "$REFS_ROOT:/refs:ro")
  fi
  action_args+=(--run-dir "/outputs/$EACBP_RUN_ID")
fi

APPTAINER="${EACBP_APPTAINER:-apptainer}"
command -v "$APPTAINER" >/dev/null 2>&1 || fail "Apptainer executable is unavailable: $APPTAINER"
cpus="${SLURM_CPUS_PER_TASK:-4}"
[[ "$cpus" =~ ^[1-9][0-9]*$ ]] || fail "SLURM_CPUS_PER_TASK must be a positive integer"
SCRATCH="${EACBP_SCRATCH_DIR:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/eacbp-${SLURM_JOB_ID:-local}-$EACBP_RUN_ID}"
[[ "$SCRATCH" == /* ]] || fail "Scratch directory must be an absolute path: $SCRATCH"
require_absolute_path EACBP_SCRATCH_DIR "$SCRATCH"
mkdir -p "$SCRATCH"
SCRATCH="$(realpath -e -- "$SCRATCH")"
bind_args+=(--bind "$SCRATCH:/scratch")

export APPTAINERENV_HOME=/scratch/home
export APPTAINERENV_XDG_CACHE_HOME=/scratch/cache
export APPTAINERENV_NUMBA_CACHE_DIR=/scratch/cache/numba
export APPTAINERENV_MPLCONFIGDIR=/scratch/cache/matplotlib
export APPTAINERENV_TMPDIR=/scratch/tmp
export APPTAINERENV_CUDA_VISIBLE_DEVICES=""
export APPTAINERENV_OMP_NUM_THREADS="$cpus"
export APPTAINERENV_MKL_NUM_THREADS="$cpus"
export APPTAINERENV_OPENBLAS_NUM_THREADS="$cpus"
export APPTAINERENV_NUMBA_NUM_THREADS="$cpus"

mkdir -p "$SCRATCH/home" "$SCRATCH/cache/numba" "$SCRATCH/cache/matplotlib" "$SCRATCH/tmp"
echo "EACBP action=$EACBP_ACTION run_id=$EACBP_RUN_ID job_id=${SLURM_JOB_ID:-manual}"
echo "Container SIF: $EACBP_IMAGE"
echo "Scratch: $SCRATCH"

exec "$APPTAINER" exec --cleanenv --no-home --pwd /outputs "${bind_args[@]}" \
  "$EACBP_IMAGE" /opt/venv/bin/eacbp "$EACBP_ACTION" "${action_args[@]}"
