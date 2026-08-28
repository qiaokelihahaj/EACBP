#!/usr/bin/env bash

# CPU-only generic Slurm template for EACBP.
# The gpu partition is used when the account is associated with the cluster GPU partition.
#SBATCH --job-name=eacbp-benchmark
#SBATCH --account=libin
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/public/home/qiaoke/eacbp_project/current/logs/slurm-%j.out
#SBATCH --error=/public/home/qiaoke/eacbp_project/current/logs/slurm-%j.err

set -euo pipefail

PROJECT_ROOT="${REMOTE_PROJECT_ROOT:-/public/home/qiaoke/eacbp_project/current}"
PYTHON="${REMOTE_PYTHON:-/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python}"

cd "$PROJECT_ROOT"
[[ -x "$PYTHON" ]] || {
  echo "Missing environment Python: $PYTHON" >&2
  exit 2
}

export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "=== [1/2] Verifying Environment ==="
"$PYTHON" scripts/verify_environment.py

echo "=== [2/2] Running Test Suite ==="
"$PYTHON" -m pytest -q
