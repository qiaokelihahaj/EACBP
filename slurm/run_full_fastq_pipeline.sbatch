#!/usr/bin/env bash

# Slurm batch script for EACBP Unified FASTQ-to-Biology Autonomous Pipeline.
#SBATCH --job-name=eacbp-fastq
#SBATCH --account=libin
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/public/home/qiaoke/eacbp_project/current/logs/slurm-fastq-%j.out
#SBATCH --error=/public/home/qiaoke/eacbp_project/current/logs/slurm-fastq-%j.err

set -euo pipefail

PROJECT_ROOT="${REMOTE_PROJECT_ROOT:-/public/home/qiaoke/eacbp_project/current}"
PYTHON="${REMOTE_PYTHON:-/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python}"
FASTQ_DIR="/public/home/qiaoke/data"

cd "$PROJECT_ROOT"
[[ -x "$PYTHON" ]] || {
  echo "Missing environment Python: $PYTHON" >&2
  exit 2
}

mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/outputs/reports"

export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export PATH="/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin:$PATH"

echo "=== [1/3] Verifying Environment ==="
"$PYTHON" scripts/verify_environment.py

echo "=== [2/3] Submitting Full FASTQ-to-Biology Autonomous Pipeline ==="
"$PYTHON" scripts/run_fastq_to_biology_pipeline.py \
  --fastq-dir "$FASTQ_DIR" \
  --target-gene "Kat8" \
  --output-dir "$PROJECT_ROOT/outputs"

echo "=== [3/3] Execution Finished. Scientific report ready in outputs/reports/ ==="
