#!/usr/bin/env bash

# Slurm batch script for EACBP Unified FASTQ-to-Biology Autonomous Pipeline.
#SBATCH --job-name=eacbp-fastq
#SBATCH --account=libin
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/public/home/qiaoke/eacbp_project/current/logs/slurm-fastq-%j.out
#SBATCH --error=/public/home/qiaoke/eacbp_project/current/logs/slurm-fastq-%j.err

set -euo pipefail

PROJECT_ROOT=""
PYTHON=""
FASTQ_DIR="/public/home/qiaoke/data/P9kat8pr单细胞数据/Data/Data/CleanData"

cd ""
[[ -x "" ]] || {
  echo "Missing environment Python: " >&2
  exit 2
}

mkdir -p "/logs" "/outputs/reports"

export PYTHONPATH=""
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=""
export MKL_NUM_THREADS=""
export OPENBLAS_NUM_THREADS=""

echo "=== [1/3] Verifying Environment ==="
"" scripts/verify_environment.py

echo "=== [2/3] Running Full FASTQ-to-Biology Autonomous Pipeline ==="
"" scripts/run_fastq_to_biology_pipeline.py --fastq-dir "" --target-gene "Kat8"

echo "=== [3/3] Execution Finished. Scientific report ready in outputs/reports/ ==="
