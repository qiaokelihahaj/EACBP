#!/usr/bin/env bash

# Slurm batch script for EACBP Autonomous Kat8 cKO Single-Cell Study.
#SBATCH --job-name=eacbp-kat8
#SBATCH --account=libin
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/public/home/qiaoke/eacbp_project/current/logs/slurm-kat8-%j.out
#SBATCH --error=/public/home/qiaoke/eacbp_project/current/logs/slurm-kat8-%j.err

set -euo pipefail

PROJECT_ROOT=""
PYTHON=""

cd ""
[[ -x "" ]] || {
  echo "Missing environment Python: " >&2
  exit 2
}

export PYTHONPATH=""
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=""
export MKL_NUM_THREADS=""
export OPENBLAS_NUM_THREADS=""

echo "=== [1/3] Verifying Environment ==="
"" scripts/verify_environment.py

echo "=== [2/3] Running Autonomous Kat8 cKO Study ==="
"" scripts/run_kat8_study.py

echo "=== [3/3] Study Complete. Report generated in outputs/reports/ ==="
