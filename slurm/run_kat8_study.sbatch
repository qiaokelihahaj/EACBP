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

PROJECT_ROOT="${REMOTE_PROJECT_ROOT:-/public/home/qiaoke/eacbp_project/current}"
PYTHON="${REMOTE_PYTHON:-/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python}"
: "${EACBP_H5AD:?Set EACBP_H5AD to the real input h5ad path}"

cd "$PROJECT_ROOT"
[[ -x "$PYTHON" ]] || { echo "Missing Python: $PYTHON" >&2; exit 2; }
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="$OMP_NUM_THREADS"
export OPENBLAS_NUM_THREADS="$OMP_NUM_THREADS"

"$PYTHON" scripts/verify_environment.py
"$PYTHON" scripts/run_kat8_study.py --mode real --data-path "$EACBP_H5AD" --output-dir "$PROJECT_ROOT/outputs"
echo "Study finished; reports are under outputs/runs/<study>/<run>/reports."
