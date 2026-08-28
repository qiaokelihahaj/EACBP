#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="${REMOTE_PROJECT_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PYTHON="${REMOTE_PYTHON:-/public/home/qiaoke/.local/share/mamba/envs/perturb-seq/bin/python}"

cd "$PROJECT_ROOT"
[[ -x "$PYTHON" ]] || {
  echo "Missing Python interpreter: $PYTHON" >&2
  exit 2
}

export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"
"$PYTHON" scripts/verify_environment.py