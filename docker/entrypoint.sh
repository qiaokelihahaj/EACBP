#!/bin/sh
set -eu

mkdir -p \
  "${HOME:-/scratch/home}" \
  "${XDG_CACHE_HOME:-/scratch/cache}" \
  "${NUMBA_CACHE_DIR:-/scratch/cache/numba}" \
  "${MPLCONFIGDIR:-/scratch/cache/matplotlib}" \
  "${TMPDIR:-/scratch/tmp}"

exec "$@"
