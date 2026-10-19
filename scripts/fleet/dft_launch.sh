#!/usr/bin/env bash
# Compatibility: global shard count, first shard, then local GPU indices.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
export PYTHONPATH="$ROOT/src"
N=${1:?need global shard count}; FIRST=${2:?need first shard}; shift 2
[ "$#" -gt 0 ] || { echo 'Need at least one GPU index' >&2; exit 2; }
printf -v selection '%s,' "$@"
selection="${selection%,}"
read -r -a extra <<< "${DFT_ARGS:-}"
exec "$PY" -m redox.workflows.launch --module dft --execute \
  --out "$ROOT/logs/runs/dft_shards_$(date +%Y%m%d_%H%M%S)_$$" \
  --gpus "$selection" --shard-count "$N" --shard-start "$FIRST" \
  --threads "${THREADS:-2}" -- "${extra[@]}"
