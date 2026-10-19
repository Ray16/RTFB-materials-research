#!/usr/bin/env bash
# Compatibility launcher: common preflight and reservations for UMA warm starts.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
export PYTHONPATH="$ROOT/src"
selection=()
if [ -n "${GPUS:-}" ]; then selection=(--gpus "$GPUS"); fi
exec "$PY" -m redox.workflows.launch --module uma --execute \
  --out "$ROOT/logs/runs/uma_$(date +%Y%m%d_%H%M%S)_$$" \
  --threads "${THREADS:-2}" "${selection[@]}" -- --model "${MODEL:-uma-s-1p2p1}" "$@"
