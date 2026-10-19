#!/usr/bin/env bash
# Compatibility launcher: all suitable free GPUs, common preflight and reservations.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
export PYTHONPATH="$ROOT/src"
if [ "${BACKEND:-gpu}" != "gpu" ]; then
  echo 'For CPU work use redox.qm.dft --backend cpu with an explicit shared-node allocation.' >&2
  exit 2
fi
selection=()
if [ -n "${GPUS:-}" ]; then selection=(--gpus "$GPUS"); fi
exec "$PY" -m redox.workflows.launch --module dft --execute \
  --out "$ROOT/logs/runs/dft_$(date +%Y%m%d_%H%M%S)_$$" \
  --threads "${GPU_THREADS:-2}" "${selection[@]}" -- "$@"
