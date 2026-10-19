#!/usr/bin/env bash
# Compatibility: detached single Python command with common resource preflight.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
export PYTHONPATH="$ROOT/src"
gpu=${1:?need GPU index}; log=${2:?need log path}; shift 2
mkdir -p "$(dirname "$log")"
run="$ROOT/logs/runs/command_$(date +%Y%m%d_%H%M%S)_$$"
setsid nohup "$PY" -m redox.workflows.gpu --gpu "$gpu" --out "$run" \
  --threads "${THREADS:-2}" -- "$PY" "$@" > "$log" 2>&1 < /dev/null &
echo "Resource controller PID $!; execution evidence and command log: $run"
