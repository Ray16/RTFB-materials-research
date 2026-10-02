#!/bin/bash
# orop_worker.sh <GPU> — claim-based, resumable one-GPU worker for the OROP experimental-
# validation sweep (re-optimize both redox states in SMD(MeCN) at our production level, add
# GFN2 thermal, reference to our Fc; compare to OROP experimental potentials).
#
# Shared-node safety uses the SYSTEM-WIDE reservation gate (~/bin/gpu_reserve): reserves the
# GPU cluster-wide before touching the queue (any resident process = busy) and releases on exit.
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve
cd "$REPO" || exit 1
GPU="${1:?need GPU index}"
FLEET=logs/fleet/orop
CLAIMS="$FLEET/claims"
ITEMS="$FLEET/items.txt"
mkdir -p "$CLAIMS" "$FLEET"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export PYTHONPATH=src

CUR=""
cleanup(){ [ -n "$CUR" ] && rmdir "$CLAIMS/$CUR" 2>/dev/null; "$GR" release "$GPU" 2>/dev/null; }
trap cleanup EXIT INT TERM

# PREFLIGHT — env must run on this node's GPU (envs are node-local) or we reserve nothing.
if ! CUDA_VISIBLE_DEVICES="$GPU" "$PY" -c "import redox, gpu4pyscf, cupy; cupy.cuda.runtime.getDeviceCount()" >/dev/null 2>&1; then
  echo "$(hostname) gpu$GPU PREFLIGHT FAIL"; exit 3
fi
# SYSTEM-WIDE RESERVATION — refuses if the GPU has any resident process or is reserved elsewhere.
# Retry a few times: on a busy shared node a GPU can flicker in/out of use between the free-list
# snapshot and this acquire; a transient busy should not permanently kill the worker.
reserved=0
for attempt in 1 2 3 4 5 6; do
  if "$GR" acquire "$GPU" --pid $$ --label orop >/dev/null 2>&1; then reserved=1; break; fi
  sleep 10
done
if [ "$reserved" -ne 1 ]; then
  echo "$(hostname) gpu$GPU BUSY/RESERVED after retries -> no-op"; exit 6
fi

ran=0; fails=0
while read -r s; do
  [ -z "$s" ] && continue
  # done when BOTH redox states are cached
  if [ -f "calcs/orop/$s/ox/result.json" ] && [ -f "calcs/orop/$s/red/result.json" ]; then continue; fi
  mkdir "$CLAIMS/$s" 2>/dev/null || continue          # claim (skip if another worker has it)
  CUR="$s"
  log="$FLEET/sys_${s}.log"
  CUDA_VISIBLE_DEVICES="$GPU" "$PY" scripts/validation/orop/run_orop_benchmark.py --only "$s" --backend gpu > "$log" 2>&1
  if [ -f "calcs/orop/$s/ox/result.json" ] && [ -f "calcs/orop/$s/red/result.json" ]; then
    ran=$((ran+1)); fails=0
  else
    rmdir "$CLAIMS/$s" 2>/dev/null; fails=$((fails+1))   # release claim so it can be retried
    [ "$fails" -ge 4 ] && { echo "$(hostname) gpu$GPU CIRCUIT-BREAKER (ran $ran)"; exit 4; }
  fi
  CUR=""
done < "$ITEMS"
echo "$(hostname) gpu$GPU finished, ran $ran"
