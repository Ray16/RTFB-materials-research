#!/bin/bash
# d3level_fleet_worker.sh <GPU> — claim-based, resumable one-GPU worker for the D3TaLES-level reorg
# validation (lc_wpbe + per-molecule tuned omega + def2-svp, gas).
#
# Shared-node safety uses the SYSTEM-WIDE reservation gate (~/bin/gpu_reserve, NOT a per-repo lock):
# reserves the GPU cluster-wide before touching the queue (a GPU with ANY resident process — another
# user's OR one of our own other sessions — is busy), and releases on exit. See CLAUDE.md.
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve
cd "$REPO" || exit 1
GPU="${1:?need GPU index}"
FLEET=logs/fleet/d3level
CLAIMS="$FLEET/claims"
CALC=results/d3tales_reorg_validation_d3level/calc
ITEMS="$FLEET/items.txt"
mkdir -p "$CLAIMS" "$FLEET"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export PYTHONPATH=src

CUR_MOL=""
cleanup(){ [ -n "$CUR_MOL" ] && rmdir "$CLAIMS/$CUR_MOL" 2>/dev/null; "$GR" release "$GPU" 2>/dev/null; }
trap cleanup EXIT INT TERM

# PREFLIGHT — env must run on this node's GPU (envs are node-local) or we reserve nothing.
if ! CUDA_VISIBLE_DEVICES="$GPU" "$PY" -c "import redox, gpu4pyscf, cupy, rdkit; cupy.cuda.runtime.getDeviceCount()" >/dev/null 2>&1; then
  echo "$(hostname) gpu$GPU PREFLIGHT FAIL"; exit 3
fi
# SYSTEM-WIDE RESERVATION — refuses if the GPU has any resident process or is reserved elsewhere.
# Owner = this worker's pid (long-lived) so the reservation is held for the worker's whole life.
if ! "$GR" acquire "$GPU" --pid $$ --label d3level >/dev/null 2>&1; then
  echo "$(hostname) gpu$GPU BUSY/RESERVED -> no-op"; exit 6
fi

ran=0; fails=0
while read -r id; do
  [ -z "$id" ] && continue
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then continue; fi
  mkdir "$CLAIMS/$id" 2>/dev/null || continue
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then rmdir "$CLAIMS/$id" 2>/dev/null; continue; fi
  CUR_MOL="$id"
  log="$FLEET/mol_${id}.log"
  CUDA_VISIBLE_DEVICES="$GPU" "$PY" scripts/validate_reorg_worker_d3tales.py --only "$id" --backend gpu > "$log" 2>&1
  if grep -q "DONE_MARKER" "$log" 2>/dev/null; then
    ran=$((ran+1)); fails=0
  else
    rmdir "$CLAIMS/$id" 2>/dev/null; fails=$((fails+1))
    [ "$fails" -ge 3 ] && { echo "$(hostname) gpu$GPU CIRCUIT-BREAKER (ran $ran)"; exit 4; }
  fi
  CUR_MOL=""
done < "$ITEMS"
echo "$(hostname) gpu$GPU finished, ran $ran"
