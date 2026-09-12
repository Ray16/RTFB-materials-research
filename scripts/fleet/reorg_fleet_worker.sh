#!/bin/bash
# reorg_fleet_worker.sh <GPU_INDEX> — claim-based, resumable, self-healing one-GPU worker for the
# Track-2 D3TaLES reorganization-energy sweep (our production SMD-opt pipeline, electron couple).
# Multi-node safe: claims are atomic `mkdir` on NFS; released on failure so another node retries.
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python   # NFS env, visible on all nodes
cd "$REPO" || exit 1
GPU="${1:?need GPU index}"

FLEET=logs/fleet/reorg
CLAIMS="$FLEET/claims"
CALC=results/d3tales_reorg_validation_prod/calc
ITEMS_FILE="$FLEET/items.txt"
mkdir -p "$CLAIMS" "$FLEET"

# CPU-thread cap so N GPU-workers per node never oversubscribe (GPU does the heavy lifting).
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 XTB_THREADS=2
export PYTHONPATH=src

# (1) PREFLIGHT — a node that cannot run must NOT claim anything.
if ! CUDA_VISIBLE_DEVICES="$GPU" "$PY" -c "import redox.dft, gpu4pyscf, cupy; cupy.cuda.runtime.getDeviceCount()" >/dev/null 2>&1; then
  echo "$(hostname) gpu$GPU PREFLIGHT FAIL -> no-op"; exit 3
fi

ran=0; fails=0
while read -r id; do
  [ -z "$id" ] && continue
  # (2) DONE-TEST: skip a molecule whose result is already ok/partial.
  # skip any molecule already ATTEMPTED (ok/partial/error) — retrying deterministic failures
  # (OOM, SMD-gradient non-convergence) just circuit-breaks workers. Errors are triaged offline.
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then
    continue
  fi
  mkdir "$CLAIMS/$id" 2>/dev/null || continue          # atomic claim; skip if taken
  # re-check after claiming (another node may have finished it between test and claim)
  # skip any molecule already ATTEMPTED (ok/partial/error) — retrying deterministic failures
  # (OOM, SMD-gradient non-convergence) just circuit-breaks workers. Errors are triaged offline.
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then
    continue
  fi
  log="$FLEET/mol_${id}.log"
  CUDA_VISIBLE_DEVICES="$GPU" "$PY" scripts/validate_reorg_worker_prod.py --only "$id" --backend gpu \
      > "$log" 2>&1
  if grep -q "DONE_MARKER" "$log" 2>/dev/null; then
    ran=$((ran+1)); fails=0                             # keep the claim
  else
    rmdir "$CLAIMS/$id" 2>/dev/null; fails=$((fails+1)) # release for another node
    if [ "$fails" -ge 3 ]; then
      echo "$(hostname) gpu$GPU CIRCUIT-BREAKER after $fails consecutive fails (ran $ran)"; exit 4
    fi
  fi
done < "$ITEMS_FILE"
echo "$(hostname) gpu$GPU finished, ran $ran"
