#!/bin/bash
# d3level_fleet_worker.sh <GPU> — claim-based, resumable, AUTO-YIELDING one-GPU worker for the
# D3TaLES-level reorg validation (lc_wpbe + per-molecule tuned omega + def2-svp, gas).
#
# Shared-node enforcement (see scripts/fleet/gpu_guard.sh, cluster.env, CLAUDE.md "Submitting GPU jobs"):
#   * ONE TASK PER GPU: atomic cluster-wide gpu_claim; a second worker on the same GPU no-ops.
#   * NEVER CO-RESIDENT: refuses to start on a GPU another user is on, and a background sentinel
#     kills the in-flight calc + yields the GPU if a foreign job appears mid-run (molecule retried).
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
cd "$REPO" || exit 1
GPU="${1:?need GPU index}"
# shellcheck disable=SC1091
source scripts/fleet/cluster.env 2>/dev/null || true
# shellcheck disable=SC1091
source scripts/fleet/gpu_guard.sh
FLEET=logs/fleet/d3level
CLAIMS="$FLEET/claims"
CALC=results/d3tales_reorg_validation_d3level/calc
ITEMS="$FLEET/items.txt"
mkdir -p "$CLAIMS" "$FLEET"
export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
export PYTHONPATH=src

CUR_MOL=""                                   # molecule claim held right now (for crash cleanup)
cleanup(){ [ -n "$CUR_MOL" ] && rmdir "$CLAIMS/$CUR_MOL" 2>/dev/null; gpu_release "$GPU"; }
trap cleanup EXIT INT TERM

# PREFLIGHT — env must run on this node's GPU (envs are node-local) or we claim nothing.
if ! CUDA_VISIBLE_DEVICES="$GPU" "$PY" -c "import redox, gpu4pyscf, cupy, rdkit; cupy.cuda.runtime.getDeviceCount()" >/dev/null 2>&1; then
  echo "$(hostname) gpu$GPU PREFLIGHT FAIL"; exit 3
fi
# ENFORCE one-task-per-GPU + no other user resident, BEFORE touching the queue.
if ! gpu_claim "$GPU"; then echo "$(hostname) gpu$GPU ALREADY-CLAIMED -> no-op"; exit 6; fi
if ! gpu_free_of_others "$GPU"; then echo "$(hostname) gpu$GPU OCCUPIED-BY-OTHER -> yield"; exit 7; fi

ran=0; fails=0
while read -r id; do
  [ -z "$id" ] && continue
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then continue; fi
  mkdir "$CLAIMS/$id" 2>/dev/null || continue
  if [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null; then rmdir "$CLAIMS/$id" 2>/dev/null; continue; fi
  CUR_MOL="$id"
  log="$FLEET/mol_${id}.log"
  CUDA_VISIBLE_DEVICES="$GPU" "$PY" scripts/validate_reorg_worker_d3tales.py --only "$id" --backend gpu > "$log" 2>&1 &
  PYPID=$!
  yielded=0
  while kill -0 "$PYPID" 2>/dev/null; do
    sleep "${SENTINEL_INTERVAL:-20}"
    if ! gpu_free_of_others "$GPU"; then          # a foreign job appeared -> AUTO-YIELD
      echo "$(hostname) gpu$GPU AUTO-YIELD (foreign job appeared) during $id"
      kill -TERM "$PYPID" 2>/dev/null; sleep 2; kill -KILL "$PYPID" 2>/dev/null
      yielded=1; break
    fi
  done
  wait "$PYPID" 2>/dev/null
  if [ "$yielded" -eq 1 ]; then
    rmdir "$CLAIMS/$id" 2>/dev/null; CUR_MOL=""    # release molecule so a clean GPU retries it
    exit 5                                         # trap releases the GPU lock
  fi
  if grep -q "DONE_MARKER" "$log" 2>/dev/null; then
    ran=$((ran+1)); fails=0
  else
    rmdir "$CLAIMS/$id" 2>/dev/null; fails=$((fails+1))
    [ "$fails" -ge 3 ] && { echo "$(hostname) gpu$GPU CIRCUIT-BREAKER (ran $ran)"; exit 4; }
  fi
  CUR_MOL=""
done < "$ITEMS"
echo "$(hostname) gpu$GPU finished, ran $ran"
