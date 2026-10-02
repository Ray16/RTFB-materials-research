#!/usr/bin/env bash
# Uniform-basis (def2-TZVPD) energy campaign worker: one shard per GPU.
#   setsid bash scripts/fleet/sp_fleet_worker.sh <N> <I> <gpu_idx>
#   CROSS_IDS="a b c" ... -> cross points only for those groups (default: every group)
# Runs (1) active-protocol single points for pipeline states, (2) reorg cross points,
# (3) the OROP benchmark states — all resumable (records are protocol-addressed files).
# The GPU is taken through the system-wide gate (gpu_reserve); CPU threads are capped.
set -uo pipefail
N=$1; I=$2; GPU=$3
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO"
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh && conda activate redox
export PYTHONPATH="$REPO/src" OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve run "$GPU" -- bash -c "
  python -m redox.qm.sp --all --shard $N:$I --backend gpu
  if [ -n \"${CROSS_IDS:-}\" ]; then
    k=0; for g in ${CROSS_IDS:-}; do
      [ \$((k % $N)) -eq $I ] && python -m redox.properties.reorg --only \$g --backend gpu
      k=\$((k+1)); done
  else
    python -m redox.properties.reorg --all-cross --shard $N:$I --backend gpu
  fi
  python -m redox.qm.sp --orop --shard $N:$I --backend gpu
"
echo "[worker $N:$I gpu$GPU] done $(date)"
