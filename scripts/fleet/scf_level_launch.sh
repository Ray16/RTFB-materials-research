#!/bin/bash
# [SET=reference|candidates] [EXTRA="--reverse --claim"] scf_level_launch.sh <N> <i0> <gpu>... — start shards i0, i0+1, ... of the SCF/level check on THIS
# node, one detached worker per GPU, each reserved through the system-wide gate (gpu_reserve run)
# with a 2-thread CPU cap. Logs: logs/validation/scf_level/<host>_gpu<g>.log
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve
cd "$REPO" || exit 1
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh && conda activate redox
N=$1; i=$2; shift 2; H=$(hostname -s); mkdir -p logs/validation/scf_level
for g in "$@"; do
  setsid nohup "$GR" run "$g" -- env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
    PYTHONPATH=src python scripts/validation/reference/scf_level_check.py --set "${SET:-reference}" --shard "$N:$i" ${EXTRA:-} \
    > "logs/validation/scf_level/${SET:-reference}${EXTRA:+_claim}_${H}_gpu${g}.log" 2>&1 < /dev/null &
  [ -n "${EXTRA:-}" ] || i=$((i+1))   # claim mode: every worker walks the full list
  sleep 2
done
