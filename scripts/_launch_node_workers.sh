#!/usr/bin/env bash
# Launch N local reorg-validation workers on this node with global stride sharding.
# args: BASE_RANK  N_LOCAL  W_TOTAL  GPU_CSV
# Each worker: unique global rank (BASE_RANK+i), stride W_TOTAL, pinned to one GPU, 4 CPU threads.
set -u
cd /nfs/lambda_stor_01/homes/rzhu/0_redox
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox
BASE=$1; NLOC=$2; W=$3; IFS=',' read -ra GPUS <<< "$4"
mkdir -p tmp/reorg_val
host=$(hostname)
for (( i=0; i<NLOC; i++ )); do
  rank=$((BASE + i)); gpu=${GPUS[$i]}
  CUDA_VISIBLE_DEVICES=$gpu OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
    PYTHONPATH=src nohup python scripts/validate_reorg_worker.py --offset $rank --stride $W \
    > tmp/reorg_val/w_${rank}_${host}.log 2>&1 &
done
# CPU sanity: report thread budget vs cores
echo "$host: launched $NLOC workers (ranks $BASE..$((BASE+NLOC-1))) on GPUs [$4]; " \
     "thread budget = $((NLOC*4)) / $(nproc) cores; loadavg=$(cut -d' ' -f1 /proc/loadavg)"
