#!/usr/bin/env bash
# Launch the D3TaLES reorg validation batch across idle GPUs on THIS node.
# One worker per idle GPU, stride-sharded so load balances and reruns skip finished molecules.
# Shared-node rules: only uses GPUs that are actually idle, caps CPU threads per worker.
#
#   bash scripts/launch_reorg_validation.sh [MAX_WORKERS]
set -u
cd "$(dirname "$0")/.."
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox

MAXW="${1:-6}"   # leave some GPUs free for other users by default
mkdir -p tmp/reorg_val

# idle GPUs: <5% util AND <500 MiB used
mapfile -t IDLE < <(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader \
                    | awk -F, '$2+0<5 && $3+0<500 {print $1}')
N=${#IDLE[@]}
(( N > MAXW )) && N=$MAXW
if (( N == 0 )); then echo "no idle GPUs; aborting"; exit 1; fi
echo "idle GPUs: ${IDLE[*]}  -> using $N workers (stride=$N)"

for (( k=0; k<N; k++ )); do
  gpu=${IDLE[$k]}
  CUDA_VISIBLE_DEVICES=$gpu OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
    PYTHONPATH=src nohup python scripts/validate_reorg_worker.py --offset $k --stride $N \
    > tmp/reorg_val/worker_$k.log 2>&1 &
  echo "  worker $k on GPU $gpu (pid $!)"
done
echo "launched $N workers. progress: ls results/d3tales_reorg_validation/calc/ | wc -l"
