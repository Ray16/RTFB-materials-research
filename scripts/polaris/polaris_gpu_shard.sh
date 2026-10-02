#!/bin/bash
# Per-rank launcher for the Polaris reorg run: one MPI rank per A100 (4/node). Maps the PALS
# local rank -> CUDA_VISIBLE_DEVICES, and uses the global rank/size as the worker's shard
# (--offset global_rank --stride total_ranks). The worker is resumable, so re-submitting after
# a walltime cutoff just continues where it left off.
LR=${PMI_LOCAL_RANK:-0}
GR=${PMI_RANK:-0}
NT=${PMI_SIZE:-1}
export CUDA_VISIBLE_DEVICES=$LR
export OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
PY=/grand/FRAME-IDP/rzhu/env/redox/bin/python
cd /grand/FRAME-IDP/rzhu/redox_reorg
echo "rank $GR/$NT on $(hostname) GPU $LR"
PYTHONPATH=/grand/FRAME-IDP/rzhu/redox_reorg/src exec $PY scripts/validation/reorg_d3tales/validate_reorg_worker_prod.py \
    --source data/molecules.csv --offset "$GR" --stride "$NT" --backend gpu
