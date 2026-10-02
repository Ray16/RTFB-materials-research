#!/bin/bash
# dft_launch.sh <N> <i0> <gpu>... — multi-node DFT+SMD fan-out with ONE global shard numbering.
# Starts shards i0, i0+1, ... of `python -m redox.qm.dft --all --shard N:i --backend gpu` on THIS
# node, one detached worker per GPU, each reserved through the system-wide gate (gpu_reserve run)
# with a CPU thread cap (THREADS, default 2). run_batch is resumable (states with result.json are
# skipped), so give each node a disjoint i-range of the same N:
#   ssh lambda1 "bash .../dft_launch.sh 15 0 0 1 2 3 4 5 6"      # shards 0-6
#   ssh lambda5 "bash .../dft_launch.sh 15 7 0 1 2 3 4 5 6 7"    # shards 7-14
# Logs: logs/fleet/dft/<host>_gpu<g>.log
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve
THREADS="${THREADS:-2}"
cd "$REPO" || exit 1
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh && conda activate redox
PY="$(command -v python)"
# gpu4pyscf-cuda12x needs libcublas.so.12 etc. from torch's pip nvidia-* packages (as run_dft.sh)
SP="$($PY -c 'import site; print(site.getsitepackages()[0])')"
NVLIB=""
for d in cublas cusolver cusparse cuda_runtime cuda_nvrtc nccl cufft curand; do
  [ -d "$SP/nvidia/$d/lib" ] && NVLIB="$SP/nvidia/$d/lib:$NVLIB"
done
N=$1; i=$2; shift 2; H=$(hostname -s); mkdir -p logs/fleet/dft
for g in "$@"; do
  setsid nohup "$GR" run "$g" -- env OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
    OPENBLAS_NUM_THREADS="$THREADS" PYTHONPATH="$REPO/src" LD_LIBRARY_PATH="$NVLIB${LD_LIBRARY_PATH:-}" \
    "$PY" -m redox.qm.dft --all --shard "$N:$i" --backend gpu \
    > "logs/fleet/dft/${H}_gpu${g}.log" 2>&1 < /dev/null &
  i=$((i+1)); sleep 2
done
