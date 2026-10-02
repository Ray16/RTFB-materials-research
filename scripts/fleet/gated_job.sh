#!/bin/bash
# gated_job.sh <gpu> <log> <python args...> — run ONE python job on THIS node through the
# system-wide reservation gate (gpu_reserve run), CPU-thread-capped (THREADS, default 2), with
# the CUDA library path gpu4pyscf needs (as run_dft.sh). Detached: launch with `ssh -n -f`.
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve
THREADS="${THREADS:-2}"
cd "$REPO" || exit 1
source /nfs/lambda_stor_01/homes/rzhu/miniforge3/etc/profile.d/conda.sh && conda activate redox
PY="$(command -v python)"
SP="$($PY -c 'import site; print(site.getsitepackages()[0])')"
NVLIB=""
for d in cublas cusolver cusparse cuda_runtime cuda_nvrtc nccl cufft curand; do
  [ -d "$SP/nvidia/$d/lib" ] && NVLIB="$SP/nvidia/$d/lib:$NVLIB"
done
g=$1; log=$2; shift 2; mkdir -p "$(dirname "$log")"
setsid nohup "$GR" run "$g" -- env OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
  OPENBLAS_NUM_THREADS="$THREADS" PYTHONPATH="$REPO/src" LD_LIBRARY_PATH="$NVLIB${LD_LIBRARY_PATH:-}" \
  "$PY" "$@" > "$log" 2>&1 < /dev/null &
