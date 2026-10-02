#!/usr/bin/env bash
# Fan the composite DFT+SMD run (r2SCAN-D4 opt // wB97X-D3 energy, in-solvent) across GPUs.
# One worker per free GPU on THIS host, reserved through the system-wide gate (gpu_reserve run,
# which sets CUDA_VISIBLE_DEVICES); modulo-sharded over manifest
# states. Resumable: states with result.json are skipped. Falls back to CPU with BACKEND=cpu.
#
#   ./scripts/pipeline/run_dft.sh                 # auto-detect idle GPUs, fan out
#   GPUS=1,2,3,4,5,6,7 ./scripts/pipeline/run_dft.sh
#   BACKEND=cpu THREADS=8 ./scripts/pipeline/run_dft.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/../.."
cd "$ROOT/src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox
export PYTHONPATH="$ROOT/src"
PY="$(command -v python)"
LOGDIR="$ROOT/calcs/dft"; mkdir -p "$LOGDIR"

BACKEND="${BACKEND:-gpu}"

if [ "$BACKEND" = "cpu" ]; then
  THREADS="${THREADS:-8}"
  NCORE="$($PY -c 'import os; print(os.cpu_count())')"
  NW=$(( NCORE / THREADS )); [ "$NW" -lt 1 ] && NW=1
  echo ">> DFT+SMD (CPU): $NW workers x $THREADS threads"
  pids=()
  for i in $(seq 0 $((NW-1))); do
    OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
      "$PY" -m redox.qm.dft --all --shard "$NW:$i" --backend cpu --nthreads "$THREADS" \
        > "$LOGDIR/shard_$i.log" 2>&1 &
    pids+=($!)
  done
  fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
  echo ">> CPU shards done (fail=$fail)"; exit $fail
fi

# --- GPU path ---
# LD_LIBRARY_PATH fix: gpu4pyscf-cuda12x needs libcublas.so.12 etc., which ship inside
# torch's pip nvidia-* packages but aren't on the loader path by default.
SP="$($PY -c 'import site; print(site.getsitepackages()[0])')"
NVLIB=""
for d in cublas cusolver cusparse cuda_runtime cuda_nvrtc nccl cufft curand; do
  [ -d "$SP/nvidia/$d/lib" ] && NVLIB="$SP/nvidia/$d/lib:$NVLIB"
done
export LD_LIBRARY_PATH="$NVLIB${LD_LIBRARY_PATH:-}"

GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve   # system-wide reservation gate (CLAUDE.md)
# Free GPUs on THIS host per the gate (zero resident processes, not reserved) unless GPUS is set.
if [ -z "${GPUS:-}" ]; then
  GPUS="$("$GR" list 2>/dev/null | awk -v h="$(hostname -s)" '$1==h{printf "%s%s",(c++?",":""),$2}')"
fi
[ -z "$GPUS" ] && { echo "!! no idle GPUs found (set GPUS=...)"; exit 1; }
IFS=',' read -ra GARR <<< "$GPUS"
NW="${#GARR[@]}"
echo ">> DFT+SMD (GPU) composite r2SCAN-D4//wB97X-D3 : $NW workers on GPUs [$GPUS]"

# CPU thread cap per GPU worker: gpu4pyscf runs the SCF on the GPU but does CPU-side integral
# screening/DF that defaults to ALL cores -> N workers would oversubscribe the shared node
# (CLAUDE.md: no CPU contention). Cap and guard (threads x workers) << nproc.
GPU_THREADS="${GPU_THREADS:-4}"
NPROC="$($PY -c 'import os;print(os.cpu_count())')"
if [ $(( GPU_THREADS * NW )) -gt $(( NPROC / 2 )) ]; then
  GPU_THREADS=$(( NPROC / (2 * NW) )); [ "$GPU_THREADS" -lt 1 ] && GPU_THREADS=1
  echo ">> capped GPU_THREADS to $GPU_THREADS to avoid CPU oversubscription (nproc=$NPROC)"
fi
echo ">> $GPU_THREADS CPU threads/worker = $((GPU_THREADS*NW))/$NPROC cores; load $(cut -d' ' -f1-3 /proc/loadavg)"
pids=()
for i in "${!GARR[@]}"; do
  g="${GARR[$i]}"
  "$GR" run "$g" -- env OMP_NUM_THREADS="$GPU_THREADS" MKL_NUM_THREADS="$GPU_THREADS" \
    OPENBLAS_NUM_THREADS="$GPU_THREADS" PYTHONPATH="$ROOT/src" LD_LIBRARY_PATH="$LD_LIBRARY_PATH" \
    "$PY" -m redox.qm.dft --all --shard "$NW:$i" --backend gpu \
      > "$LOGDIR/gpu_${g}.log" 2>&1 &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> GPU workers done (fail=$fail)"; exit $fail
