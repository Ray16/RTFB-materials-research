#!/usr/bin/env bash
# Launch DFT+SMD for the experimental validation-core ids, one id per idle GPU.
# These ids are disjoint from the decorated-monomer workers already running on GPU6/7,
# so there is no result.json race. `ferrocene` is intentionally EXCLUDED here because the
# GPU7 worker already owns it (the internal reference is a single shared id).
#
#   GPUS=1,2,3,4,5 ./scripts/run_dft_valcore.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/.."
cd "$ROOT/src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox
export PYTHONPATH="$ROOT/src"
PY="$(command -v python)"
LOGDIR="$ROOT/calcs/dft"; mkdir -p "$LOGDIR"

# LD_LIBRARY_PATH: gpu4pyscf needs libcublas.so.12 etc. from torch's pip nvidia-* wheels.
SP="$($PY -c 'import site; print(site.getsitepackages()[0])')"
NVLIB=""
for d in cublas cusolver cusparse cuda_runtime cuda_nvrtc nccl cufft curand; do
  [ -d "$SP/nvidia/$d/lib" ] && NVLIB="$SP/nvidia/$d/lib:$NVLIB"
done
export LD_LIBRARY_PATH="$NVLIB${LD_LIBRARY_PATH:-}"

IDS=(methyl_viologen tempo_parent phenothiazine_parent anthraquinone_parent methylpyridinium)
GPUS="${GPUS:-1,2,3,4,5}"
IFS=',' read -ra GARR <<< "$GPUS"
[ "${#GARR[@]}" -lt "${#IDS[@]}" ] && { echo "!! need >= ${#IDS[@]} GPUs, got ${#GARR[@]}"; exit 1; }
echo ">> valcore DFT+SMD: ${#IDS[@]} ids on GPUs [$GPUS]"

pids=()
for i in "${!IDS[@]}"; do
  id="${IDS[$i]}"; g="${GARR[$i]}"
  echo ">>   $id -> GPU $g"
  CUDA_VISIBLE_DEVICES="$g" \
    "$PY" -m redox.dft --only "$id" --backend gpu \
      > "$LOGDIR/valcore_${id}.log" 2>&1 &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> valcore DFT workers done (fail=$fail)"; exit $fail
