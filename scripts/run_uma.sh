#!/usr/bin/env bash
# Fan UMA relaxation across IDLE GPUs only, one worker per GPU, each CPU-thread-capped so we
# NEVER oversubscribe the shared node (see CLAUDE.md: no GPU *or* CPU contention — it impacts
# other users). MLIP work runs on the GPU, so a small CPU thread cap is correct.
# Resumable: states with an existing result.json are skipped.
#
#   ./scripts/run_uma.sh                       # auto-pick idle GPUs, 2 threads/worker
#   MODEL=uma-s-1p2p1 THREADS=2 ./scripts/run_uma.sh
#   GPUS=0,1,2 ./scripts/run_uma.sh            # force a GPU set (still thread-capped)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE/../src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox

MODEL="${MODEL:-uma-s-1p2p1}"
THREADS="${THREADS:-2}"      # CPU threads PER worker; keep SMALL — the compute is on the GPU

# Pick IDLE GPUs only (free_gpus.py skips GPUs another job is using) unless GPUS is set.
if [ -z "${GPUS:-}" ]; then
  GPUS="$(python "$HERE/free_gpus.py" 2>/dev/null || true)"
fi
[ -z "$GPUS" ] && { echo "!! no idle GPU free — wait, or set GPUS=... explicitly"; exit 1; }
IFS=',' read -ra GARR <<< "$GPUS"
NW="${#GARR[@]}"

# Hard guard against CPU oversubscription: (threads x workers) must stay well under nproc.
NPROC="$(nproc)"
if [ $(( THREADS * NW )) -gt $(( NPROC / 2 )) ]; then
  echo "!! THREADS($THREADS) x workers($NW) = $((THREADS*NW)) would oversubscribe nproc=$NPROC."
  echo "   Lower THREADS or narrow GPUS."; exit 1
fi
echo ">> UMA ($MODEL): $NW workers on GPUs [$GPUS], $THREADS threads each = $((THREADS*NW))/$NPROC cores"
echo "   load before launch: $(cut -d' ' -f1-3 /proc/loadavg)"

pids=()
for i in "${!GARR[@]}"; do
  g="${GARR[$i]}"
  CUDA_VISIBLE_DEVICES="$g" OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
    OPENBLAS_NUM_THREADS="$THREADS" \
    python -m redox.uma --model "$MODEL" --device cuda --shard "$NW:$i" \
      > "$HERE/../calcs/uma/shard_${g}.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> all shards done (fail=$fail)"
exit $fail
