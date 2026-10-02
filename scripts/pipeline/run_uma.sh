#!/usr/bin/env bash
# Fan UMA relaxation across free GPUs on THIS host, one worker per GPU, each reserved through
# the system-wide gate (gpu_reserve run) and CPU-thread-capped so we
# NEVER oversubscribe the shared node (see CLAUDE.md: no GPU *or* CPU contention — it impacts
# other users). MLIP work runs on the GPU, so a small CPU thread cap is correct.
# Resumable: states with an existing result.json are skipped.
#
#   ./scripts/pipeline/run_uma.sh                       # auto-pick idle GPUs, 2 threads/worker
#   MODEL=uma-s-1p2p1 THREADS=2 ./scripts/pipeline/run_uma.sh
#   GPUS=0,1,2 ./scripts/pipeline/run_uma.sh   # restrict to these GPUs (still gated + thread-capped)
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE/../../src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox

MODEL="${MODEL:-uma-s-1p2p1}"
THREADS="${THREADS:-2}"      # CPU threads PER worker; keep SMALL — the compute is on the GPU

GR=/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve   # system-wide reservation gate (CLAUDE.md)
# Free GPUs on THIS host per the gate (zero resident processes, not reserved) unless GPUS is set.
if [ -z "${GPUS:-}" ]; then
  GPUS="$("$GR" list 2>/dev/null | awk -v h="$(hostname -s)" '$1==h{printf "%s%s",(c++?",":""),$2}')"
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
  "$GR" run "$g" -- env OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
    OPENBLAS_NUM_THREADS="$THREADS" PYTHONPATH="$HERE/../../src" \
    python -m redox.qm.uma --model "$MODEL" --device cuda --shard "$NW:$i" \
      > "$HERE/../../calcs/uma/shard_${g}.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> all shards done (fail=$fail)"
exit $fail
