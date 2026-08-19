#!/usr/bin/env bash
# Fan the UMA relaxation of all manifest states across every available GPU.
# Each GPU runs one process (one model load) handling a shard of the states.
# Resumable: states with an existing result.json are skipped.
#
#   MODEL=uma-s-1p1 ./scripts/run_uma.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE/../src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox

MODEL="${MODEL:-uma-s-1p1}"
NGPU="$(python -c 'import torch; print(max(1, torch.cuda.device_count()))')"
echo ">> fanning UMA ($MODEL) across $NGPU GPU(s)"

pids=()
for i in $(seq 0 $((NGPU-1))); do
  CUDA_VISIBLE_DEVICES="$i" python -m redox.uma --model "$MODEL" --device cuda \
      --shard "$NGPU:$i" > "$HERE/../calcs/uma/shard_$i.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> all shards done (fail=$fail)"
exit $fail
