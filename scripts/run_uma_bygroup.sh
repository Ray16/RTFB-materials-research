#!/usr/bin/env bash
# Fan the UMA conformer x spin-multiplicity SCAN across GPUs, one process per redox group
# (one model load each; each process relaxes all conformers x mults for its group's states).
# Robust single-process-per-GPU parallelism. Resumable (skips states with result.json).
#
#   MODEL=uma-s-1p2p1 ./scripts/run_uma_bygroup.sh
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/.."
source ~/miniforge3/etc/profile.d/conda.sh
RUN="conda run -n redox --no-capture-output"

MODEL="${MODEL:-uma-s-1p2p1}"
NGPU="$($RUN python -c 'import torch; print(max(1, torch.cuda.device_count()))')"
$RUN python -c "import sys; sys.path.insert(0,'$ROOT/src'); from redox.uma import ensure_registered; ensure_registered('$MODEL')"

# groups straight from the filesystem (no CSV parsing)
GROUPS=()
for d in "$ROOT"/library/*/conformers; do GROUPS+=("$(basename "$(dirname "$d")")"); done
echo ">> groups: ${GROUPS[*]}  | model=$MODEL | $NGPU GPUs"

i=0; pids=()
for g in "${GROUPS[@]}"; do
  gpu=$(( i % NGPU ))
  ( cd "$ROOT/src" && PYTHONUNBUFFERED=1 CUDA_VISIBLE_DEVICES="$gpu" \
    $RUN python -m redox.uma --only "$g" --model "$MODEL" --device cuda ) \
      > "$ROOT/calcs/uma/bygroup_$g.log" 2>&1 &
  pids+=($!); i=$((i+1))
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> all groups done (fail=$fail)"