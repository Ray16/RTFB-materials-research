#!/usr/bin/env bash
# UMA pre-opt for the merrifield_multielectron batch: one worker per candidate id, each on a
# GPU reserved through the system-wide gate (gpu_reserve), each CPU-thread-capped.
#
# NOTE: on lambda5 $HOME is /homes/rzhu, which is a DIFFERENT filesystem from the NFS repo
# home — the `redox` env lives under the NFS path, so conda is sourced from there explicitly
# rather than from ~ (this is why scripts/pipeline/run_uma.sh does not work here).
set -uo pipefail
NFS=/nfs/lambda_stor_01/homes/rzhu
REPO="$NFS/0_redox"
GR="$NFS/bin/gpu_reserve"
cd "$REPO/src"
source "$NFS/miniforge3/etc/profile.d/conda.sh" && conda activate redox

IDS="${IDS:-bisviologen aq_benzyloxy nq_benzyloxy dtbc_phenol aq_benzylamino}"
THREADS="${THREADS:-2}"
LOGDIR="$REPO/calcs/uma"; mkdir -p "$LOGDIR"

echo ">> load before launch: $(cut -d' ' -f1-3 /proc/loadavg) on $(nproc) cores"
mapfile -t SLOTS < <("$GR" list 2>/dev/null | awk -v h="$(hostname)" '$1==h {print $2}')
echo ">> $(echo "$IDS" | wc -w) ids, ${#SLOTS[@]} free GPU slots on $(hostname)"
[ "${#SLOTS[@]}" -eq 0 ] && { echo "!! no free GPU"; exit 1; }

pids=(); i=0
for id in $IDS; do
  g="${SLOTS[$((i % ${#SLOTS[@]}))]}"; i=$((i+1))
  OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
    "$GR" run "$g" -- python -m redox.uma --only "$id" --device cuda \
      > "$LOGDIR/batch2_${id}.log" 2>&1 &
  pids+=($!); echo "   $id -> gpu $g (log: calcs/uma/batch2_${id}.log)"
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> UMA batch2 done (fail=$fail)"
exit $fail
