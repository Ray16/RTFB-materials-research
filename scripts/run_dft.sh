#!/usr/bin/env bash
# Fan the DFT+SMD single points across CPU cores (PySCF is CPU-bound).
# Partitions cores into workers; each worker runs a shard of states with a fixed thread
# slice (no oversubscription). Resumable: states with result.json are skipped.
#
#   XC=b3lyp BASIS=def2-svp THREADS=4 ./scripts/run_dft.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE/../src"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox

XC="${XC:-b3lyp}"; BASIS="${BASIS:-def2-svp}"; THREADS="${THREADS:-4}"
NCORE="$(python -c 'import os; print(os.cpu_count())')"
NW=$(( NCORE / THREADS )); [ "$NW" -lt 1 ] && NW=1
echo ">> DFT+SMD $XC/$BASIS : $NW workers x $THREADS threads on $NCORE cores"

pids=()
for i in $(seq 0 $((NW-1))); do
  OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" \
    python -m redox.dft --all --shard "$NW:$i" --xc "$XC" --basis "$BASIS" \
      --nthreads "$THREADS" > "$HERE/../calcs/dft/shard_$i.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo ">> all DFT shards done (fail=$fail)"
exit $fail
