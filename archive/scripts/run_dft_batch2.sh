#!/usr/bin/env bash
# DFT+SMD for the merrifield_multielectron batch, sharded ONE STATE PER GPU.
#
# Why per-state and not per-molecule: states of the same molecule are independent, so
# molecule-level sharding serializes 2-3 DFT optimizations behind each other and makes the
# wall-clock the SUM of them. Per-state, wall-clock is the slowest SINGLE state.
#
# Long-lived workers: each holds one GPU for its lifetime via `gpu_reserve acquire --pid $$`
# and releases on trap (the pattern CLAUDE.md prescribes), pulling tasks from an mkdir-based
# claim queue so the run is resumable and safe to re-launch.
#
# Dianions (red2, diffuse def2-tzvpd, slow SCF) are queued FIRST — longest-job-first keeps
# them off the critical path at the end.
set -uo pipefail
NFS=/nfs/lambda_stor_01/homes/rzhu
REPO="$NFS/0_redox"
GR="$NFS/bin/gpu_reserve"
cd "$REPO/src"
source "$NFS/miniforge3/etc/profile.d/conda.sh" && conda activate redox

RUN="${RUN:-batch2}"
THREADS="${THREADS:-2}"          # per worker; (threads x workers) must stay well under nproc
QDIR="$REPO/logs/fleet/$RUN"; CLAIMS="$QDIR/claims"; LOGS="$QDIR/logs"
mkdir -p "$CLAIMS" "$LOGS"
TASKS="$QDIR/tasks.txt"

# --- build the task list (pending states only; red2/dianions first) ------------------------
python - "$TASKS" <<'PY'
import csv, sys
from pathlib import Path
ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
IDS = ["bisviologen","aq_benzyloxy","nq_benzyloxy","dtbc_phenol","aq_benzylamino",
       "aq_benzyloxy_sa","nq_benzyloxy_sa","dtbc_phenol_sa","aq_benzylamino_sa"]
rows = list(csv.DictReader((ROOT/"library/manifest.csv").open()))
# A state already CLAIMED by a concurrently-running wave has no result.json yet, so the
# result-exists test alone would re-queue it and duplicate the work. Exclude every state
# claimed by any wave (claims live in logs/fleet/<run>/claims/<id>__<state>).
claimed = {d.name for d in ROOT.glob("logs/fleet/*/claims/*") if d.is_dir()}

# only states whose UMA pre-opt geometry exists (DFT seeds from it) and that are not done
pend = [r for r in rows if r["id"] in IDS
        and not (ROOT/f"calcs/dft/{r['id']}/{r['state']}/result.json").exists()
        and (ROOT/f"calcs/uma/{r['id']}/{r['state']}/relaxed.xyz").exists()
        and f"{r['id']}__{r['state']}" not in claimed]
if claimed:
    print(f"[queue] {len(claimed)} state(s) already claimed by a running wave - skipped")
n_wait = sum(1 for r in rows if r["id"] in IDS
             and not (ROOT/f"calcs/uma/{r['id']}/{r['state']}/relaxed.xyz").exists())
if n_wait:
    print(f"[queue] {n_wait} state(s) still awaiting UMA - re-run this script to pick them up")
# longest-job-first: most negative charge (diffuse basis, hardest SCF) goes first
pend.sort(key=lambda r: (int(r["charge"]), r["id"]))
Path(sys.argv[1]).write_text("".join(f"{r['id']}:{r['state']}\n" for r in pend))
print(f"[queue] {len(pend)} pending states -> {sys.argv[1]}")
for r in pend[:4]:
    print(f"        first: {r['id']}/{r['state']} q={r['charge']}")
PY

N_TASKS=$(wc -l < "$TASKS")
[ "$N_TASKS" -eq 0 ] && { echo ">> nothing pending"; exit 0; }

# --- reserve GPUs, one long-lived worker each ---------------------------------------------
NPROC=$(nproc)
mapfile -t FREE < <("$GR" list 2>/dev/null | awk -v h="$(hostname)" '$1==h {print $2}')
MAXW=$(( NPROC / (2*THREADS) ))                       # CPU-oversubscription guard
NW=${#FREE[@]}; [ "$NW" -gt "$MAXW" ] && NW=$MAXW
[ "$NW" -gt "$N_TASKS" ] && NW=$N_TASKS
[ "$NW" -eq 0 ] && { echo "!! no free GPU"; exit 1; }
echo ">> $N_TASKS tasks, $NW workers x $THREADS threads = $((NW*THREADS))/$NPROC cores"
echo ">> load before launch: $(cut -d' ' -f1-3 /proc/loadavg)"

worker() {
  local gpu="$1"
  "$GR" acquire "$gpu" --label "dft-$RUN" --pid $$ >/dev/null 2>&1 || return 0
  trap '"$GR" release "'"$gpu"'" >/dev/null 2>&1' EXIT INT TERM
  while read -r task; do
    [ -z "$task" ] && continue
    local id="${task%%:*}" st="${task##*:}"
    mkdir "$CLAIMS/${id}__${st}" 2>/dev/null || continue   # someone else has it
    echo "[gpu$gpu] $id/$st start $(date +%H:%M:%S)"
    OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS" OPENBLAS_NUM_THREADS="$THREADS" \
      CUDA_VISIBLE_DEVICES="$gpu" \
      python -m redox.dft --only "$id:$st" --backend gpu --nthreads "$THREADS" \
        > "$LOGS/${id}__${st}.log" 2>&1
    local rc=$?
    if [ -f "$REPO/calcs/dft/$id/$st/result.json" ]; then
      echo "[gpu$gpu] $id/$st DONE  $(date +%H:%M:%S)"
    else
      echo "[gpu$gpu] $id/$st FAIL rc=$rc $(date +%H:%M:%S)"
      rmdir "$CLAIMS/${id}__${st}" 2>/dev/null   # release so a retry can pick it up
    fi
  done < "$TASKS"
}

pids=()
for ((i=0; i<NW; i++)); do
  worker "${FREE[$i]}" & pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

DONE=$(python - <<'PY'
import csv
from pathlib import Path
ROOT=Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
IDS=["bisviologen","aq_benzyloxy","nq_benzyloxy","dtbc_phenol","aq_benzylamino",
     "aq_benzyloxy_sa","nq_benzyloxy_sa","dtbc_phenol_sa","aq_benzylamino_sa"]
rows=[r for r in csv.DictReader((ROOT/"library/manifest.csv").open()) if r["id"] in IDS]
print(sum((ROOT/f"calcs/dft/{r['id']}/{r['state']}/result.json").exists() for r in rows))
PY
)
echo ">> DFT $RUN finished: $DONE / 25 states complete"
