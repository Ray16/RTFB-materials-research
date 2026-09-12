#!/bin/bash
# gpu_watchdog.sh [--dry-run] [--once] [--heal] [--run d3level] — cluster-wide contention enforcer.
#
# Every WATCHDOG_INTERVAL seconds it:
#   1) ENFORCE (auto-yield): on each ALLOWED host, kill OUR jobs on any GPU shared with another user;
#      on each RESERVED host, kill ALL our GPU jobs. Never touches other users' processes.
#   2) RELEASE orphans: drop GPU locks whose owner is dead, and molecule claims that are neither
#      complete nor currently running on any live node (so yielded work is retried elsewhere).
#   3) HEAL (only with --heal): relaunch workers on genuinely-free ALLOWED GPUs so yielded molecules
#      resume on clean hardware. A worker self-claims its GPU, so relaunching is idempotent.
#
# --dry-run reports what it WOULD kill/launch and changes nothing (use it to validate). --once runs a
# single cycle. Detach for continuous enforcement:  setsid nohup bash scripts/fleet/gpu_watchdog.sh
# --heal > logs/fleet/d3level/watchdog.log 2>&1 < /dev/null &
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO" || exit 1
# shellcheck disable=SC1091
source scripts/fleet/cluster.env 2>/dev/null || true
PY=/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python
SSH=(ssh -o BatchMode=yes -o ConnectTimeout=4 -o StrictHostKeyChecking=no)

DRY=0; ONCE=0; HEAL=0; RUN=d3level
while [ $# -gt 0 ]; do case "$1" in
  --dry-run) DRY=1;; --once) ONCE=1;; --heal) HEAL=1;; --run) RUN="$2"; shift;; esac; shift; done

FLEET="logs/fleet/$RUN"; CLAIMS="$FLEET/claims"; CALC="results/d3tales_reorg_validation_d3level/calc"
ITEMS="$FLEET/items.txt"; WORKER="scripts/fleet/${RUN}_fleet_worker.sh"
dryflag=""; [ "$DRY" -eq 1 ] && dryflag="--dry-run"

log(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

enforce(){
  local h out
  for h in $ALLOWED_HOSTS; do
    out=$("${SSH[@]}" "$h" "bash -s -- $dryflag" < scripts/fleet/gpu_enforce_host.sh 2>/dev/null | grep KILL)
    [ -n "$out" ] && log "$h: $out"
  done
  for h in $RESERVED_HOSTS; do
    out=$("${SSH[@]}" "$h" "bash -s -- $dryflag --reserved" < scripts/fleet/gpu_enforce_host.sh 2>/dev/null | grep KILL)
    [ -n "$out" ] && log "$h(reserved): $out"
  done
}

# ids currently running on ANY live node (so we never release a live molecule claim). Only ALLOWED
# hosts — we never have jobs on RESERVED (unusable) hosts, so skipping them avoids dead-host timeouts.
active_ids(){
  local h
  { pgrep -u "$(id -un)" -af 'validate_reorg_worker_d3tales.py' 2>/dev/null
    for h in $ALLOWED_HOSTS; do
      "${SSH[@]}" "$h" "pgrep -u $(id -un) -af validate_reorg_worker_d3tales.py" 2>/dev/null
    done
  } | grep -oE -- '--only [A-Z0-9]+' | awk '{print $2}' | sort -u
}

release_orphans(){
  # stale GPU locks (owner process dead on its host)
  local d o opid ohost id act; act=$(active_ids)
  for d in "$REPO"/logs/fleet/gpu_locks/*/; do
    [ -d "$d" ] || continue
    o=$(cat "$d/owner" 2>/dev/null); opid=$(echo "$o"|awk '{print $1}'); ohost=$(echo "$o"|awk '{print $3}')
    [ -z "$opid" ] && continue
    if ! "${SSH[@]}" "$ohost" "kill -0 $opid 2>/dev/null"; then
      [ "$DRY" -eq 1 ] && log "would release stale gpu-lock $(basename "$d")" || { rm -rf "$d"; log "released stale gpu-lock $(basename "$d")"; }
    fi
  done
  # orphaned molecule claims (not complete AND not running anywhere)
  local n=0
  for d in "$CLAIMS"/*/; do
    id=$(basename "$d"); [ "$id" = "*" ] && continue
    [ -f "$CALC/$id.json" ] && grep -q '"status":' "$CALC/$id.json" 2>/dev/null && continue
    echo "$act" | grep -qx "$id" && continue
    if [ "$DRY" -eq 1 ]; then n=$((n+1)); else rmdir "$d" 2>/dev/null && n=$((n+1)); fi
  done
  [ "$n" -gt 0 ] && log "$([ "$DRY" -eq 1 ] && echo would-release || echo released) $n orphaned molecule claims"
  return 0
}

remaining(){
  local total done; total=$(grep -c . "$ITEMS" 2>/dev/null || echo 0)
  done=$(grep -l '"status":' "$CALC"/*.json 2>/dev/null | wc -l)
  echo $(( total - done ))
}

heal(){
  [ "$HEAL" -eq 1 ] || return 0
  local rem; rem=$(remaining)
  [ "$rem" -le 0 ] && { log "heal: sweep complete, nothing to launch"; return 0; }
  local launched=0 host idx H
  while read -r host idx; do
    [ -z "$host" ] && continue
    if [ "$DRY" -eq 1 ]; then log "would launch $RUN worker on $host gpu$idx"; launched=$((launched+1)); continue; fi
    "${SSH[@]}" "$host" "cd $REPO && setsid nohup bash $WORKER $idx \
        > logs/fleet/$RUN/worker_\$(hostname -s)_gpu${idx}.log 2>&1 < /dev/null &" 2>/dev/null \
      && { log "launched $RUN worker on $host gpu$idx"; launched=$((launched+1)); }
    sleep 3
  done < <("$PY" scripts/free_gpus.py --all -n 8 2>/dev/null)
  [ "$launched" -eq 0 ] && log "heal: no free ALLOWED GPUs available right now"
}

cycle(){ enforce; release_orphans; heal; }

log "watchdog start (dry=$DRY once=$ONCE heal=$HEAL run=$RUN interval=${WATCHDOG_INTERVAL:-120}s)"
if [ "$ONCE" -eq 1 ]; then cycle; exit 0; fi
while true; do
  cycle
  # Auto-exit when there is nothing left to guard: sweep complete AND no workers of ours anywhere.
  if [ "$(remaining)" -le 0 ] && [ -z "$(active_ids)" ]; then
    log "sweep complete and no active workers -> watchdog exiting (nothing left to guard)"
    exit 0
  fi
  sleep "${WATCHDOG_INTERVAL:-120}"
done
