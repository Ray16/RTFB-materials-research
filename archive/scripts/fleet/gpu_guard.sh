#!/bin/bash
# gpu_guard.sh — sourced helpers that ENFORCE the shared-node GPU policy (see cluster.env + CLAUDE.md
# "Submitting GPU jobs"). Source this from a fleet worker AFTER `cd`-ing to the repo root:
#     source scripts/fleet/gpu_guard.sh
# Provides:
#   gpu_free_of_others <idx>  -> 0 if NO other user is resident on physical GPU idx, else 1
#   gpu_claim <idx>           -> atomic cluster-wide reservation (one task per GPU). 0 got it / 1 taken
#   gpu_release <idx>         -> release a reservation we own
_GUARD_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_GUARD_REPO="$(cd "$_GUARD_DIR/../.." && pwd)"
# shellcheck disable=SC1090
source "$_GUARD_DIR/cluster.env" 2>/dev/null || true
GPU_LOCKS="$_GUARD_REPO/logs/fleet/gpu_locks"
mkdir -p "$GPU_LOCKS" 2>/dev/null

_guard_host(){ hostname -s 2>/dev/null || hostname | cut -d. -f1; }

# 0 if physical GPU <idx> has NO process owned by another user; 1 if a foreign job is resident.
gpu_free_of_others(){
  local idx="$1" line foreign
  line=$(bash "$_GUARD_DIR/gpu_probe.sh" 2>/dev/null | awk -F, -v i="$idx" '$1==i{print; exit}')
  [ -z "$line" ] && return 0                 # GPU not visible -> nothing to collide with
  foreign=$(echo "$line" | cut -d, -f4)      # 4th field = foreign owners (empty = clean)
  [ -z "$foreign" ]
}

# Atomic cluster-wide, one-task-per-GPU reservation via NFS mkdir. Reclaims a lock whose owner
# process on THIS host is dead (crash cleanup). Returns 0 if we hold it, 1 if another live task does.
gpu_claim(){
  local idx="$1" d o opid ohost
  d="$GPU_LOCKS/$(_guard_host)_gpu${idx}"
  if mkdir "$d" 2>/dev/null; then
    printf '%s %s %s\n' "$$" "$(date +%s)" "$(_guard_host)" > "$d/owner" 2>/dev/null
    return 0
  fi
  o=$(cat "$d/owner" 2>/dev/null)
  opid=$(echo "$o" | awk '{print $1}'); ohost=$(echo "$o" | awk '{print $3}')
  if [ "$ohost" = "$(_guard_host)" ] && [ -n "$opid" ] && ! kill -0 "$opid" 2>/dev/null; then
    printf '%s %s %s\n' "$$" "$(date +%s)" "$(_guard_host)" > "$d/owner" 2>/dev/null   # reclaim stale
    return 0
  fi
  return 1
}

# Release a reservation, but only if WE own it (never delete another task's lock).
gpu_release(){
  local idx="$1" d opid
  d="$GPU_LOCKS/$(_guard_host)_gpu${idx}"
  opid=$(awk '{print $1}' "$d/owner" 2>/dev/null)
  if [ -z "$opid" ] || [ "$opid" = "$$" ] || ! kill -0 "$opid" 2>/dev/null; then
    rm -rf "$d" 2>/dev/null
  fi
}
