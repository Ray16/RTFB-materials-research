#!/bin/bash
# gpu_enforce_host.sh [--dry-run] [--reserved] — runs ON a host (locally or via `ssh host bash -s`).
# AUTO-YIELD: kill OUR GPU compute jobs that violate the shared-node policy —
#   * --reserved : this host is off-limits (e.g. lambda13) -> yield ALL our GPU jobs here
#   * otherwise  : yield only our jobs on a GPU that ANOTHER user is also using (collision)
# Never touches other users' processes. Prints one line per victim:
#   KILL <pid> gpu<idx> <reason>        (DRYRUN KILL ... when --dry-run: reports only, kills nothing)
set -u
DRY=0; RES=0
for a in "$@"; do case "$a" in --dry-run) DRY=1;; --reserved) RES=1;; esac; done
self=$(id -un)

declare -A IDX
while IFS=',' read -r i b; do IDX[$(echo "$b"|tr -d ' ')]=$(echo "$i"|tr -d ' '); done \
  < <(nvidia-smi --query-gpu=index,gpu_bus_id --format=csv,noheader 2>/dev/null)

declare -A OWN
while IFS=',' read -r b p _m; do
  b=$(echo "$b"|tr -d ' '); p=$(echo "$p"|tr -d ' '); [ -z "$p" ] && continue
  [ -n "${IDX[$b]+x}" ] || continue
  u=$(ps -o user= -p "$p" 2>/dev/null | tr -d ' '); [ -z "$u" ] && u="unknown"
  OWN[$b]="${OWN[$b]:-} ${u}:${p}"
done < <(nvidia-smi --query-compute-apps=gpu_bus_id,pid,used_memory --format=csv,noheader 2>/dev/null)

for b in "${!OWN[@]}"; do
  list="${OWN[$b]}"; idx="${IDX[$b]:-?}"
  hasself=$(echo "$list" | tr ' ' '\n' | grep -c "^${self}:")
  hasother=$(echo "$list" | tr ' ' '\n' | grep -v "^${self}:" | grep -c ':')
  reason=""
  if [ "$RES" -eq 1 ] && [ "$hasself" -ge 1 ]; then reason="reserved-host"
  elif [ "$hasself" -ge 1 ] && [ "$hasother" -ge 1 ]; then reason="collision"; fi
  [ -z "$reason" ] && continue
  for up in $list; do
    u="${up%%:*}"; pid="${up##*:}"
    [ "$u" = "$self" ] || continue
    if [ "$DRY" -eq 1 ]; then
      echo "DRYRUN KILL $pid gpu$idx $reason"
    else
      echo "KILL $pid gpu$idx $reason"
      kill -TERM "$pid" 2>/dev/null
      ( sleep 3; kill -KILL "$pid" 2>/dev/null ) &
    fi
  done
done
