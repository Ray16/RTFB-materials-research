#!/bin/bash
# gpu_probe.sh [self_user] — per-GPU ownership-aware status for THIS host, one CSV line per GPU:
#     idx,util,mem,foreign_owners
# where foreign_owners is a space-joined, de-duplicated list of OTHER users (not self_user) with a
# compute process resident on that GPU (empty = nobody else). self_user defaults to the current
# user, so over ssh it correctly evaluates as the remote user (same account on the shared NFS env).
#
# This is the single source of "is GPU X free of other users" used by free_gpus.py, gpu_guard.sh,
# and gpu_watchdog.sh. A GPU with a foreign owner is NEVER free, regardless of util/mem (an idle-
# but-resident foreign job — e.g. 314 MiB at 0% — must not be treated as available).
set -u
self="${1:-$(id -un)}"

declare -A INFO      # bus -> "idx,util,mem"
declare -A FOREIGN   # bus -> " owner owner ..."  (pre-initialized so `set -u` is version-safe)
while IFS=',' read -r i b u m; do
  b=$(echo "$b" | tr -d ' ')
  INFO[$b]="$(echo "$i" | tr -d ' '),$(echo "$u" | tr -d ' '),$(echo "$m" | tr -d ' ')"
  FOREIGN[$b]=""
done < <(nvidia-smi --query-gpu=index,gpu_bus_id,utilization.gpu,memory.used \
                    --format=csv,noheader,nounits 2>/dev/null)

while IFS=',' read -r b p _m; do
  b=$(echo "$b" | tr -d ' '); p=$(echo "$p" | tr -d ' ')
  [ -z "$p" ] && continue
  [ -n "${INFO[$b]+x}" ] || continue
  o=$(ps -o user= -p "$p" 2>/dev/null | tr -d ' '); [ -z "$o" ] && o="unknown"
  [ "$o" != "$self" ] && FOREIGN[$b]="${FOREIGN[$b]} $o"
done < <(nvidia-smi --query-compute-apps=gpu_bus_id,pid,used_memory \
                    --format=csv,noheader 2>/dev/null)

for b in "${!INFO[@]}"; do
  IFS=',' read -r idx util mem <<< "${INFO[$b]}"
  f=$(echo "${FOREIGN[$b]}" | tr ' ' '\n' | grep -v '^$' | sort -u | tr '\n' ' ' | sed 's/ *$//')
  echo "$idx,$util,$mem,$f"
done | sort -t, -k1,1n
