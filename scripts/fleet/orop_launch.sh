#!/bin/bash
# orop_launch.sh <gpu>... — start one detached OROP-sweep worker per GPU on THIS node.
# Each worker reserves its GPU through gpu_reserve (no-ops if busy), so this is safe to run
# on several nodes. Keep the count CONSERVATIVE (a few GPUs/node) to leave headroom.
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO" || exit 1
H=$(hostname); mkdir -p logs/fleet/orop
for g in "$@"; do
  setsid nohup bash scripts/fleet/orop_worker.sh "$g" \
      > "logs/fleet/orop/worker_${H}_gpu${g}.log" 2>&1 < /dev/null &
  sleep 4
done
echo "launched $# orop workers on $H (gpus: $*)"
