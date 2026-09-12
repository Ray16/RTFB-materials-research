#!/bin/bash
# d3level_launch.sh <gpu>... — start one detached D3TaLES-level worker per GPU on THIS node.
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO" || exit 1
H=$(hostname); mkdir -p logs/fleet/d3level
for g in "$@"; do
  setsid nohup bash scripts/fleet/d3level_fleet_worker.sh "$g" \
      > "logs/fleet/d3level/worker_${H}_gpu${g}.log" 2>&1 < /dev/null &
  sleep 4
done
echo "launched $# d3level workers on $H (gpus: $*)"
