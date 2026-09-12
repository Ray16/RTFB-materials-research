#!/bin/bash
# reorg_launch.sh <gpu> [<gpu> ...] — start one detached reorg fleet worker per GPU on THIS node.
# Staggered to avoid an NFS import stampede. Invoke per node via ssh:
#   ssh <host> "cd /nfs/lambda_stor_01/homes/rzhu/0_redox && bash scripts/fleet/reorg_launch.sh 0 1 2"
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
WORKER=scripts/fleet/reorg_fleet_worker.sh
cd "$REPO" || exit 1
H=$(hostname)
mkdir -p logs/fleet/reorg
for g in "$@"; do
  setsid nohup bash "$WORKER" "$g" > "logs/fleet/reorg/worker_${H}_gpu${g}.log" 2>&1 &
  sleep 4
done
echo "launched $# reorg workers on $H (gpus: $*)"
