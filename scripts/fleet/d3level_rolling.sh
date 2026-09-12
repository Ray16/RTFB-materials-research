#!/bin/bash
# d3level_rolling.sh — refresh the D3TaLES-level comparison + both figures once an hour until the
# fleet has attempted all 468 molecules (then a final refresh and stop). Detached; logs to
# logs/fleet/d3level/rolling.log. Safety cap of 12 hours.
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO" || exit 1
CALC=results/d3tales_reorg_validation_d3level/calc
LOG=logs/fleet/d3level/rolling.log
echo "[$(date)] rolling updater started" >> "$LOG"
for i in $(seq 1 12); do
  sleep 3600
  bash scripts/fleet/d3level_refresh.sh >> "$LOG" 2>&1
  done=$(ls $CALC/*.json 2>/dev/null | wc -l)
  echo "[$(date '+%H:%M:%S')] hourly refresh #$i done ($done/468 attempted)" >> "$LOG"
  if [ "$done" -ge 468 ]; then
    echo "[$(date)] all 468 attempted -> final refresh complete, stopping rolling updater" >> "$LOG"
    break
  fi
done
