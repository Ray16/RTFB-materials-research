#!/usr/bin/env bash
# Unattended finisher: wait for all DFT+SMD systems, then finalize + plot + push.
#
# Polls until every manifest state has a result.json AND no redox.dft workers remain,
# then runs the finalization chain, the validation-error bar plot, and commits/pushes
# the small result artifacts. Designed to survive user logout (launch with nohup).
#
#   nohup ./scripts/wait_and_finalize.sh > calcs/dft/wait_and_finalize.log 2>&1 &
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/.."
cd "$ROOT"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox
export PYTHONPATH="$ROOT/src"

POLL="${POLL:-60}"          # seconds between checks
MAX_HOURS="${MAX_HOURS:-24}"
DEADLINE=$(( $(date +%s) + MAX_HOURS*3600 ))

# expected state count = data rows in the manifest (minus header)
EXPECTED=$(( $(wc -l < library/manifest.csv) - 1 ))
echo "$(date '+%F %T') >> waiting for $EXPECTED DFT+SMD systems (poll ${POLL}s, max ${MAX_HOURS}h)"

while :; do
  DONE=$(find calcs/dft -name result.json | wc -l)
  WK=$(pgrep -fc "redox.dft" || true)
  echo "$(date '+%F %T')    done=$DONE/$EXPECTED  workers=$WK"
  # Complete when every system has a result AND no worker is still churning.
  if [ "$DONE" -ge "$EXPECTED" ] && [ "${WK:-0}" -eq 0 ]; then
    echo "$(date '+%F %T') >> all systems complete"
    break
  fi
  if [ "$(date +%s)" -ge "$DEADLINE" ]; then
    echo "$(date '+%F %T') !! deadline hit at done=$DONE/$EXPECTED workers=$WK — finalizing what exists"
    break
  fi
  sleep "$POLL"
done

echo "$(date '+%F %T') == finalize chain =="
bash scripts/finalize_after_dft.sh || echo "!! finalize_after_dft.sh returned non-zero"

echo "$(date '+%F %T') == validation-error bar plot =="
python scripts/plot_validation_error.py || echo "!! plot_validation_error.py returned non-zero"

echo "$(date '+%F %T') == commit + push small artifacts =="
# Only small, reproducible artifacts + code — never calc outputs, weights, or tokens.
git add -A scripts config src results/*.csv results/figures/*.png 2>/dev/null || true
if ! git diff --cached --quiet; then
  git commit -m "DFT+SMD batch complete: E° table, descriptors, validation error bars

Auto-finalized by wait_and_finalize.sh after all $EXPECTED systems converged.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>" \
    && git push && echo "$(date '+%F %T') >> pushed"
else
  echo "$(date '+%F %T') >> nothing new to commit"
fi

echo "$(date '+%F %T') == DONE =="
