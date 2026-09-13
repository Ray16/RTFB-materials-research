#!/bin/bash
# d3level_refresh.sh — one-shot: re-aggregate the D3TaLES-level sweep + regenerate BOTH figures
# from whatever results exist so far (safe to run repeatedly while the fleet fills in).
set -u
REPO=/nfs/lambda_stor_01/homes/rzhu/0_redox
cd "$REPO" || exit 1
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox 2>/dev/null
export PYTHONPATH=src
CALC=results/d3tales_reorg_validation_d3level/calc
n=$(ls $CALC/*.json 2>/dev/null | wc -l)
echo "[$(date '+%H:%M:%S')] refresh: $n/468 d3level results"
python scripts/aggregate_d3tales_reorg_prod.py 2>&1 | grep -E "aggregated|MATCH|production|shift"
python scripts/finalize_reorg_qc.py 2>&1 | grep -E "excluded|reliable-only"   # fold in unbound-anion QC
python scripts/plotting/plot_d3tales_reorg_prod_parity.py 2>&1 | grep -v findfont | tail -1
python scripts/plotting/plot_d3level_hole_electron.py 2>&1 | grep -v findfont | tail -1
