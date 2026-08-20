#!/usr/bin/env bash
# Post-DFT finalization chain. Run ONLY after all DFT+SMD result.json exist.
# Fast/CPU-light bookkeeping: reference -> E° table -> descriptors -> figures.
#
#   ./scripts/finalize_after_dft.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/.."
cd "$ROOT"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox
export PYTHONPATH="$ROOT/src"

echo "== [1/4] level-matched Fc/Fc+ reference from ferrocene DFT+SMD =="
python scripts/set_fc_reference.py

echo "== [2/4] redox E° table (DFT+SMD columns + validation-core rows) =="
python -m redox.redox

echo "== [3/4] structure descriptors (RMSD from DFT geometries) =="
python -m redox.descriptors

echo "== [4/4] discussion figures =="
python scripts/plot_results.py

echo "== finalization complete =="
