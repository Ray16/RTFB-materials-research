#!/usr/bin/env bash
# Post-DFT finalization chain. Run ONLY after the DFT+SMD result.json files exist.
#
# This used to stop after descriptors+figures, which is exactly how a stale
# results/reorganization.csv could coexist with newer QC code in redox.properties.reorg (FINDINGS #20):
# the tables that feed ranking were never regenerated. It now runs the WHOLE chain, in
# dependency order, ending in the scorecard / Pareto / spreadsheet that people actually read.
#
#   ./scripts/pipeline/finalize_after_dft.sh   # consume existing calculations
# Compute missing cross points separately through a reservation-aware workflow.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$HERE/../.."
cd "$ROOT"

# Finalization owns derived products, not quantum calculations. Reject the old
# compute switch before changing any files or starting the environment.
if [ "${SKIP_REORG_CROSS:-1}" != "1" ]; then
  echo "Cross-point computation moved out of finalization. Use the reserved run_motifs.py workflow." >&2
  exit 2
fi

# $HOME is NOT the NFS repo home on every lambda node (on lambda5 it is /homes/rzhu, a
# different filesystem with no `redox` env). Source conda from the NFS path explicitly.
NFS=/nfs/lambda_stor_01/homes/rzhu
source "$NFS/miniforge3/etc/profile.d/conda.sh" && conda activate redox
export PYTHONPATH="$ROOT/src"

step() { echo; echo "== [$1] $2 =="; }

# FAIL-FAST: every step must succeed. A failed step used to print "!!" and continue, which
# is how stale-but-plausible tables reached the scorecard. Missing active-protocol energies
# are reported per couple as INCOMPLETE by the modules themselves (never silently zero).

step T "data-contract tests (thermo identities, caches, capacity paths, Pareto)"
python -m pytest -q tests

step 0 "all declared states: geometry, active-protocol energy, thermal coverage"
python -m redox.validation.readiness --strict

step 1 "level-matched Fc/Fc+ reference (same G definition + protocol as molecules)"
python scripts/pipeline/set_fc_reference.py

step 2 "redox E° table (active protocol; INCOMPLETE couples flagged)"
python -m redox.properties.potentials
python scripts/validation/orop/run_orop_benchmark.py --aggregate all   # OROP E° benchmark (feeds sigma prior)

step 3 "using existing reorganization cross points (no quantum calculations)"

step 4 "reorganization aggregation (lambda_O, lambda_R, 4-point sum + QC flag)"
python -m redox.properties.reorg --aggregate

step 5 "state integrity (bound + intact; necessary condition, not reversibility)"
python -m redox.properties.integrity

step 6 "stability: disproportionation (dG_disp) per intermediate"
python -m redox.properties.stability
python -m redox.validation.stability      # dG_disp vs experimental wave spacing -> sigma_disp
python scripts/validation/benchmark/score_benchmark_mecn.py   # E1/E2/spacing vs sourced MeCN set

step 7 "outer-sphere lambda_o consolidation (provenance-labelled)"
python scripts/pipeline/consolidate_lambda_outer.py

step 8 "capacity, SA score, solubility proxy"
python -m redox.properties.capacity_and_proxies


step 10 "scorecard (per-pool contiguous paths; INCOMPLETE propagation) -> Pareto"
python -m redox.screening.scorecard
python -m redox.screening.pareto

step 11 "candidate spreadsheet + result-dependent figures"
python scripts/pipeline/fill_candidates_xlsx.py
python scripts/plotting/candidates/plot_candidate_analysis.py
python scripts/plotting/validation/plot_validation_error.py
python scripts/plotting/validation/plot_orop_benchmark.py
python scripts/plotting/candidates/plot_candidates_2x2.py
python scripts/plotting/candidates/plot_merrifield_summary.py
python scripts/plotting/candidates/plot_pareto_front.py
python scripts/plotting/pipeline/plot_pipeline.py
python scripts/plotting/validation/plot_benchmark_mecn.py
python scripts/plotting/candidates/plot_discovery.py
python scripts/plotting/candidates/plot_leads_vs_existing.py

step 12 "independent audit: every published number re-derived from raw records"
python scripts/validation/audit/recompute_axes.py

echo; echo "== finalization complete =="
