# Candidate, grafting and oligomer evidence audit

This scoped audit replays cached arithmetic and applies the current prospective
route gate. It does not certify all computed properties or launch new quantum
calculations. [Manifest](manifest.json) records source paths and SHA-256 hashes.
The archived material remains at `/nfs/lambda_stor_01/homes/rzhu/archive/0_redox`;
its original absolute paths may need explicit relocation when used by code.

- The [default 75-candidate feed](default_gate.json) has 75 holds and zero
  prospective proceeds under the current evidence registry. This feed does not
  include the separately assessed TEMPO exact-support case; its [archived review](../../../archive/0_redox/results/grafting/tempo_exact_support_2026-10-26/README.md)
  passes the route gate, while its overall material decision is pending.
- The [AQ carboxylate panel](aq_gate.json) and [PMDI linker panel](pmdi_gate.json)
  each contain four held proposals. These are previously prepared exact graphs,
  not newly calculated candidates. Missing evidence is listed per candidate.
- The [portable cached replay](cached_replay.json) reproduces the capped PMDI
  first reduction as −1.075790864 V versus computed Fc/Fc+, HCl-reference
  grafting ΔG as −0.291571 kJ/mol, and DBU separated-ion grafting ΔG as
  −103.068826 kJ/mol. This checks saved-state arithmetic, identity and balance;
  it does not rerun quantum calculations or independently validate physical QC.
  The base-assisted value remains conditional on the specified ion/solvent model.
- The [archived exact PMDI cycle](../../../archive/0_redox/results/reviews/2026-10-25/pmdi_neighbor_cycle_checked/README.md)
  reports −1.8616 kJ/mol conditional nonadditivity at 50 cm⁻¹ for four
  single-minimum states. The [local phthalimide interaction](../../../archive/0_redox/results/reviews/2026-10-26/local_motif_interactions/README.md)
  reports −24.76 kJ/mol in MeCN for a frozen capped-fragment triplet. These are
  different chemical models and quantities; neither establishes converged
  populations or electronic coupling. Counterpoise and broader sampling remain
  pending. The separate-model [60-state replay](../publication_native_replay_2026-10-28/README.md)
  cannot be mixed into the primary PMDI cycle.

Reproduce with `PYTHONPATH=src python scripts/pipeline/screen_grafting_routes.py`
using the manifest's feed, evidence and assessment (the latter is the script's
project-root default), with a fresh `--out`; supply each panel using
`--candidates`. Run the portable arithmetic with
`python3 -I -S examples/cached_thermochemistry/reproduce.py`. The copied default
feed contains historical absolute references to the old project root; use the
archive's original records and explicit source relocation instead of treating
those paths as current raw data. No inference of actual running jobs is made
from historical campaign status files.
