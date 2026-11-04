# Exact viologen: finite stationary electronic occupancy cycles

Eleven independently checked stationary components provide three complete
atom/charge-balanced electronic cycles for the exact viologen design. All terms
use the same r2SCAN-D4/def2-SVP, smooth-surface SMD/acetonitrile model and identical
corrected-screening source/kernel identities. Both matched empty references were
freshly optimized and independently checked. Old seed-selection energies do not
enter these sums. See [inputs and cycle receipts](cycles.json).

| Controlled configuration | Graft pattern | Electronic nonadditivity (kJ/mol) |
| --- | --- | ---: |
| All-same | Adjacent sites 2/3 | -3.461 |
| Alternating | Adjacent sites 2/3 | +1.202 |
| Alternating | One spacer, sites 2/4 | -13.482 |
| All-same | One spacer, sites 2/4 | Pending one stationary component |

Each value is E(both)+E(empty)-E(first)-E(second), from one stationary conformer
per graph. These are whole-chain electronic nonadditivities. They include
backbone/cap and deformation changes; they are not isolated motif binding
energies or electronic couplings. The signs and differences are not converged
ensemble or experimental predictions. Unknown graft loading and tacticity
probabilities are not inferred.

All delta-delta-G fields remain blank. No nuclear/spin correction is borrowed
from another structure or the preceding implementation, and no E value is
substituted for missing G. The ongoing native component pilot and later
conformer sampling remain required for thermodynamic interpretation and
material selection. None of these results establishes feasibility along all axes.

The analyzer checks each actual independent source binding, exact graph/state,
model, ion convention, implementation source/kernel and integer atom/charge
balance. Missing components remain explicitly pending. The [44 focused tests](tests.xml)
include bad stoichiometry, nonfinite inputs, absent components, changed kernel
bytes and rejection of E-to-G promotion. The [manifest](manifest.json) freezes
this scope and its as-used code; the raw calculation inputs remain at the bound
namespace locations.
