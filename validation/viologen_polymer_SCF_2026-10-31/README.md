# Exact viologen: charged polymer solver qualification

The exact route-bound baseline viologen now has a checked initial
**137-atom, +4-charge singlet polymer SCF and full gradient**. It is the
alternating-backbone, sites-2/4 two-graft scenario from the twelve audited
initial graphs. The open five-repeat backbone, unreacted benzyl chloride,
stereochemistry and benzyl--pyridinium N bonds retain their original audits.
No experimental graft loading or ion distribution is inferred.

The model contains the bare charged polymer with omitted separated spectator
ions in MeCN. It uses the separately declared experimental SWIG-CDS surface,
r2SCAN-D4/def2-SVP, the corrected grid-screening implementation, SCF thresholds
1e-12 Ha/1e-8 gradient and full integration-grid response. The PF6 dry-salt
composition reference remains separate; no explicit PF6 ions are introduced.

The first direct V100 attempt ran out of memory during SMD matrix construction,
before SCF iteration. Its [failed record](original_memory_failure.json) remains
preserved. The existing compact-buffer/contraction adapter keeps the same
surface, basis, precision and model equations. Two new namespaces compare it
against the direct calculation on a larger-memory card:

| Execution | GPU peak (MiB) | Initial maximum gradient (Ha/bohr) |
| --- | ---: | ---: |
| Compact V100 | 29,964 | 0.01396579005 |
| Direct H100 | 87,026 | 0.01396579001 |

The independent [compact](compact_review.json) and [direct](direct_review.json)
reviews check the actual SCF checkpoint energy, charge/spin, coordinates,
gradient-component closure, exact route gate, exclusive residency, source
bindings and host/GPU resources. The
[comparison](execution_agreement.json) gives a 3.64e-12-Ha energy difference
and 2.94e-10-Ha/bohr maximum force difference, within the 1e-7-Ha/1e-6-Ha/bohr
execution-agreement criteria. This is a resource/numerical control, not an
interaction energy or grafting free energy.

The initial force is about 31 times the 4.5e-4-Ha/bohr stationarity threshold.
Solvent optimization is therefore required. The versioned optimizer saves every
evaluated geometry, energy and gradient, a rolling SCF checkpoint, and the
optimizer-returned coordinates even when convergence is false. A fresh SCF and
force evaluation at the returned geometry prevents an earlier scanner energy
from being attached to a different structure. Topology, stereochemistry,
minimum/Hessian and conformer checks remain necessary before thermochemistry.

The first optimizer invocation rejected the GPU gradient object before any
geometry evaluation. The corrected invocation uses its gradient scanner, as in
the original geometry driver. That failed attempt remains preserved. The
current optimization namespace is:

`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-10-31_viologen_polymer_solution_optimization_v2/`

The compact/reference source namespaces are adjacent under
`2026-10-31_viologen_polymer_compact_pilot/` and
`2026-10-31_viologen_polymer_large_GPU_reference/`. Their plans freeze exact
input structures, prospective route evidence, model and resource settings.
Bulk checkpoints stay on shared storage; small checked records are versioned.

Use `python -m redox.workflows.charged_polymer_pilot` for preparation and the
shared guarded screening launcher for node probe/reservation/execution.
`--optimization-agreement` requires the checked paired resource result before
preparing an optimization. `review_charged_polymer_pilot.py` supplies separate
resource review and execution comparison.

No charged-polymer minimum, conformer ensemble, nuclear partition, neighboring
graft penalty, operating-state coupling or all-axis material pass follows from
this initial-state control. M1 and M4 remain open.
