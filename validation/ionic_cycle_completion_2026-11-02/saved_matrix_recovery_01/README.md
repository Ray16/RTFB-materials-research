# Fourth ionic cycle: saved-matrix validation recovery

The alternating adjacent double-graft calculation stopped during a displaced
SCF checkpoint write with HDF5 file-lock error 11. The original worker, wrapper
and budget observer are confirmed closed. Its failure, full Hessian, full solvent
H1, nine successfully solved response blocks (411 rows), exact orbital/grid frame
and nineteen completed geometry evaluations remain preserved in their namespace.
No G is admitted from that interrupted attempt.

The [source evidence](source_evidence.json) checks the exact ownership/resource
chain, absence of foreign residency or budget stop, response residuals, solvent
H1 frame hashes, full finite symmetric FP64 matrix and stationary center.
A separate CPU projection reconstructs six rigid modes and a lowest internal
frequency of 9.23351 cm⁻¹, with no unstable internal mode below the gate. This
projection alone does not qualify the matrix or a thermodynamic component.

The [new namespace](namespace.json) and [plan](plan.json) retain the original
method, graph, charge +4, singlet state, MeCN model, physical thresholds and full
H100/host/storage resource class. A fresh center must reproduce the source
energy, gradient and exact geometry before the saved matrix can be loaded.
The recovery must perform all 25 center/displaced SCF evaluations in the new namespace;
all sampled columns, actual wavefunctions, modes, symmetry and own nuclear
terms still require independent review. The saved response is labeled explicitly
as reused, never newly solved. No HDF5 locking setting is disabled.

Preparation is `python -m redox.workflows.charged_hessian_recovery --help`.
The bounded controller is `python scripts/fleet/run_saved_ionic_validation.py
--help`; it reuses the existing node probe, reservation, CPU/RAM allocation and
aggregate storage dispatcher. It records live ownership, preserves uncertain
attempts, stops at scientific failure, bounds prelaunch refusals and only resumes
a controller after exact closure. The append-only shared storage roster includes
this recovery. The raw controller status is linked in [current.json](current.json).
A submission receipt is distinct from an actual SCF start and independent QC.

The first preflight caught a neutral-only 33-GiB storage rule being applied to
the charged 38-GiB envelope, before any scientific worker started. The launcher
now chooses the explicit scope and still validates the complete charged resource
class. The refused placement/log and stopped controller are retained. The
corrected controller was explicitly resumed after verified closure. The same
source change stopped new electronic-queue admissions as designed; that exact
controller was also confirmed closed and resumed, preserving running workers.

The [focused checks](focused_checks.xml) pass 91 tests, including matrix/state
mismatch rejection, incomplete or failed response rows, uncertain/live ownership,
foreign residency, wrong failure type, resource bounds and the distinct neutral
versus ionic storage rules. Actual frozen-source preparation additionally passes
source hashes and original response-frame checks. The completed slide/manuscript
build bindings remain intact. The fourth G cycle, occupied ensembles, electronic
coupling, salt speciation, material admission and JCIM submission remain pending.
