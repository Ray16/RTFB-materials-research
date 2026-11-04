# Neutral phenothiazine: strict signed-mode geometry refinement

All three strict geometry refinements now pass independent review, including the
[neutral plus branch](../phenothiazine_oxidized_geometry_2026-11-02/neutral_plus_geometry_review.json).
The [fresh electronic campaign](../refined_electronic_2026-11-02/README.md) retains
both neutral branches and the oxidized state; one actual V100 start is verified.
Own electronic/native qualification remains pending. Earlier observations below
retain their original waiting/running counts.

The [minus branch](../phenothiazine_oxidized_geometry_2026-11-02/neutral_minus_geometry_review.json)
now passes independent review of all 53 saved SCF frames and the strict final
gradient (3.50e-7 Ha/Bohr). Its partner has a separately verified actual start;
the [current record](current.json) retains both. High-level electronic scoring
and own native/Hessian/S-column checks at the new geometry remain pending, and
neither old-geometry energy nor nuclear terms transfer.

The [current record](current.json) links both exact neutral-singlet signed
branches and their live bounded controller. The source first-native calculation
failed its minimum gate; [independent reconstruction](../phenothiazine_primary_2026-11-02/neutral_failed_minimum_diagnosis.json)
finds a -25.47-cm-1 internal mode, verifies nine actual SCF frames and repeats
the two standard sampled columns. No native G or potential is admitted from that
saddle. Its source geometry, full Cartesian Hessian and failed attempt remain
unchanged.

The two seeds follow opposite signs of the independently reconstructed mode,
with maximum atom displacement 0.1 Angstrom. Both preserve the exact
N-(4-methylbenzyl)phenothiazine graph and neutral singlet assignment. They share
the existing strict geometry calculation and review code with chloride controls,
while their own admission profile retains the correct charge, multiplicity,
39-atom/397-AO dimension and validated 32-GB device class. The three existing T4
chloride calculations keep their original source-frozen workers and continue.

The calculation uses corrected r2SCAN-D4/def2-SVP in primary SMD MeCN. The
optimizer requires energy 1e-9 Ha, gmax 1e-5 and grms 5e-6 Ha/Bohr, with drms
1e-4 and dmax 2e-4. Every geometry/SCF evaluation and optimizer return is saved
on durable shared storage. Converged optimizer output still requires a fresh
final SCF/gradient, exact graph check, spin check and independent reconstruction
of actual wavefunctions and strict final stationarity.

The existing quantum, resource, source-guard and ownership checks are retained.
A small explicit adapter lets the shared remote dispatcher call the maintained
geometry launcher instead of adding another scheduler. Other workflow rows keep
their original launcher. Three admission coordinators were gracefully stopped
for that shared-code update; their exact owners were confirmed closed before
explicit resume. No scientific worker was signaled. The
[handoff records](controllers_actual_resume_verified.json) preserve that earlier
observation. The current v2 paths and owners replace the earlier neutral queue:
archive preflight found hidden checkout evidence before either scientific start.
The [recovery account](../phenothiazine_oxidized_geometry_2026-11-02/README.md)
records closure of those controllers, unchanged evidence bytes and successful
private-namespace probes. Both neutral branches are retained. The separately
diagnosed oxidized-doublet gradient failure now has its own strict refinement
queue and one advisory pending-input priority; it cannot reserve an extra slot.

The [85 focused checks](focused_checks.xml) cover exact state/method/resource
admission, rejection of T4 substitution or arbitrary launcher modules, legacy
launcher compatibility, queue ownership and the existing chloride contracts.
Prepared, submitted, optimized and minimum-qualified remain separate states.
After a new geometry passes independent review, its corrected electronic energy
and complete own native Hessian/S-column/partition checks must be computed at
that geometry. The archived-geometry energy and any prior nuclear terms do not
transfer. Both branches are retained until those checks support comparison.

Run `python -m redox.workflows.route_geometry --help` to inspect preparation,
launch and queue commands. Independent geometry review remains
`python scripts/validation/review_counterion_geometry.py --help`, whose output
explicitly identifies the molecular route-geometry scope rather than an ion-pair
association calculation. Native G and oxidation potential stay blank here.
