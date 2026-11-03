# Neutral phenothiazine: strict signed-mode geometry refinement

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
[handoff records](controllers_actual_resume_verified.json) verify the resumed
controllers and the new geometry waiting process separately from scientific
starts. The original native oxidized-state priority remains ahead of these
refinements; no extra slot is reserved for them.

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
