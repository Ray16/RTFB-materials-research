# Exact unsubstituted phenothiazine primary comparison

Fresh strict refinements now supply two independently checked stationary
geometries: [neutral minus](../phenothiazine_oxidized_geometry_2026-11-02/neutral_minus_geometry_review.json)
and [oxidized doublet](../phenothiazine_oxidized_geometry_2026-11-02/oxidized_geometry_review.json),
with all 53/43 actual SCF frames reviewed. The second neutral branch is running.
These are separate from the original-geometry electronic energies below.
New-geometry high-level scoring, Hessians/minimum and sulfur-column checks remain
pending; native G and oxidation potential stay blank.

Both exact phenothiazine electronic states now pass independent corrected-
integration, actual-wavefunction, spin, protocol and closed-resource review:
[neutral](neutral_electronic_review.json) and [oxidized](oxidized_electronic_review.json).
Corrected-minus-stored energies differ by at most 5.19e-11 Ha. These are fixed-
geometry electronic comparisons, not a redox potential or a native free energy.

The first neutral native calculation completed its Hessian and nine actual SCF
frames but failed the minimum gate. [Independent diagnosis](neutral_failed_minimum_diagnosis.json)
reconstructs a -25.47-cm-1 internal mode and the two standard sampled columns
(maximum relative error 7.75e-6). Both molecular and independent projections use
the same symmetric matrix for their comparison. Extra sulfur checks were not
reached, and no G is admitted. The [two signed seeds](neutral_repair_seeds.json)
retain the exact graph and await fresh optimization, corrected electronic
scoring and full native/S-column review. The original geometry and attempt are
unchanged.

The first oxidized native startup was interrupted during CUDA initialization
when another process became resident on its assigned GPU. The
[closure record](oxidized_startup_interruption.json) verifies every original
owner closed and no scientific program, SCF or Hessian was started. No foreign
process was signaled and no reservation was force-released. A fresh namespace
was prepared while that interrupted attempt remains preserved. That fresh attempt
completed its actual SCF but failed the unchanged stationarity gate:
0.00045735997 Ha/Bohr exceeds 0.00045. The [independent diagnosis and strict
refinement](../phenothiazine_oxidized_geometry_2026-11-02/README.md) preserve its
wavefunction and gradient; no Hessian or G was reached. Both neutral signed
branches and the oxidized state therefore require geometry refinement and new
own electronic/native calculations. The multi-input
coordinator now continues independent ready states while retaining a partner's
failure; it never automatically retries the failed calculation.

During review of the added sulfur probes, a unit-conversion mismatch was caught:
the pinned native engine uses 0.52917721067 Angstrom/Bohr, whereas the library's
constant is 0.52917721092. Added probes now use the actual recorded native center.
The first neutral failure preceded those probes. The unused oxidized preparation
and its original source bytes remain preserved; the selected recovery uses the
corrected coordinate convention. This is a validation-extension correction;
no nuclear or potential value was promoted from the earlier attempts.

[Current record](current.json) links the observed states and current refinement records;
[108 focused checks](all_checks.xml) cover the workflow, scope, source, ownership,
extra-column reconstruction and coordinate regression. Only validated 32-GB-
class devices run these stages; free T4 availability does not qualify it for
these molecular Hessians. Potential and G remain blank.

The exact product is N-(4-methylbenzyl)phenothiazine, a local site model for
post-functionalization of a preformed polymer. Its [prospective route
gate](route_gate.json) passes using the separately reviewed unsubstituted bead
precedent. This admission does not extend to the F/CN proposals, including the
newly reviewed 2-cyano precursor listing. The baseline's two equivalent registry
declarations use different SMILES spellings but the same exact structure and
neutral-singlet / +1-doublet hypotheses. A shared reader now accepts those
repeated sources while rejecting conflicting chemistry or charge/spin states.
All original declarations remain in the unified registry.

The archive contains optimized geometries and electronic records for these two
states, but no qualified own native thermochemistry. Their immutable source
records are preserved. Fresh installed/corrected electronic comparisons use the
existing worker, exact primary omegaB97M-V/def2-TZVPD/SMD MeCN protocol and
independent actual-wavefunction reviewer. Each state has 39 atoms and 976
high-level basis functions. The namespace, plan and source hashes are frozen.
No source-parent measurements or historical thermal terms are substituted.

The [calculation contract](calculation_contract.json) keeps the next gates
explicit: fresh stationarity and covalent/hydrogen topology; full solution
Hessian, rigid-mode projection and minima; both-step derivative checks; occupied
spin and frame checks; own nuclear/spin terms at 25/50/100-cm-1 rotor cutoffs;
and a context-matched corrected Fc/Fc+ reference. This sulfur-containing core
also requires S-centered x/y/z sampled columns, supplementing the standard two
columns. The intended 21 unique center/displaced native points are a prospective
requirement, not completed calculations.

Oxidation is defined as G(graft+)+G(Fc)-G(graft)-G(Fc+), one electron, in V versus
Fc/Fc+. The balanced two-versus-two solute difference cancels the common 1-M
translation conversion. Potential and native G remain blank. No second oxidation,
actual NaH/DMF grafting G, bead loading, coupling, transport or material admission
is inferred from this first stage.

Use `python -m redox.workflows.route_molecular_panel --help` for prepare,
one-pass dispatch and bounded queue entry points. The queue shares the same
owner-aware coordinator as saved-ionic validation, plus the existing atomic GPU
reservation, node dependency probe, CPU/RAM/storage admission and resource-only
retry accounting. Its advisory priority exposes one next pending input at a
time and cannot hold an extra slot for a started state. Actual started work is
preserved for independent review rather than automatically retried.

The [initial 61 focused checks](focused_checks.xml) cover exact state declarations,
held-route rejection, source namespaces, alias identity, shared queue states,
reservation-refusal accounting and one-active-input priority behavior. This is
preparation/execution infrastructure and scientific planning; it supplies no
new qualified redox potential or native free energy.
