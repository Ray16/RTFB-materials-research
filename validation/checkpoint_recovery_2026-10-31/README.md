# Checkpoint recovery and six stationary charged-polymer geometries

Six single-site viologen polymer geometries now pass independent fresh final-SCF,
stationarity, full indexed connectivity, stereochemistry, resource and ownership
checks. They cover sites 2, 3 and 4 in both controlled backbone configurations.
The open five-repeat graph retains ten backbone carbons, nine backbone C--C
bonds and two ends. The [summary](summary.json) links every reviewed geometry to
its exact candidate alias and frozen optimizer return. Each is one conformer;
Hessian minimum, ensemble and intersite interaction/free-energy checks remain
incomplete. None is a material admission.

The six optimizations converged before a GPU allocation failure in the separate
final solve. The retained optimizer/scanner held its mean field while another
was allocated. New source code releases those references first. For the closed
original runs, the new finalization entry point evaluates the full-precision
returned coordinates in a separate process. It does not repeat or replace the
converged optimization. All source evaluations, original failure records and
returned coordinates remain unchanged. Fresh final forces and actual SCF
checkpoint coordinates/energies must agree, and actual connectivity/stereochemistry
must be retained before the stationary-geometry record can be collected.

The geometry-level model is r2SCAN-D4/def2-SVP with the explicitly separate
smooth-surface SMD model in acetonitrile. Coordinates contain bare charged
polymers with omitted separated spectator ions. Dry-salt mass bookkeeping and
primary-model reaction thermochemistry are separate. These geometry energies
are not new redox potentials, grafting G or polymer occupancy G.

The three full-Hessian controls saved finite symmetric matrices before their
resource stops: one V100 reached its original four-hour wall bound and the
restored V100/H100 runs failed during the next derivative calculation. A retained
bound gradient object kept the restored mean field alive; prospective source now
releases it. The new continuation accepts only a proven closed resource failure
with a complete source matrix and all original response rows passing. It binds
the exact state/model and matrix hash, performs a fresh center solve, and runs
the original sampled derivatives, projection and minimum/thermodynamic checks.
The new phase has a prospective one-GiB output bound, two-hour wall allowance,
80-GiB host bound and 33/32-GiB initial/live storage margins. It still uses the
shared gate and disjoint CPU/thread allocation.

Two implementation faults in the recovery layer are preserved in
[source_failures.json](source_failures.json): a solvent-wrapper shallow copy
inherited an instance override and recursed; after fixing that, an inherited
adapter incorrectly demanded a new CPHF solve in a saved-matrix phase. The copy
behavior now has a regression test. The saved-matrix phase uses the same bounded
operators and explicitly records no new response solve, while the independent
reviewer checks the complete source response records. Full recomputation retains
its original new-response assertion. No physical derivative, stationarity, spin
or minimum threshold is relaxed.

The two completed derivative panels from the adapter-assertion attempts each
contain eight checks; their maximum relative errors are below 9e-6. They remain
failed execution records. The corrected phase-specific controls are separate
namespaces, listed in [current_validation_manifest.json](current_validation_manifest.json).
They must complete normally and pass independent reconstruction before a final
reproduction comparison or native G claim is made.

The [test report](tests.xml), source snapshots and [manifest](manifest.json)
freeze this scope. Raw quantum records and all original attempts remain at the
locations in the summary/manifests. Restore a closed converged optimizer return
with `redox.workflows.charged_polymer_pilot --finalize-from`; prepare post-Hessian
validation with `redox.workflows.screening_control prepare-Hessian-validation`.
The ordinary node probe and bounded gated launcher still perform submission.
