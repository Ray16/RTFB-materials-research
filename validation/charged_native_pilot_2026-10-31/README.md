# Charged-polymer native thermochemistry pilot

Nine exact viologen-polymer structures now pass independent stationary-geometry
checks: six single-site structures, both adjacent-pair configurations and one
one-spacer configuration. The [collection](summary.json) preserves the reviewed
SCF, topology, stereochemistry, geometry identity and resource results. One
remaining one-spacer optimization is separate. None of these geometry records
establishes Hessian minima, converged ensembles or interaction free energies.

The first charged native component pilot is running at the independently checked
114-atom site-2 geometry, charge +2 and singlet multiplicity. Actual process/start,
CPU affinity, six native thread caps, exclusive GPU residency, SCF iteration and
SCF charge/multiplicity are checked in [startup](startup.json). Its 1,100 AO and
213 occupied dimensions are below the checked neutral 1,193/235 envelope. The
aligned response workspace is three GiB; the prospectively bounded full job uses
80 GiB host RAM, a 31,744-MiB observed GPU ceiling, four GiB durable output, an
eight-hour wall bound and 36/32-GiB initial/live shared-storage margins. Only this
first charged pilot is admitted until its full physical/resource review passes;
aggregate shared-volume budgets also restrict further full native starts.

The maintained driver now declares charge/spin explicitly, retains neutral
control defaults, and requires a fresh full Hessian for this first ionic native
pilot. It reuses the corrected screening implementation, identical bounded
operators and original numerical settings. Sampled derivatives include the
protocol coordinates, largest-force coordinates and the actual graft N/benzyl-C
coordinates mapped from graph to XYZ order. SCF/gradient, full matrix, sampled
finite differences and modes remain subject to the unchanged gates.

Exact-geometry symmetry is obtained through the pinned xTB CLI from the projected
DFT Hessian; the CLI Hessian/mode filter is not used to supply nuclear free
energies. The full internal-mode qRRHO kernel provides all three rotor cutoffs.
Component G is assembled only from the fresh geometry-level electronic parent,
its own nuclear term and the explicit singlet spin term. The independent
reviewer checks actual SCF HDF geometry/state/energy, source identity, full
response rows, fresh gradients, projection, masses, symmetry and independent
nuclear equations before admitting component thermochemistry.

These are alternate smooth-surface SMD/acetonitrile component calculations. They
cannot be mixed with primary-model redox/grafting components, and a single
component is not a conditional occupancy cycle, ensemble or material decision.
The running pilot is unqualified until normal completion and independent review.
The [manifest](manifest.json) seals the source/plan/startup receipts and code.

Prepare a further eligible source with
`redox.workflows.charged_polymer_pilot --native-from-review REVIEW --out NEW_NS`,
then use the ordinary node probe and bounded gated launcher. Review normal
completion with `scripts/validation/review_charged_polymer_native.py`. Current
state and remaining requirements are recorded in the results index.
