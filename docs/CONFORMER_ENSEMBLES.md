# Finite conformer ensembles for conditional polymer comparisons

An ensemble free energy requires native free energies for the included minima,
not just electronic energies or the number of times a search found a structure.
For one chemical graph, backbone configuration, solvent, temperature and protocol,
the implemented finite partition is

\[
G_{\mathrm{set}}=G_{\min}-k_BT\ln\sum_i
\exp[-(G_i-G_{\min})/(k_BT)].
\]

Each unique minimum contributes once. The per-state free energies already contain
the electronic spin partition and the detected rotational symmetry number; these
factors are not applied again. The normalized terms describe conditional model
populations within the included set. They are not measured polymer populations,
metadynamics trajectory weights or weights for an experimentally unknown sequence.

The [partition implementation](../src/redox/properties/conformer_ensemble.py)
uses a shifted log-sum to handle large absolute quantum energies. New arithmetic
uses exact SI constants for eV-to-kJ/mol conversion. Missing native G leaves the
ensemble blank. The calling workflow must independently check electronic and
thermal parents, exact graph, stereochemistry and protocol consistency.

## Duplicate minima and sampling sensitivity

Coordinate comparisons first map the declared graph to the actual XYZ atom order.
Symmetry-aware heavy-atom RMSD uses proper rotations. A complete-link group requires
every pair to have RMSD at most 0.10 Å and electronic energies within 0.20 kJ/mol.
Complete-link grouping prevents a chain of small pairwise differences from merging
distinct endpoint structures. One minimum-G member represents each duplicate group;
rediscovery does not create an artificial entropy contribution.

These are numerical duplicate heuristics. The collector repeats the grouping at
0.05, 0.10 and 0.20 Å and reports the resulting cycle sensitivity. The treatment
does not resolve a full hindered-rotor partition or guarantee separation of all
conformational basins. Adding finite wells to a qRRHO treatment remains an
approximate conformer model, particularly for soft backbone motions.

For the wider searches, the predeclared diagnostic compares the first 10 and 20
DFT starts: minimum electronic energy, minimum G and ensemble G must each change
by at most 2 kJ/mol. Missing records cannot satisfy this condition. This is a
within-search comparison, not proof of global convergence. The two backbone
configurations remain separate; no experimental tacticity weights are assigned.
The search software's own conformer filtering is a separate operation; see the
[CREST documentation](https://crest-lab.github.io/crest-docs/page/documentation/keywords.html).

## Two-pose baseline and broader searches

The original phthalimide five-repeat panel contains six occupancy graphs with
two checked DFT geometries each. Its six previously selected native G values are
preserved. The [all-pose native campaign](../results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2/README.md)
calculates the other six G values. It cannot establish adequate sampling by itself.
The broader 12 CREST searches cover all six graphs at two backbone configurations
and remain an independent source of additional DFT starts.

The new collector reconstructs each available native G from its active electronic
record and projected nuclear/spin partition. It requires the same geometry as the
energy parent, current native admission and one scientific partition protocol.
It verifies graph/stereochemistry before using a geometry in the RMSD matrix.
Previous checked values are copied as source-bound receipts, not recomputed or
silently replaced.

For a single backbone configuration the neighbor cycle is

\[
\Delta\Delta G=G_{ij}+G_0-G_i-G_j.
\]

Here each term may be a complete declared finite-set partition. The quantity is
conditional site nonadditivity. It includes conformational and scaffold effects
and does not isolate motif binding or determine an intersite electron-transfer
rate. Reagent terms cancel from this balanced cycle; an absolute grafting free
energy still requires the defined reaction and qualified reagent terms.

## Execution and limits

The all-pose worker uses the existing qualified bounded PCM implementation,
finer CDS Hessian derivative and independent worst-column checks. A new native
pilot gates expansion. Every GPU job uses the shared reservation gate, two disjoint
CPU cores and matching native thread limits. The conservative host allowance is
256 GiB per concurrent Hessian worker: a bounded 128-GiB curvature cache and
128-GiB headroom. Available host RAM therefore constrains concurrency in addition
to suitable free GPUs. This is a resource bound, not a molecular accuracy claim.

`scripts/pipeline/polymer_conformer_thermals.py` prepares, collects or runs the
finite baseline. `coordinate_conformer_thermals.py` gives a bounded continuation
after the pilot was submitted. It preserves failures and never relaxes the
physical gates or retries an invalid state automatically.

## Distinct smooth-surface continuation

The separately identified smooth-surface MeCN model is tested as an independent
protocol; it is not substituted into the primary-model reaction cycles. The
current refinement declaration covers twenty separated starts for each of twelve
graph/configuration sets (240 total). Geometry, hybrid electronic energy and
native thermochemistry are separately checked stages. Each stage retains the
exact source geometry and model identity. Occupied graphs require their own
hybrid and native resource pilots before expansion.

`scripts/pipeline/smooth_polymer_native.py` writes full baseline and displaced
SCF/gradient records, the native Hessian and the nuclear partition at three rotor
cutoffs. Its finite-difference panel includes the protocol's two columns plus
the two largest-force coordinates, with duplicates removed, at both declared
steps. It uses the geometry's xTB-detected rotational symmetry number. A separate
CPU reviewer reconstructs the column errors, internal modes and nuclear sums.
These checks establish a state within the named model, not a complete ensemble
or predictive solvation accuracy.

`prepare_smooth_native_batch.py` transfers only a checked execution class to
other exact-parent states of the same constitution, basis, charge/spin, model,
settings and runtime. `coordinate_smooth_native_batch.py` waits for the independent
native pilot before declaring the cohort. It then waits for each electronic
review before preparing a new native state. An electronic RAM-recovery option
must refer to the same geometry; no energy or thermal term is borrowed. Only an
identical completed reference is reused. Actual failures and uncertain ownership
remain preserved; resource contention does not create a scientific retry.

The current native execution bound is 80 GiB host RSS plus 32 GiB free-RAM
headroom and a GPU-peak limit of 31,744 MiB on a minimum 32-GB-class card.
This differs from the historical baseline adapter's allocation policy above.
Until its large-state pilot passes, the new bound remains prospective and cannot
certify other conformers. CPU-only dependency probes defer GPU-library imports
to the reserved worker. Available GPUs, disjoint CPU slots and remaining host
memory promises jointly determine concurrency.
