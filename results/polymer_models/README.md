# Polymer-fragment and interaction models

All result categories: [results index](../README.md).

**Open the [rotatable 3D interaction report](index.html)** for the current pilot.
[Interaction progress](interaction_status.json) and [reviewed results](interaction_checked.json)
are separate from the [559-graph inventory](current.csv) and its [status](status.json).

The prepared inventory contains explicitly hydrogen-capped polystyrene fragments
with one, two or three repeat units. Backbone stereochemical alternatives are
enumerated; existing pendant stereochemistry is retained. These are covalent
oligomers, not two whole polymer beads. The one-repeat reference contains the same
backbone construction as the higher oligomers and differs from the earlier
methyl-capped graft-site proxy.

For representative quinone, AQ, NDI, PMDI, Cata-NDI and viologen motifs, matched
single-site and unfunctionalized-backbone controls are also prepared. They support
an atom/charge-balanced site-cooperativity cycle on the same two-site backbone:

`deltaDeltaG = G(P_AA) + G(P_ClCl) - G(P_ACl) - G(P_ClA)`.

This measures nonadditivity of grafting two sites, including conformational and
local interaction effects. It is not a pure pair binding energy. Directly
subtracting twice a hydrogen-capped monomer energy from a covalent dimer is not
a balanced interaction reaction.

Two separate polymer-segment complexes will be treated as a distinct interchain
model, with separation/orientation sampling and an association standard state.
Entire-bead contacts, swelling and crosslinked-network behavior need a larger-scale
model. The current requested screen targets local pendant-site stacking and
contacts between separate polymer segments; it does not simulate whole beads.

Redox-site charge distributions, spin states, counterions, chain conformation and
stereochemistry must be tracked. First and second molecular redox waves are not
interchangeable with the first and second electrons added to a multi-site cluster.
Statistical site degeneracy must be separated from genuine electronic/cooperative
splitting before interpreting oligomer potential shifts. Charged viologen oligomers
require an explicit counterion/ionic-environment treatment before physical ranking.

## Active stacking-method pilot

The [frozen finite roster](campaigns/2026-10-19_pi_stacking/plan.json) contains
27 neutral electronic-structure cases for NDI, PMDI and phthalimide controls:

- One-repeat monomers, with the same H-capped backbone construction as the oligomers.
- Two-repeat chains, starting in open and stacked arrangements.
- Three-repeat chains, retaining two backbone diastereomers and both starting arrangements.
- Pairs of separate one-repeat segments, with parallel and antiparallel starting poses.

These are method controls rather than synthesis-shortlist admissions. The new
[route gate](../../docs/SYNTHESIS_ROUTE.md) applies prospectively to new candidate
property campaigns. BQ/AQ production downselection remains on hold.

The conformer search includes unbiased, open, parallel, slipped and antiparallel
seeds, then GFN2-xTB/ALPB(MeCN) relaxation with the geometric biases removed.
The successful initial ensembles are preserved. A diagnosed trimer reference-
embedding failure was resolved using RDKit's standard random-coordinate fallback;
earlier records and source bytes remain in the campaign's preserved attempts.

The DFT layer uses the production SMD(MeCN) geometry/energy protocol. Small
monomer controls run on reserved 32-GB-class GPUs. Parallel larger-model work
depends on an H100 energy/gradient memory pilot on the largest 149-atom case.
Its success qualifies energies and gradients, not Hessians. The finite
[coordinator](campaigns/2026-10-19_pi_stacking/coordinator.json) waits for suitable
free devices without preemption and expands onto available suitable H100s after
that dependency. It preserves failures and stops after a bounded resource wait.
The [start observation](campaigns/2026-10-19_pi_stacking/start_verification.json)
distinguishes actual resident workers from queued cases.

The first comparisons are **electronic energies**, not free energies:

- Covalent folding: `E(stacked) - E(open)` on the same exact backbone and stereochemistry.
- Separate segments: `E(pair) - 2 E(monomer)` with identical atoms and charges,
  including monomer deformation. BSSE and association entropy remain unqualified.

Independent solution-gradient and connectivity checks qualify stationary
structures. No Hessian-minimum or conformational-population claim is made at
this stage. If both starts converge to the same contact class, the report does
not invent a distinct stacked/open comparison.

The [independent structural review](interaction_review.json) checks backbone
stereochemistry, connectivity, declared charge/spin and active SP records before
energy comparisons appear in the plots. The [raw stationary collection](interaction_results.json)
is preserved separately. Warm starts for covalent oligomers are also checked for
stereochemistry before the conformer feeder publishes their DFT inputs.

The geometry report records centroid distance, angle between fitted core-plane
normals, normal-projected separation, lateral offset and closest core-atom contact.
For nonparallel cores, the normal projection is not a unique interplane distance.
The initial stack flag uses angle at most 30 degrees, normal projection 2.8–4.2 Å,
offset at most 3.5 Å and near-planar cores. These are explicit geometric heuristics,
not universal aggregation or performance thresholds. The selected starting pose
does not establish the final orientation.

Raw outputs live at `calcs/polymer_models/2026-10-19_pi_stacking/<case>/rest/`
from the project root. Input poses, xTB logs, source hashes, owners and individual
receipts stay in the campaign directory. The three-dimensional browser labels
prepared/xTB previews separately from checked DFT structures. It and the two
[maintained figures](figures/) refresh after completed electronic-structure cases.
Charged-state interactions, spin pairing, counterions, redox shifts and free-energy
sampling form the next qualification stage rather than being filled from this
neutral pilot.

## Literature context

NDI aggregation can change with redox state; a study distinguishes precipitation-
forming extended stacks from reversible short-range clusters
([Singh et al., 2023](https://doi.org/10.1002/adma.202210859)). Thus the target is
persistent or performance-limiting aggregation, rather than rejection of every
pi contact. The [rylene handcuff study](https://doi.org/10.1039/C9SC00167K) also
motivates examining multiple relative arrangements instead of a single dimer pose.
Neither aqueous or dye-model result is transferred numerically to our MeCN
polymer fragments.
