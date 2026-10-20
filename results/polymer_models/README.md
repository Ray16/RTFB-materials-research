# Polymer-fragment and interaction models

All result categories: [results index](../README.md).

**Open the [rotatable 3D interaction report](index.html)** for the current pilot.
[Interaction progress](interaction_status.json) and [reviewed results](interaction_checked.json)
are separate from the [559-graph inventory](current.csv) and its [status](status.json).

**Screening objective (user clarification, 2026-10-20): prefer weakly interacting
redox sites and polymer segments.** Strong persistent pi stacking and undesired
intersite electron transfer are design risks. A more negative association or
folding energy is a warning to investigate, not a selection benefit. Track
stacking propensity separately from electronic coupling and electron-transfer
rates; the existing neutral electronic-energy pilot has not computed those rates.
No numerical interaction cutoff is validated. Preserve the current results as
method controls, without treating stacking as a candidate advantage.

The current finite DFT pilot attempted 27 cases: 17 structurally reviewed
stationary structures and 10 preserved failures, with no active queue claims.
It includes monomers, covalent dimers and trimers, and two separate-segment
complexes per motif. No tetramers were prepared or computed in this pilot.
Interaction free energies, charged-state coupling and redox shifts remain unfinished.
The [bounded geometry recovery](campaigns/2026-10-20_geometry_recovery/plan.json)
retries only those ten unconverged cases; original outputs remain unchanged.
Its current results join the original slots only after receipt and structural QC.
Follow [interaction status](interaction_status.json) for actual recovery progress.

## Reproducible pipeline and checkpoints

| Stage | Versioned implementation | Recorded outputs and checks |
| --- | --- | --- |
| Build H-capped polystyrene chains and attach exact motifs | `src/redox/build/oligomers.py` | Repeat count, site occupancy, exact SMILES and backbone stereoisomers; current pilot uses 1–3 repeats |
| Generate open/stacked poses and conformer ensembles | `scripts/pipeline/prepare_pi_stacking.py`, `src/redox/build/polymer_contacts.py` | Unbiased/open/parallel/slipped/antiparallel seeds, then unconstrained GFN2-xTB/ALPB(MeCN); ring-core atom mapping retained |
| Stage immutable warm starts | `scripts/pipeline/stage_polymer_conformers.py` | Connectivity, stereochemistry, input XYZ and source hashes; separate-segment pairs retain two intact molecular graphs |
| Optimize and evaluate energy | `scripts/pipeline/run_pi_stacking.py`, `src/redox/qm/dft.py`, `src/redox/qm/sp.py` | r2SCAN-D4/SVP/SMD(MeCN) geometry; active wB97M-V/TZVPD/SMD electronic energy; independent solution gradient |
| Save optimization progress | `src/redox/workflows/geometry_checkpoints.py` | Per-evaluation XYZ trajectory, energies, gradients, SCF convergence and latest coordinates on NFS; trial evaluations are explicitly unqualified |
| Audit physical geometry | `scripts/validation/audit_polymer_geometries.py`, `scripts/validation/review_polymer_interactions.py` | Composition, bond-length plausibility, severe clashes, full connectivity, stereochemistry, core contacts and active SP provenance |
| Recover nonconvergence | `scripts/pipeline/recover_pi_stacking.py` | Separate frozen ten-case campaign; 400-step budget and tighter geometry SCF (1e-11), unchanged physical acceptance criteria; optimizer-returned coordinates saved even without convergence |
| Collect and compare | `scripts/pipeline/collect_pi_stacking.py` | Reuse verified successes; admit only identity-matched recovery receipts; same-composition folding and atom-balanced segment association remain distinct |
| Present results | `scripts/plotting/polymer_models/update_interaction_figures.py` | Maintained 3D viewer and energy comparisons; weak interactions are the design objective |

The [geometry audit](../reviews/2026-10-20/polymer_geometry_recovery/geometry_audit.json)
checked all 27 starts and 17 completed DFT geometries without flagged sanity issues.
All ten original failures exhausted 100 optimizer steps. Their final temporary
coordinates were not saved by the earlier runner and were no longer recoverable;
the new attempts therefore use the preserved, checked original seeds. The driver
now checkpoints in future runs as well as in recovery. Checkpoints preserve
restart coordinates; convergence, topology, spin and physical QC still control
whether a result can be used. No Hessian, conformer-population or rate claim is
made from an unconverged checkpoint.

The recovery raw root is `calcs/polymer_models/2026-10-20_geometry_recovery/`.
Each case has `rest/geometry_checkpoints/evaluated_trajectory.xyz`,
`evaluations.jsonl`, `last_evaluated.xyz` and `latest.json`; after optimizer return,
`optimizer_returned.xyz` and its convergence record are retained. The NFS campaign
folder owns queue claims, failed attempts, logs, ownership and completion receipts.
The [recovery registry](recovery_campaigns.json) binds the frozen plan; the collector
refuses a changed solvent/method, a different stereoisomer or replacement of an
already verified original case.

Read-only review and collection (these commands submit no GPU calculations):

```bash
PYTHONPATH=src python scripts/validation/audit_polymer_geometries.py --out /tmp/polymer_geometry_audit.json
PYTHONPATH=src python scripts/pipeline/collect_pi_stacking.py
```

The campaign preparation, launch and worker entrypoints are deliberately separate.
Do not rerun an existing preparation or clear a failed attempt in place. Inspect
its ownership, input hashes, last checkpoint and error before planning a bounded
continuation. Launch recovery through `recover_pi_stacking.py launch --role queue_large`
only after node preflight; the entrypoint uses the system GPU reservation gate.

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
