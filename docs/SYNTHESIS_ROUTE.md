# Current synthesis route

The user corrected the route on 19 October 2026: **functionalize an existing
polymer with redox groups**. For the Merrifield-type system, the working picture
is attachment through reactive benzyl sites on the preformed polymer.

The existing 4-methylbenzyl-capped structures represent local attachment
chemistry. They remain molecular models: polymer constraints, grafting density,
swelling, accessibility, counterions and charge transport require separate
assessment. Molecular nominal capacity does not establish usable polymer capacity.

Candidate downselection should combine redox properties with precursor and
attachment feasibility, attainable loading, accessibility of grafted sites and
redox-center retention. Inhibition during vinyl polymerization is no longer the
primary route constraint.

Grafting feasibility is a selection gate. If a balanced reaction remains strongly
unfavorable under defined experimental conditions, that route should not enter
the synthesis shortlist regardless of capacity or potential. Positive free
energies for the current neutral HCl-producing reference do not yet quantify a
base-assisted experimental reaction. The completed DBU/DIPEA pilot gives favorable
modeled reactions for both highlighted BQ designs and AQ-1-O: **base required;
retained for evaluation**. AQ-2-O/6-OMe has unresolved precursor QC. Attachment
selectivity and yield remain untested. Their computed redox descriptors remain
available, separately labeled. A near-zero reference
value for another candidate does not establish a feasible or selective route.
The current decisions are in
[the synthesis assessment](../results/grafting/selection_assessment.json).

For primary-source attachment precedents, their actual conditions and the limits
of transfer to our structures, see [practical grafting](PRACTICAL_GRAFTING.md).

## Prospective evaluation order

1. Check precursor access and a credible attachment reaction before expensive
   property calculations. Review the actual reactive form, solvent, base,
   counterions, steric access, site selectivity and preservation of the redox core.
2. Assess the balanced grafting reaction under those conditions. Reaction free
   energy and kinetic ease are distinct; use known chemistry or targeted barrier
   calculations where needed. Assign proceed, hold or exclude with evidence.
3. Advance candidates with reviewed workable routes to the full redox-property
   campaign. Retain unsuccessful designs and their reasons instead of silently
   removing them from the record.
4. Compare monomer, dimer and trimer models for persistent stacking and changes
   in redox behavior. Include open/stacked conformers, backbone stereochemistry
   and separate-segment contacts. Charged states need their own spin and ionic-
   environment models.

The oligomer selection objective is **weak interaction between redox sites and
polymer segments** (user clarification, 2026-10-20). Persistent strong pi stacking
and undesired intersite electron transfer are risks. Strongly negative modeled
association energies are concerns to investigate, not benefits to maximize.
Geometric contacts and electronic energies do not by themselves establish an
electron-transfer rate; coupling and charged-state behavior require separate
evaluation. Whole-bead contacts remain beyond the current molecular models.

Post-grafting loading and spacing may be heterogeneous and remain unmeasured.
Separate parent reactive-handle density from graft conversion; compare occupancy
and conformational scenarios without assigning experimental probabilities. Current
fully occupied oligomers are dense local controls. An unoccupied modeled benzyl-
chloride site is distinct from a repeat unit lacking a chloromethyl handle.
Use matched, padded-backbone controls to assess adjacent and separated motifs;
sequence spacing does not determine their three-dimensional separation.

Run `scripts/pipeline/screen_grafting_routes.py` to produce the prospective
[route screen](../results/grafting/route_screen.json) from the evidence file.
At its initial checkpoint no exact-candidate route is automatically certified:
the family-level literature sources have not yet been converted into the required
candidate-specific review. Existing 75-candidate results are retained. Explicit
method controls can be modeled to establish the interaction protocol without
representing a synthesis-priority decision.

The earlier vinyl-monomer branch is preserved under
`results/discovery/vinyl_monomers/`. Its numerical results belong to those explicit
structures. No result is transferred between source parents, capped graft-site
models, vinyl structures or polymers.
