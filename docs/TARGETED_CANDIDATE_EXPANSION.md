# Focused expansion toward defensible candidate decisions

The existing 837-design library gives systematic coverage within viologen, PMDI
and NDI chemistry. Its size is adequate for a focused initial screen. The next
expansion should fill chemical and polymer-design gaps while connecting every
structure to the same all-axis evaluator.

The predeclared [comparison plan](../config/targeted_expansion_2026-10-23.json)
contains **34 predeclared design slots**, originally specified before exact-graph
construction and route review. A later [PMDI matrix](../results/discovery/pmdi_linker_matrix_2026-10-26/README.md)
prepares four exact graph hypotheses, two reused and two new; shielding and
route evidence remain unresolved, and formal coverage of the 34 slots is not
yet established. Existing graphs are not counted as new molecules.

Separately, one [AQ-2-carboxylate local graft proposal](../results/discovery/aq_carboxylate_2026-10-26/README.md)
has now been built and screened against the route gate. It has a primary-source
polymer attachment precedent, an exact precursor/product graph and a mapped bond;
its route is held and it has no new redox calculations. It is not assigned to a
shielding/linker comparison slot until a matched panel is defined and reviewed.

A separate [amino-TEMPO local graft proposal](../results/discovery/tempo_amino_2026-10-26/README.md) compares an N-linked attachment handle with the existing O-linked model. Its route is held; it is not a matched shielding or linker-rigidity comparison, and has no calculated properties.

A [source-reported hydroxyethyl PMDI precursor](../results/discovery/pmdi_hydroxyethyl_source_control_2026-10-26/README.md) provides a separate exact O-graft proposal. It retains a competing imide N–H, and selective attachment remains unverified; it does not validate the ethyl-capped linker matrix precursors.

An existing oxidation-active [TEMPO local graft](../results/grafting/tempo_exact_support_2026-10-26/README.md)
now has an exact-support attachment precedent and a checked local-model
oxidation potential. Its route screen passes, but its material decision remains
pending. It is reused historical chemistry, not a newly enumerated product or a
completed shielding/linker matrix.

| Component | Design question | Planned coverage |
| --- | --- | --- |
| Core breadth | Does the route-aware screen behave consistently across different redox chemistry? | Viologen, PMDI, NDI, BQ, AQ, phenothiazine and TEMPO |
| Shielding × linker panel | Can a design reduce unwanted contacts/coupling while retaining accessible grafting and useful capacity? | Baseline/shielded × flexible/constrained linker, four slots per core: 28 total |
| Complementary electronics | Do donor/withdrawer perturbations alter the usable redox path and stability without breaking attachment? | Two additional electronic perturbations for BQ, AQ and phenothiazine: six slots |
| Positive attachment control | Can the workflow reconstruct a literature-demonstrated support functionalization? | Ferrocene bead attachment, pending exact model/source identity review |

Every slot needs an exact precursor, product, attachment atom/new bond, declared
state path, baseline and hypothesis. Ring position and precursor selectivity
must be explicit. Preserve all symmetry-equivalent generation recipes, but count
one exact graph once. A failed route slot is reported rather than silently
replaced by a more convenient positional isomer.

## Chemical constraints

BQ/AQ attachment must preserve both quinone carbonyls; grafting a required redox
oxygen into an irreversible ether would change the chemistry. TEMPO designs
must preserve the N–O radical and use a distinct supported handle. One accessible
electron may be appropriate for a nitroxide or ferrocene control; no family gets
two useful electrons merely from its name or a proposed state ladder.

Constrained and flexible linkers should retain attachment chemistry and similar
reach where possible. Record residual formula, mass, polarity and conjugation
changes. Rigidification may favor persistent stacking, while shielding may hinder
attachment or reduce gravimetric capacity. Those are testable competing effects,
not assumed improvements.

Prior demonstrations establish the relevance of ferrocene/viologen bead
functionalization and of substituent-dependent phenothiazine oxidation stability.
They do not validate every proposed derivative or medium. Sources:
[Kim et al.](https://doi.org/10.1002/chem.202200149) and
[Kowalski et al.](https://doi.org/10.1039/C7TA05883G).

## Evaluation order

Cheap graph/precursor/route triage precedes new expensive redox campaigns. A
documented route can proceed to the next stated calculation stage; it is not a
full material endorsement. Selected exact cases then receive state/path and
kinetic checks, defined-route G, and conditional polymer environment assessment.
Every required axis must meet its acceptance criterion before shortlisting.

Use a few predeclared representative cases for complete end-to-end assessment,
retaining negative and pending results. The manuscript should explain which
gates affect decisions and which conclusions change with model assumptions.
Its claim is a focused, route-aware study, not comprehensive coverage of all
redox chemistry. A larger enumerated count cannot substitute for those complete
case studies.
