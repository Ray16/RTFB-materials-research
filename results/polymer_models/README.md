# Polymer-fragment and interaction models

All result categories: [results index](../README.md).

**Stable inventory:** [current.csv](current.csv). **Progress:** [status.json](status.json).

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
model. The intended meaning of “two beads” is being clarified with the user.

Redox-site charge distributions, spin states, counterions, chain conformation and
stereochemistry must be tracked. First and second molecular redox waves are not
interchangeable with the first and second electrons added to a multi-site cluster.
Statistical site degeneracy must be separated from genuine electronic/cooperative
splitting before interpreting oligomer potential shifts. Charged viologen oligomers
require an explicit counterion/ionic-environment treatment before physical ranking.

This is a structure inventory; quantum calculations and interaction values remain
pending. The initial representative panel is recorded in the status file.
