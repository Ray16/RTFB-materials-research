# Exact baseline viologen route review

The baseline N-methyl viologen design now passes the prospective attachment
gate. The source is Kim, Sanford, Vaid and McNeil, *Chemistry—A European Journal*
2022, DOI [10.1002/chem.202200149](https://doi.org/10.1002/chem.202200149),
including its primary Supporting Information. This is a route admission for
property calculations, not a material-feasibility decision.

SI S3–S4 establishes access to the mono-methyl bipyridinium iodide precursor.
SI S5 demonstrates its reaction with preformed Merrifield resin in refluxing
MeCN for 48 hours. The available pyridyl nitrogen is the attachment site;
the other nitrogen is already methyl-quaternized. SI S6 and the main article
report separate counterion exchange and retained viologen electroactivity.
The rendered precursor and bead-attachment drawings were checked against the
declared graphs. The source's 75% bead conversion is an approximate assumption
with inconsistent elemental-analysis constraints, not an exact measured yield.

The local product is `Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1`, InChIKey
`PBNXFHGCSKYNHC-UHFFFAOYSA-N`. Its 4-methylbenzyl cap models the polymer site;
the exact capped molecule is not the demonstrated support. The atom-mapped
new bond and neutral pyridyl precursor are reconstructed by the pipeline, not
asserted by candidate ID. The balanced model reaction has charges
`(+1) + 0 -> (+2) + (-1)` and releases chloride. It requires neither an HCl
byproduct nor the proton-capture/base offset used for N–H/O–H substitution.
Precursor iodide cancels only under a separated spectator-ion convention;
ion association and the subsequent PF6 exchange are separate unresolved terms.

Both `search_viologen_b7a4dca7a5e1` and `viologen` identify this one exact graph.
[prospective_gate.json](prospective_gate.json) contains their actual gate runs.
[registry_comparison.json](registry_comparison.json) verifies that only this
graph's route status/decisions change in the 904-graph registry. All identities,
recipes, declarations and property records are retained. Route counts become
two proceed and 902 hold. Other substituted viologens remain unreviewed.
No quantum calculation or old thermochemistry is qualified by this route review.
In particular, the archived viologen records have no checked native Hessians;
they cannot supply the new full-property campaign's free energies.

The maintained publisher accepts an explicit updated route catalog and snapshots
its bytes. It checks that neither the catalog nor the assessment changes during
assembly. The archived catalog remains untouched. Sixteen focused tests include
the actual baseline pass, charge/atom balance, rejection of an ethyl-cap transfer
and an unreviewed alias, an archive-preserving registry publication and a change during snapshot creation.

[primary_source_review.json](primary_source_review.json) records source sections
and file hashes. The primary article XML, original SI PDF and retrieval receipts
are stored privately under `calcs/literature/2026-10-30_viologen_route/`; these
third-party source files are not redistributed as project results. The current
full registry and its input snapshot are in
`results/candidates/pool/2026-10-30_viologen_checked/`. The paper's existing frozen
all-axis cohort remains a separate historical assessment. The new route evidence
is prepared for the next scientific manuscript revision; complete end-to-end
candidate evaluation under M1 remains open.
