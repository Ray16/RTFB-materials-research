# Exact AQ-carboxylate graft proposal

[Mapped candidate](candidates_mapped.json) · [Route screen](route_screen.json) ·
[Identity source](pubchem_identity.json) · [Checks](summary.json) · [Tests](tests.xml)

This proposes one capped-site product formed by linking anthraquinone-2-carboxylic
acid through its carboxylate oxygen to the benzyl carbon of a 4-methylbenzyl
chloride model. The acid identity matches PubChem CID 67030. A reported attachment
of that acid to soluble poly(4-chloromethylstyrene) supplies a close polymer
precedent ([Oyaizu and coauthors](https://doi.org/10.1002/pat.1968); [thesis](https://waseda.repo.nii.ac.jp/record/19599/files/Honbun-5904.pdf)).
The capped product is a model of the connection, not an isolated product or a
demonstration of attachment to crosslinked beads.

The two quinone carbonyls are retained in the graph, and the mapped heavy atoms
identify the new O–C bond. The neutral-reference reaction is atom and charge
balanced, but its HCl byproduct does not represent the actual DBU-assisted route.
Solvent, counterion and base speciation, selectivity, redox-core stability,
bead accessibility and electrochemical properties remain to be evaluated.
The proposed neutral and singly reduced charge/spin states are hypotheses,
not computed state assignments.

The prospective route gate returns **hold**: no full redox campaign is admitted.
This one proposal does not complete the predeclared shielding-by-linker matrix,
and is not a material shortlist. No experimentally measured conversion or
candidate performance is inferred from the graph or the literature precedent.
The pre-mapping prototype is retained as `candidates.json` for traceability.

The route registry has since been extended with TEMPO. This snapshot's original
registry hash is preserved and reconciled in the [append-only registry audit](../../grafting/tempo_exact_support_2026-10-26/registry_reconciliation.json);
the AQ entry itself is unchanged.
