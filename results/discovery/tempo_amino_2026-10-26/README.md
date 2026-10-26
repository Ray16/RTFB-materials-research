# Proposed N-linked amino-TEMPO model

[Exact mapped graph](candidates.json) · [Source identity](pubchem_identity.json) · [Route screen](route_screen.json) · [Checks](summary.json) · [Registry reconciliation](registry_reconciliation.json) · [Focused tests](tests.xml) · [Regression](tests_combined.xml)

This is a proposed 4-amino-TEMPO attachment through its exocyclic nitrogen to a 4-methylbenzyl chloride site model. [PubChem CID 550942](https://pubchem.ncbi.nlm.nih.gov/compound/550942) identifies the precursor. The stored atom mapping specifies the benzyl C–N bond. The reference reaction balances atoms and charge, and the precursor and product graphs each retain one nitroxide radical.

The earlier [4-hydroxy-TEMPO route](../../grafting/tempo_exact_support_2026-10-26/README.md) is an O-linked support precedent, not experimental evidence for this N-linked graft. This is a comparison of attachment handles, not a controlled shielding or linker-rigidity comparison. The neutral HCl-reference reaction is stoichiometric bookkeeping; it does not specify the actual reaction conditions.

The prospective route gate returns **hold**. Precursor procurement, conditions, overalkylation and selectivity, radical retention under reaction conditions, actual-route thermodynamics, and electrochemical and polymer properties are unverified. No new quantum jobs were submitted, no full-axis selection is implied, and this candidate is not included in the previously frozen 879-graph ledger.

Rebuild with `PYTHONPATH=src python scripts/pipeline/prepare_amino_tempo_candidate.py --identity results/discovery/tempo_amino_2026-10-26/pubchem_identity.json --out NEW_CANDIDATE_FILE.json`, followed by `scripts/pipeline/screen_grafting_routes.py` using the current route registry.
