# Source-reported hydroxyethyl phenothiazine precursor

[Precursor and proposed graft](../../figures/candidates/ptz_hydroxyethyl_source.png) · [Vector figure](../../figures/candidates/ptz_hydroxyethyl_source.svg) · [Mapped candidate](candidates.json) · [Source review](source_review.json) · [Route gate](route_screen.json) · [Material decision](all_axis_decision.json) · [Checks](summary.json) · [Tests](tests_verified.xml)

[Persson’s 2012 thesis, printed page 7](https://api.mountainscholar.org/server/api/core/bitstreams/afbf9593-a8d1-4410-9c5a-c9f0b9fe441a/content), reports preparation of 2-(10H-phenothiazin-10-yl)ethanol and proton NMR for the precursor. [PubChem CID 352306](https://pubchem.ncbi.nlm.nih.gov/compound/352306) matches the represented C₁₄H₁₃NOS graph. This is one source report; it does not establish independent replication or present procurement.

The proposed local model attaches the alcohol oxygen to a benzyl carbon. Its atom-mapped HCl-reference reaction is balanced, but that bookkeeping does not define actual bead-reaction conditions. The related [TEMPO bead attachment](../../grafting/tempo_exact_support_2026-10-26/README.md) concerns a different precursor and does not prove this graft’s selectivity or redox retention. This O-linked product is a different graph from the existing N-benzyl phenothiazines.

The route remains **hold**, and the shared all-axis decision remains **pending**. Grafting conditions, selectivity, actual-route thermodynamics, exact-graft redox properties, association, coupling and accessibility are unverified. One oxidation state change is proposed as a computational hypothesis; reversibility and usable capacity are not inferred. No quantum jobs were submitted, and this exact product is absent from the frozen 879-graph ledger.

Reproduce from the project root with fresh output names:

```bash
export PYTHONPATH="$PWD/src"
python scripts/pipeline/prepare_ptz_hydroxyethyl_candidate.py --identity results/discovery/ptz_hydroxyethyl_source_2026-10-26/pubchem_identity.json --out NEW_CANDIDATE.json
python scripts/pipeline/screen_grafting_routes.py --candidates NEW_CANDIDATE.json --out NEW_SCREEN.json
python scripts/pipeline/assess_source_route_candidate.py --candidate NEW_CANDIDATE.json --gate NEW_SCREEN.json --policy config/candidate_feasibility.json --out NEW_DECISION.json
```

The [plotting script](../../../scripts/plotting/candidates/plot_source_precursor_route.py) draws the stored mapped product and attachment bond.
