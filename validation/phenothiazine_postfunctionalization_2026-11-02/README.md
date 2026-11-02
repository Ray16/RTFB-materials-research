# Direct N-linked phenothiazine: route and pool expansion

The unified registry now contains 912 exact local-graft graphs, preserving every
prior graph, alias and reported evaluation. Eight new F/CN-substituted designs
share a direct N-benzylic attachment with the existing phenothiazine baseline.
All new electrochemical properties remain pending; no quantum jobs were launched.
The [preservation check](registry_preservation.json) verifies all 904 previous
graphs and all previous evaluations unchanged.

A newly checked primary source supports further computation of the unsubstituted
baseline: Wright and Jin, *Polymer-supported cation radicals*, J. Org. Chem. 1989,
54, 965–968, [DOI 10.1021/jo00265a043](https://doi.org/10.1021/jo00265a043).
Scheme I and the experimental section report N-attachment of phenothiazine to
chloromethylated polystyrene beads using NaH in DMF at ambient temperature for
six hours. General procedures identify commercial phenothiazine. The washed
polymer's nitrogen analysis supports the reported product assignment. Subsequent
chemical oxidation and reduction support retained redox function. The paper
also reports that K2CO3 did not modify this polymer with phenothiazine; that
negative result must not be replaced by a general assumption about bases.

The source is a bead-attachment precedent for this unsubstituted connection.
Our 4-methylbenzyl molecule is a local model, and our actual crosslinking,
accessibility, loading and conversion are unknown. The paper's reported yield
and composition are not measurements of our beads. Chemical redox retention
in the source is not proof of MeCN battery accuracy, cycling, transport or weak
operating-state interactions. A separate primary polymer comparison, Golriz et
al., [DOI 10.1039/C4RA17107A](https://doi.org/10.1039/C4RA17107A), reports poorer
stability for its short-spacer phenothiazine polymer; its monomer-first synthesis
is separate from our route and its measured potentials are not imported.

The [reviewed exact route](reviewed_route.json) passes the maintained prospective
gate for the unsubstituted baseline. Its [NaH stoichiometry](condition_stoichiometry.json)
is atom- and charge-balanced with NaCl and H2 as bookkeeping products. Salt
phases/activities and hydrogen chemical potential have not been calculated;
there is no actual-condition reaction G here. The generator's separate HCl
reference equation remains explicitly labeled as a reference, not the reported
NaH synthesis.

The eight substituted precursors remain held for access, selectivity, condition
transfer and redox-retention review. Neither the direct baseline's route nor an
O-linked hydroxyethyl route is inherited. The [nine-row gate panel](route_decisions.json)
includes the existing baseline plus eight held new graphs. The full 912-graph
registry has three route passes and 909 holds, with zero full material admissions.
Thus this is an additional route-supported computation priority, not a newly
validated battery material or an established descriptor leader.

The [source review](source_review.json) records original article page scans and
publisher DOI metadata, with durable shared-storage paths and hashes. The
publisher PDF request returned 403; the publicly available original page scans
were consulted directly rather than accepting a third-party abstract. The
2025 styrofoam paper is additional post-functionalization context; its inaccessible
SI is not used to assert experimental conditions. Original scan images and
calculation outputs stay on shared storage.

Reproduce the design panel with the maintained declarative generator:

```bash
python scripts/candidates/expand_ring_panel.py \
  --spec config/candidate_ring_expansion_direct_ptz.json \
  --out <new-proposal-file.json>
```

Each generated design contains its exact precursor/product pair, mapped new
bond, retained symmetry-equivalent recipes and one-electron state hypothesis.
The shared attachment validator is also used by the prospective route gate.
The maintained merger reads `config/candidate_registry_sources.json` and writes
a fresh union inside `results/candidates/pool/`; `current.json` is its stable
pointer. The older 904-graph manuscript snapshot remains a distinct frozen
assessment until its next evidence update.

The [panel checks](panel_checks.xml) pass 18 tests, including exact mapping,
registry preservation logic and existing route-stage controls. The [route
checks](route_checks.xml) pass seven tests covering the expanded mapped panel,
the NaH stoichiometry and rejection of transferring the baseline's review to a
ring-substituted or ether-linked structure. These overlap in shared controls;
they are not additional quantum calculations or experimental validation.
