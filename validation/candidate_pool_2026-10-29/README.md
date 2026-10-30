# Unified candidate pool and publication revision

[Current pool pointer](current.json) · [Source inventory and counts](summary.json)
· [Recorded property coverage](property_coverage.json)
· [Independent manuscript arithmetic](manuscript_arithmetic.json)

The maintained pool is now under `results/candidates/pool/`, with a single
`current.json` pointer, a graph-deduplicated registry and byte-frozen source
inputs together in its publication directory. It contains **904 exact local
graft graphs, 921 aliases and 944 source declarations**. The baseline union had
896 graphs; a predeclared F/CN phenothiazine ring panel adds eight new graph
hypotheses. Its baseline and eight derivatives arise from 17 site recipes;
symmetry-equivalent sites are retained as recipes rather than separate products.
These derivatives have no calculated redox or grafting properties.

The registry merges the 837-graph systematic set, the historical 75-candidate
feed and one comparison control, seven later proposal sources, the original
local-model registrations, and the new ring panel. Each graph retains all source
roles, aliases, precursor/attachment declarations, recipes and proposed state
ladders. Source sections are explicit. A candidate ID reused for a different
graph, a conflicting declared product, or inherited parent properties in an
exact-graft evaluation fails assembly. Route screening is applied to every
actual declaration; conflicting proceeds/exclusions trigger review. There are
903 held graphs and one route proceed (TEMPO). Material decisions remain pending.
This is the union of the declared source inventory, not a claim to have found
all conceivable or unreviewed historical designs. Historical frozen ledgers and
manuscript cohort counts are preserved.

Seventy-six graphs link to the historical evaluation feed. Its reported-value
coverage is 64 first potentials, 63 second potentials and 60 capacity bookkeeping
values. Both solution reorganization fields are blank throughout that feed;
separate gas diagnostics do not fill them. These counts are field presence,
not raw physical revalidation or experimental validation. Existing TEMPO work
outside the feed remains in its separately reviewed archive and is not silently
promoted into these counts.

Four regression tests pass, including rejection of identity conflicts and parent
property inheritance and preservation of the new panel's O-linked linker.
Independent manuscript checks reconstruct 39 benchmark endpoints, 20
base-assisted values, finite occupancy quantities and additional recorded
controls. This is arithmetic against a frozen evidence snapshot; numerical
replay, physical QC, uncertainty and external accuracy remain separate gates.

The manuscript's title, abstract, introduction, opening discussion and conclusions
now emphasize exact-identity balanced thermodynamic comparisons across local
molecular, reaction and polymer occupancy models. Achazi et al.'s prior multistep
screening (DOI 10.1002/jcc.27299) is explicitly acknowledged; new-reference
metadata are checked against Crossref. The existing 48-reference frozen audit
is preserved, and the added reference is in a separate bibliography. The paper
and SI compile without unresolved references or citations; edited pages 1, 2
and 14 were visually inspected. The build report identifies the exact PDFs.
The main paper is 20 pages and remains a working draft. Broad conformer sampling,
full reaction accuracy, complete material criteria and public deposition remain
open. The local restoration declares no Overleaf connection and no Overleaf
remote is configured; no Overleaf synchronization was possible or claimed.

## Repeat and extend

Run from the project root with `PYTHONPATH=src` and native thread pools capped:

```bash
python scripts/candidates/expand_ring_panel.py --out NEW_PANEL.json
python scripts/candidates/merge_registry.py \
  --archive /nfs/lambda_stor_01/homes/rzhu/archive/0_redox \
  --out NEW_REGISTRY_DIRECTORY
python manuscript/check_claims.py
python manuscript/build.py
```

The first command takes `config/candidate_ring_expansion.json`; the merger takes
`config/candidate_registry_sources.json`. Use fresh output names, add new panel
sources to that inventory, and publish a new pointer only after identity, route
and input-hash checks pass. Register hypotheses before expanding. Use matched
baseline/perturbation pairs, retain exact attachment sites and all symmetry
recipes, and report mass, electronics and stereochemistry changes. Route holds
remain ahead of expensive full-property work; explicit method controls have a
separate scope. The historical evaluated feed remains an evaluation source;
`results/candidates/pool/current.json` is the unified design entry point.

A self-contained reviewer/Overleaf upload package is compiled at `calcs/validation/publication_revision_2026-10-29/reviewer_package.zip`. Its main and SI text match the checked local PDFs (20 and 22 pages); `manuscript/overleaf_export.json` records hashes and contents. It is an export only and has no remote connection. Frozen evidence bytes, including original CSV line endings, are preserved.
