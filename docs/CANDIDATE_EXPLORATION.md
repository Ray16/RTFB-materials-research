# Systematic candidate exploration

## Current results: 2026-10-19

The [candidate workbook](../results/candidates/latest.xlsx) contains 27 additional
exact grafts: **24 have admitted native potential pairs**. The first family
comparison panel and the OMe/Me panel have completed their computational stages;
read the [matched comparisons](../results/candidates/comparisons.csv) for current
values under the corrected thermal protocol and Fc reference. Completion of a
pipeline is separate from each field's QC. The hexyl-NDI pair remains excluded
because neither refined radical geometry passed its derivative gates.

Three phenothiazine accessible capacities (two in the additional batch) are withheld by the evidence-bound
[model-selection policy](../config/model_validation_limits.json). Their computed
potentials and nominal capacities remain available as model estimates; the checked
free-parent second-wave error is +0.547 V. Exact-graft experimental performance
and solution kinetics remain unvalidated. A larger voltage shift alone does not
establish a better battery candidate.

The source catalog and systematic pool below remain the basis for further
exploration beyond D3TaLES. The nine source-directed aryl-viologen grafts remain preparation work. The
[new cata-imide branch](../results/discovery/cata_core_expansion_2026-10-19/README.md)
now has four curated parent graphs and two prepared exact grafts, with 214 checks
and 18 verified initial conformers. All six gas UMA states passed
their warm-start contracts and their workers are terminal. The cata-NDI dianion
requires a DFT spin comparison. Both designs are now registered and their six solution-state calculations have
been submitted. Actual starts and physical admission are recorded in the
[solution campaign](../results/validation/cata_solution_2026-10-19/README.md);
no new solution potential or kinetic result is implied. Source-parent measurements retain their own conditions and
are never assigned to an exact graft.

## Earlier exploration checkpoints

The following dated preparation and execution descriptions preserve their original
counts and protocol. Their statements that calculations are pending are historical;
use the current candidate feed for completed properties.


Validation of existing properties and exploration run in parallel. The six new
motifs already computing are a pilot. The broader search starts with viologens,
PMDIs and NDIs, supported by multiple data sources and explicit enumeration.

## Prepared search space

[`config/family_search.py`](../config/family_search.py) declares three core types,
nine N-substituents, three graft anchors, and H/F/CN/OMe/Me ring variants. Ring
substitution is at one site at a time, with symmetry-equivalent graphs merged
while retaining all generation recipes. The present caps are aliphatic;
aryl-substituted viologen chemistry is an additional literature-guided branch,
not already covered by this first enumeration.

| Core | Distinct exact grafts | New relative to registry |
| --- | ---: | ---: |
| Viologen | 459 | 455 |
| PMDI | 135 | 128 |
| NDI | 243 | 236 |
| Total | 837 | 819 |

The [pool](../results/discovery/systematic_families_2026-10-18/search_space.csv)
contains molecular identities, generation recipes and bookkeeping descriptors.
No candidates were registered and no quantum jobs submitted by enumeration.
Two-electron ladders are hypotheses; each usable wave, spin state and stability
gate must be established. Exact precursor/graft synthesis remains unverified.

```bash
PYTHONPATH=src python -m redox search --out results/discovery/new_family_search
```

## First systematic comparison batch

The subsequent [12-candidate batch](../results/discovery/family_comparisons_2026-10-18/uma_completion.json)
registers four explicit contrasts per core against an unchanged ethyl control:
propyl cap length, one ring F, one ring CN, and a two-carbon ether graft linker.
F/CN variants use the same template site. Exact graphs are deduplicated against
the registry; already registered contrasts are reported rather than replaced by
an unrelated structure. The pool's 819-new count above is the original enumeration
checkpoint, before registering these twelve.

All 36 new charge states passed fixed-model UMA input, graph and convergence
checks. All nine corresponding control states also have verified current-model
records. The older ethylviologen cache lacked sufficient provenance, so its three
controls were recomputed in an isolated directory while preserving the originals.
The [gas electronic comparisons](../results/discovery/family_comparisons_2026-10-18/gas_matched_diagnostics.csv)
suggest ring CN shifts both steps positively across these cores. These differences
are preliminary gas diagnostics, not MeCN potentials or evidence of improved
battery performance. Two PMDI doubly reduced states require alternate-spin checks.

Solution DFT and native thermochemistry have begun for the selected batch; the
[execution records](../results/discovery/family_comparisons_2026-10-18/continuation_launch.json)
distinguish submitted calculations from admitted fields. Exact routes and
experimental graft properties remain unverified. The selector command is:

```bash
PYTHONPATH=src python -m redox search \
  --select-from results/discovery/systematic_families_2026-10-18 \
  --out results/discovery/new_comparison_batch
```

The saved batch contains the selected identities and recipes. These twelve are
already registered; running this command against the current registry reports
that those contrasts are no longer new. It does not register duplicate grafts.

The six earlier pilot candidates continue separately. The first completed pilot,
difluorophenothiazine, has model-admitted oxidation potentials of approximately
0.708 and 1.715 V versus Fc. Its result does not establish stability or kinetics
for the other substituents or families.

All twelve selected family pipelines have now started. The first complete
systematic ladder, the ring-CN viologen, gives approximately -0.297/-0.973 V
versus Fc under the native model. Its ethyl control is receiving matched native
thermochemistry: the older GFN2 correction cannot support a direct shift against
these new values. Partial ladders retain blanks rather than inferred second waves.

## Complementary sources

The next [six-candidate ring panel](../results/discovery/donor_ring_comparisons_2026-10-18/README.md)
adds OMe and Me to all three cores at the same site as F/CN. The exact ethyl
controls and graft anchors are unchanged. Eighteen states were registered and
built; UMA gas preoptimization started on all 12 free qualified T4 devices.
This is a mechanistic contrast motivated by the initial potential/spacing shifts,
not a prediction of improvement. Native solution calculations remain pending.
The selector now rejects an already registered ring site rather than silently
substituting another positional isomer. Named comparison panels are declared in
`config/family_search.py`.

The machine-readable [source catalog](../config/discovery_sources.json) records
scope and admission status. Raw downloads, repository revisions and checksums
are preserved in `data/raw/discovery_sources/2026-10-18/`.

| Source | Intended contribution | Important boundary and current status |
| --- | --- | --- |
| [Gao et al., 2025](https://doi.org/10.1016/j.jpowsour.2024.236035) | Experimental compilation with 336 nonaqueous quinone first reductions in DMF, plus 195 aqueous records | [Institutional full text](https://engineering.lehigh.edu/sites/engineering.lehigh.edu/files/_DEPARTMENTS/ise/pdf/tech-papers/25/25T_012.pdf) archived. Section 2.2.3 traces the nonaqueous subset to Prince et al. 2022; it does not fill the viologen/imide gap. Supplementary rows remain unavailable. |
| [ReSolved](https://github.com/grynova-ccc/ReSolved) | Downloaded 19,785 computed molecular records with five solvent targets, including acetonitrile | All supplied SMILES have formal charge zero. Useful for neutral compounds, not direct coverage of dicationic viologen reductions. Absolute source predictions remain separate from project Fc-relative values. |
| [RedDB](https://www.nature.com/articles/s41597-022-01832-2) | Additional structures; downloaded v2 table with 15,882 reaction records | Aqueous two-proton/two-electron chemistry and predicted aqueous solubility require their own scope. They cannot replace MeCN sequential electron-transfer labels. |
| [Aryl-viologen study](https://doi.org/10.1021/acsaem.5c02255), [raw data](https://zenodo.org/records/17671488) | Seven named parent structures curated; nine source-directed O-grafted hypotheses prepared | Aqueous pH14/AgAgCl data remain separate from MeCN/Fc. Raw archive CC-BY-4.0; supporting-information record CC-BY-NC-4.0. Exact-graft synthesis, speciation and performance remain unverified. |
| [NDI ion-pairing study](https://doi.org/10.1021/jacs.3c11176) | Five parent components, nine salts, 12 MeCN formal potentials, nine solubilities and seven diffusion/rate records | The earlier `.s003` record has coordinates. The newly retrieved [main SI, `.s001`](https://acs.figshare.com/articles/25208990), contains measured tables. [Curation and 82 checks](../results/discovery/ndi_electrolyte_2026-10-18/README.md) preserve electrolyte/aggregation conditions; these are not bare-ion calibration labels. |
| [Cata diimides](https://doi.org/10.3762/bjoc.20.155) | Four compounds and eight half-wave potentials; alternative naphthalene/anthracene diimide cores | [Source rows](../data/curated/discovery/cata_diimides_2024.json) are in CH2Cl2/Fc, distinct from MeCN. [Four source graphs](../data/curated/discovery/cata_diimides_2024_graphs.json) are now curated and two graft hypotheses registered. Original literature comparison rows were not counted again. |
| [ORedOx159](https://github.com/ANRMoMoPlasm/ORedOx159) | Additional computational benchmark structures and energies | Computational PBE0 reference; not an experimental battery dataset. Repository inspected, data not imported. |

The earlier [coverage audit](../results/discovery/systematic_families_2026-10-18/source_coverage.json)
used a broad two-imide-group rule. The corrected
[core-scope audit](../results/discovery/source_quality_2026-10-18/imide_scope_coverage.json)
recognizes two PMDI/NDI cores and one other multiple-imide structure in D3TaLES,
and no viologens. RedDB's 253 apparent bis-imide hits are separate imide units,
including N–N-linked phthalimides; none match the PMDI/NDI core patterns. Only
12 of those structures occur in reactions with both optimizer grades zero.
ReSolved has no recognized viologens or PMDI/NDI cores. More rows alone therefore
do not solve the family imbalance. These are pattern recognition results, not
proof that every other redox motif is absent.

The new NDI condition series is one original study, not twelve independent
studies. Its whole ammonium-functionalized molecule has charge +2/+1/0 across
the reductions, while its NDI core is described as 0/−1/−2. For the ethyl parent,
replacing 0.5 M K(TFSI) with 0.5 M Li(TFSI) shifts the first/second formal
potentials by +0.13/+0.40 V and narrows the spacing by 0.27 V. Explicit ion pairing
and aggregation therefore belong in a separate validation branch; a bare-core
prediction must not be fitted to these values. The separate cata-imide study
motivates exploring core topology in addition to the current cap/linker/ring axes.

Discovery retains unfamiliar multiple-imide structures for review. It gives them
no PMDI/NDI wave-spacing evidence, established electron count or nominal capacity.
For registered molecules, nominal capacity instead follows the explicitly
declared charge ladder; usable capacity still has separate physical gates.

Source quality also requires inspecting the fields: ReSolved has 1,653 rows with
`same_sml` marked `no` or `idk`; their meaning needs resolution before adopting
labels. [RedDB's source definitions](https://www.nature.com/articles/s41597-022-01832-2)
identify zero as the best optimizer grade and four as nonoptimal geometry. The
[quality audit](../results/discovery/source_quality_2026-10-18/summary.json)
finds 4,503 reaction rows with both grades zero, 999 involving grade four, and
10,380 intermediate-grade rows held for review. A grade is not a Boolean
convergence flag; even the best source grade does not establish prediction accuracy.
No downloaded labels have been admitted to project calibration or substituted
for project calculations.

```bash
PYTHONPATH=src python -m redox source-audit \
  --source data/raw/discovery_sources/2026-10-18 \
  --out results/discovery/new_source_audit
```

Gao's reported table count of 19 core structures is retained as source metadata,
not evidence of 19 different battery families. Its SHE-converted DMF values also
require their original electrode provenance. They cannot be pooled directly with
MeCN/Fc labels; deduplicate against the underlying Prince study before validation.

## Selection and computation

The [new diaryl-viologen branch](../results/discovery/aryl_viologen_2026-10-18/README.md)
retains both N–aryl bonds and explores O-attachment through carboxyl or phenol
handles. This is a source-directed design choice, not a transfer of the parent's
alkaline performance. The isolated acid/phenol graphs and explicit dichloride
salts are preserved separately; pH14 speciation is not inferred. All nine exact
mono-O-grafts regenerate through existing rules. Only nominal descriptors are
computed: no registration or QM has occurred. Polyfunctional precursors require
site-selectivity/monoattachment review, and the 70-atom maximum requires a new
UMA memory pilot before screening beyond the earlier qualification scope.

1. Curate structures and condition metadata from each source, preserving original
   values and provenance. Do not silently discard salt components; curating the
   active ion is a separate recorded step.
2. Deduplicate exact structures and redox events across sources. Separate parent
   evidence from exact grafts, and group validation splits by original study and
   molecular core to reduce overlap between training and evaluation.
3. Select a batch stratified by family, cap, linker and ring substitution, retaining
   matched comparisons. Raw pool size must not give one family all compute slots.
   Nominal capacity and SA are proxies, not sufficient electrochemical objectives.
4. Apply conformer and charge/spin-aware cheap screening. Qualify runtime/memory
   for each new scope before extending the current UMA resource qualification.
5. Compute solution DFT, active SPs and native thermochemistry for the selected
   batch. Rank admitted potentials, accessible electrons, integrity, stability and
   reorganization on their separate gates; keep unavailable fields blank.

The first 837-structure pool is an explicit starting search space, not an exhaustive
chemical universe or a proposal to run all structures through expensive DFT.
Quinones and phenothiazines remain supported additional families. Dataset and
literature curation should guide further core, linker and substituent extensions.
