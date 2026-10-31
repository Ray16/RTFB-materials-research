# Working scientific manuscript

Open [the paper](manuscript.pdf) (21 pages) and
[Supporting Information](supporting_information.pdf) (27 pages).
Editable sources are [manuscript.tex](manuscript.tex),
[supporting_information.tex](supporting_information.tex),
[references.bib](references.bib) and [novelty_references.bib](novelty_references.bib).
No Overleaf project is connected. A self-contained
[upload ZIP](overleaf_upload_2026-10-31_viologen_grafting.zip) contains both documents,
bibliographies, required figures and frozen numerical TeX inputs. Both exported documents compile with
text identical to their checked local counterparts; this does not claim upload
or remote synchronization.

The [current evidence snapshot](evidence/2026-10-31_viologen_grafting/manifest.json)
preserves the primary benchmark, grafting and PMDI occupancy evidence and adds
the duplicate-aware empty-backbone comparison. Both configurations have twenty
checked native inputs. Their 20/19 distinct structures give nominal 10-to-20
ensemble changes of −1.03/−1.23 kJ/mol, and both meet the declared batch criterion
across the nine rotor/duplicate settings. SI Sections 5.1 and 8.1 describe the
scientific comparison and source reconstruction. These alternate-model controls
do not replace primary-cycle components or establish global convergence.

A separate [portable raw-input archive](../validation/native_cache_2026-10-30/README.md)
now reproduces all 75 native-state Hessians, sampled derivative checks and nuclear
terms after extraction in a fresh NumPy-only environment, with original source
paths blocked. Its historical executable/scratch closure and public deposition
remain outside that numerical-replay scope.

The evidence includes 75 reconstructed native states: forty empty, 18/16 site-2
states and one separate site-4 resource point. Occupied ensembles and the complete
material assessment remain unfinished. The public raw-data DOI and final author
information also remain open; this is a working draft, not a submission-ready
paper. The new numerical macros are generated from the checked results and
frozen with their evidence, rather than hand-transcribed into the prose.

The new exact-candidate TEMPO comparison independently checks six MeCN
electronic calculations and their density records. Directional solution
relaxation costs are 1.477 eV for reduction and 1.639 eV for oxidation. Their
10.4% difference exceeds the declared 5% equal-curvature diagnostic, so both
directions are retained and a qualified symmetric parameter remains blank.
These composite-model descriptors are separate from grafting free energy,
electronic coupling and measured rates. The exact candidate remains pending
under the all-axis policy. Methods, numerical components and the figure are
included in the paper and SI; see the
[calculation and validation record](../validation/tempo_solution_pair_2026-10-30/README.md).

Run `manuscript/build.py` with the redox interpreter, `PYTHONPATH=src`, and native
CPU thread pools capped at one. It verifies frozen evidence and reconstructs
reported arithmetic before building both PDFs. The new checker also rejects
changed energy/coverage/configuration records, nonfinite frequencies, altered
figure bytes and an unsupported global-convergence flag. Eighty-eight focused
checks pass, including duplicate/geometry, frozen-claim and self-contained export
checks. The SI also states
the tested Nitro force-repeatability failure rather than implying production
qualification. [Compiler facts](build_check.json),
[numerical checks](claim_check.json) and the
[changed-page visual inspection](visual_review.json) remain separate.

The restored workspace is deliberately small; the original manuscript history,
bulk calculation records and earlier upload bundles remain in the archive
identified by [restoration.json](restoration.json). The current figure snapshot
is kept at `results/figures/manuscript/2026-10-31_viologen_grafting/`. Its small
PDFs are versioned as publication inputs, separately from runtime results.
Do not modify existing evidence snapshots; append a new snapshot for verified
scientific additions. Once an Overleaf project is connected, pull/merge before
each edit and again before pushing the checked revision.

The current literature pass adds the directly relevant 2025 chloromethylated-PS
phenothiazine/TEMPO study and the 2026 phthalimide-polymer kinetics/MD study.
Their primary metadata are checked during each build. These context additions
do not change the frozen numerical results or promote held candidate routes.

The exact enumerated viologen now has three independently checked native states.
The paper reports −0.608/−1.208 V versus matched Fc in MeCN and +57.91 kJ/mol
for modeled monomer disproportionation. The SI gives its minimum-admission table
and redox figure. Original source saddles were repaired before acceptance.
The maintained 904-graph identity registry is explicitly distinguished from the
frozen 878-candidate assessment. Actual-condition grafting, occupancy/coupling
and full material eligibility remain incomplete. Thirty publication-focused
tests pass; every compiled figure is now snapshot-bound and rechecked through
the end of the build.

The same exact viologen now has a checked local quaternization cycle: −11.49 kJ/mol at 298.15 K in MeCN. Three independently checked precursor seeds span 1.50 kJ/mol; rotor settings give −13.30 to −11.49 kJ/mol. The SI states the separated-ion convention and remaining salt, sampling and polymer limits. See the [reproducible calculation record](../validation/viologen_grafting_2026-10-31/README.md).
