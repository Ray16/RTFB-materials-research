# Working scientific manuscript

Open [the paper](manuscript.pdf) (20 pages) and
[Supporting Information](supporting_information.pdf) (24 pages).
Editable sources are [manuscript.tex](manuscript.tex),
[supporting_information.tex](supporting_information.tex),
[references.bib](references.bib) and [novelty_references.bib](novelty_references.bib).
No Overleaf project is connected. A self-contained
[upload ZIP](overleaf_upload_2026-10-30_empty_backbones.zip) contains both documents,
bibliographies, required figures and frozen numerical TeX inputs. Both exported documents compile with
text identical to their checked local counterparts; this does not claim upload
or remote synchronization.

The [current evidence snapshot](evidence/2026-10-30_empty_backbones/manifest.json)
preserves the primary benchmark, grafting and PMDI occupancy evidence and adds
the duplicate-aware empty-backbone comparison. Both configurations have twenty
checked native inputs. Their 20/19 distinct structures give nominal 10-to-20
ensemble changes of −1.03/−1.23 kJ/mol, and both meet the declared batch criterion
across the nine rotor/duplicate settings. SI Sections 5.1 and 8.1 describe the
scientific comparison and source reconstruction. These alternate-model controls
do not replace primary-cycle components or establish global convergence.

The evidence includes 75 reconstructed native states: forty empty, 18/16 site-2
states and one separate site-4 resource point. Occupied ensembles and the complete
material assessment remain unfinished. The public raw-data DOI and final author
information also remain open; this is a working draft, not a submission-ready
paper. The new numerical macros are generated from the checked results and
frozen with their evidence, rather than hand-transcribed into the prose.

Run `manuscript/build.py` with the redox interpreter, `PYTHONPATH=src`, and native
CPU thread pools capped at one. It verifies frozen evidence and reconstructs
reported arithmetic before building both PDFs. The new checker also rejects
changed energy/coverage/configuration records, nonfinite frequencies, altered
figure bytes and an unsupported global-convergence flag. Sixty-eight focused
checks pass, including duplicate/geometry, frozen-claim and self-contained export
checks. The SI also states
the tested Nitro force-repeatability failure rather than implying production
qualification. [Compiler facts](build_check.json),
[numerical checks](claim_check.json) and the
[changed-page visual inspection](visual_review.json) remain separate.

The restored workspace is deliberately small; the original manuscript history,
bulk calculation records and earlier upload bundles remain in the archive
identified by [restoration.json](restoration.json). The current figure snapshot
is kept at `results/figures/manuscript/2026-10-30_empty_backbones/`. Its twenty small
PDFs are versioned as publication inputs, separately from runtime results.
Do not modify existing evidence snapshots; append a new snapshot for verified
scientific additions. Once an Overleaf project is connected, pull/merge before
each edit and again before pushing the checked revision.
