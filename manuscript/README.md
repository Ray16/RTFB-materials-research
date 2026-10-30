# Working scientific manuscript

Open [the paper](manuscript.pdf) (20 pages) and
[Supporting Information](supporting_information.pdf) (24 pages).
Editable sources are [manuscript.tex](manuscript.tex),
[supporting_information.tex](supporting_information.tex),
[references.bib](references.bib) and [novelty_references.bib](novelty_references.bib).
No Overleaf project is connected. A self-contained
[upload ZIP](overleaf_upload_2026-10-30_polymer_replay.zip) contains both documents,
bibliographies and all required figure PDFs. Both exported documents compile with
text identical to their checked local counterparts; this does not claim upload
or remote synchronization.

The [current evidence snapshot](evidence/2026-10-30_polymer_replay/manifest.json)
preserves the preceding numerical evidence and adds the sixty-state polymer
physical replay and its finite sampling analysis. The primary benchmark, grafting
and PMDI occupancy values remain at their original methods. New smooth-surface
MeCN values are a separate phthalimide control and cannot replace primary-cycle
components. SI Sections 5.1 and 8.1 report the observed sampling coverage and
saved-record reconstruction. The main text explains their implication for
component-specific convergence; it does not describe runtime history.

The occupied-state conformer sets remain incomplete. In the all-same site-2 set,
a late-index minimum lies 3.84 kJ/mol below the best checked early-index state,
while only seven first-ten inputs have checked G. No complete-batch occupied
convergence or global ensemble claim follows. All required candidate criteria,
the public raw-data DOI and final author information remain open; this is a
working draft, not a submission-ready paper.

Run `manuscript/build.py` with the redox interpreter, `PYTHONPATH=src`, and native
CPU thread pools capped at one. It verifies frozen evidence and reconstructs
reported arithmetic before building both PDFs. The new checker also rejects
changed energy/coverage/configuration records, nonfinite frequencies, altered
figure bytes and an unsupported global-convergence flag. Twenty-five focused
tests pass. [Compiler facts](build_check.json),
[numerical checks](claim_check.json) and the
[changed-page visual inspection](visual_review.json) remain separate.

The restored workspace is deliberately small; the original manuscript history,
bulk calculation records and earlier upload bundles remain in the archive
identified by [restoration.json](restoration.json). The current figure snapshot
is kept at `results/figures/manuscript/2026-10-30_polymer_replay/`. Its twenty small
PDFs are versioned as publication inputs, separately from runtime results.
Do not modify existing evidence snapshots; append a new snapshot for verified
scientific additions. Once an Overleaf project is connected, pull/merge before
each edit and again before pushing the checked revision.
