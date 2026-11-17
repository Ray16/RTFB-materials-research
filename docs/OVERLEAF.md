# Shared redox manuscript on Overleaf

Project: https://www.overleaf.com/project/6aea2dc8096db15da57980a0

The active collaboration checkout is `/homes/rzhu/0_redox/manuscript/overleaf`
(canonical path `/nfs/lambda_stor_01/homes/rzhu/0_redox/manuscript/overleaf`).
Its Git remote is `https://git@git.overleaf.com/6aea2dc8096db15da57980a0`,
with local and remote branch `main`.

The initial upload uses the checked 3 November 2026 manuscript in
`/homes/rzhu/0_redox/manuscript`, rather than the older
`/homes/rzhu/archive/0_redox/manuscript` snapshot. The Overleaf starter history
is retained in the merge. The article is `main.tex`; select it as Overleaf's
main document and use pdfLaTeX. Select `supporting_information.tex` to compile
the separate SI. Both documents compile locally and their extracted page text
matches the checked source PDFs: 21 article pages and 36 SI pages. Computational
inputs and figures are included; full raw calculations remain on shared storage.

The initial uploaded commit is `9aee5adb2a2c83a86637642ecb3478aaa606344e`.
The machine-readable connection and build record is [overleaf_connection.json](overleaf_connection.json).

## Current reviewed revision

Overleaf commit `5f50c614060e1470d32a104cbf2a1ed9624b5aea` clarifies the calculation levels and Fc scale comparison.
The earlier language edit remains preserved at `b881af03810dc49baa83338f90f2c4356d93a8d1`.
The first systematic narrative revision remains preserved at `2dea0497b637de6bdd2012522e4148e882629b73`.
The article is 15 pages and the SI is 41 pages. The
[manuscript entry point](../manuscript/README.md) links the checked PDFs,
self-contained source snapshot and claim/build/visual records.
No new quantum calculations were performed; the original numerical inputs
and fixed 39-endpoint benchmark are preserved.

## Editing and synchronization

Fetch and merge Overleaf immediately before every manuscript edit pass, preserving
existing local work. Work in the collaboration checkout. After editing, build
and check both documents, fetch and merge again, rebuild merged changes, then
push to `overleaf/main`. Never force-push or replace collaborators' edits with
a regenerated export. Report the pushed commit after each completed pass.

The scientific pipeline and frozen evidence remain in `/homes/rzhu/0_redox`.
Map `main.tex` back to `manuscript/manuscript.tex` when reconciling prose;
exported `inputs/` and `figures/` paths differ from the scientific source paths.
Update calculated values through the evidence pipeline and merge those updates
with the collaborative text. The earlier export and source snapshots remain
historical records; they must not overwrite this shared document.

Git authentication uses username `git` and an Overleaf Git token entered at a
terminal password prompt. The collaboration checkout uses an eight-hour
in-memory credential cache; another host or an expired cache may require login.
Keep credentials out of versioned files.

This synchronization establishes collaboration access, not submission readiness.
Scientific revision, final author metadata, archival deposition and a permanent
code/data release remain outstanding.
