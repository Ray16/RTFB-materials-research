# Redox polymer manuscript

Shared project: https://www.overleaf.com/project/6aea2dc8096db15da57980a0

Main document: `main.tex`. Supporting Information: `supporting_information.tex`.
The current narrative revision compiles to a 15-page article and 41-page SI.
The subsequent clarity edit is recorded in `language_revision.json`;
`language_author_notes.md` records their resolution in `method_clarification.json`.
See `narrative_revision.json`, `narrative_claim_check.json`,
`narrative_build_check.json` and `narrative_visual_review.json` for the changes
and checks. The initial export/build records below describe the preserved starting version.
Compiler: pdfLaTeX. Select the SI as the main document when building it separately.

This collaboration checkout was initialized from the checked 3 November 2026
article (21 pages) and SI (36 pages) in `/homes/rzhu/0_redox/manuscript`.
Figures, bibliography and frozen numerical TeX inputs are included. The article
is a working draft; this upload does not establish submission readiness.

Fetch and merge Overleaf immediately before editing. After each edit pass,
build both documents, fetch and merge again, build any merged changes, then
push to `overleaf/main`. Preserve collaborator edits and never force-push.
Use `/homes/rzhu/0_redox/manuscript/overleaf` for collaborative prose edits. Scientific
code and evidence remain in `/homes/rzhu/0_redox`; computed values and figures
must be updated through that pipeline and reconciled with the shared text.
Do not replace the collaboration checkout with a fresh export.

Git authentication uses the username `git` and an Overleaf Git token entered
in the terminal. Never commit tokens. This checkout uses an eight-hour
in-memory credential cache, so authentication may be needed on another host
or after the cache expires.
