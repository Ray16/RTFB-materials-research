# Shared redox manuscript

This is the dedicated Overleaf collaboration checkout for the current manuscript.
The initial text and all numerical inputs come from the checked 3 November 2026
export of `/homes/rzhu/0_redox/manuscript`, not the older archive tree.

- Fetch and merge Overleaf immediately before editing; preserve uncommitted work.
- Build both `main.tex` and `supporting_information.tex` with pdfLaTeX and BibTeX.
- After editing, fetch and merge again, rebuild merged changes, and push without
  force. Preserve collaborator changes and report the synchronized commit.
- `main.tex` is the Overleaf main document. The SI builds separately by selecting
  `supporting_information.tex` as the main document.
- Keep required figure PDFs in Git; ignore only generated document PDFs.
- `inputs/` and computational figures are evidence-bound outputs. Change numerical
  results through the scientific pipeline and refresh the corresponding evidence,
  rather than editing computed values in this collaboration checkout.
- Once connected, do not overwrite this checkout with a fresh export. Reconcile
  prose changes with the scientific repository while preserving its path mappings.
- This draft is not submission-ready. Author metadata, public deposition, a
  permanent code release and the scientific review issues remain outstanding.
