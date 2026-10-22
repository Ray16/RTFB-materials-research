# Working in silico manuscript

Open [manuscript.pdf](manuscript.pdf). Edit [manuscript.tex](manuscript.tex);
the [outline](outline.md) records the intended argument and remaining analyses.
The [early-review notes](review_notes.md) identify questions for a colleague.
This is a local working draft, with no Overleaf project configured.

The draft contains Introduction, Methods, verified Results, Discussion and
Conclusions, with citations and six figures including a working graphical
abstract. It uses a frozen [evidence snapshot](evidence/2026-10-22/manifest.json)
so live calculation updates cannot silently change a figure while leaving its
caption and text unchanged. Figures are under `results/figures/manuscript/`;
the original maintained figures and raw data remain at their existing locations.

The current paper is methods-focused. It reports 38/39 admitted free-molecule
potentials, defined DBU/DIPEA grafting comparisons and conditional finite-polymer
models. All six selected spacing free energies and both conditional cycles are
checked. Pending trimer, spin, method-sensitivity, external-validation and
release work is marked in the manuscript. These labels must be resolved before
submission. Authorship, affiliations, funding and competing-interest statements
are awaiting author confirmation.

Build using the existing figure environment, with native thread pools capped at one:

```bash
PYTHONPATH=src .venv-figures/bin/python manuscript/build.py
```

The build requires `pdflatex`, `bibtex`, `pypdf` and the project resource helper.
It checks the frozen evidence, independently reconstructs the reported benchmark
and reaction-cycle arithmetic, compiles citations and cross-references, and
writes [build_check.json](build_check.json). Regeneration after new scientific
results requires a new evidence snapshot and an intentional text/figure update.

## Overleaf upload

The [self-contained ZIP](overleaf_upload.zip) contains `main.tex`, the bibliography
and all five external figure PDFs. The working graphical abstract is drawn by
LaTeX. The exported project compiles to the same 13-page text as the checked
draft; every archive member is hash-verified in [overleaf_export.json](overleaf_export.json).
No Overleaf project has been created or connected yet.

In Overleaf, choose **New Project → Upload project** and upload this ZIP. Set
`main.tex` as the main document if it is not detected automatically.
Once a Git connection is configured, fetch and merge collaborators’ changes
before editing and again before pushing. The isolated local export lives outside
the computational repository so its Overleaf Git history can remain separate.
