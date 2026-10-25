# Research presentation

All result categories: [results index](../results/README.md).

Open [redox_workflow.pdf](redox_workflow.pdf).

The main story presents the background, post-grafting route, computational
workflow, method, key candidate results, polymer controls and synthesis implications. Detailed
methods, validation tables, supplementary molecular drawings and historical
comparisons follow the **Backup** divider.

The workflow now screens precursor/attachment feasibility before full redox
properties, then compares dimers, trimers and stacking. The main story retains
the BQ candidates as **base required** under conditional base-assisted models. DBU capture differs from the MeCN pKa anchor; DIPEA lacks an exact checked experimental anchor. Attachment selectivity and yield remain untested.

The current 74-page deck has 16 main slides and 58 backup pages. New first-site PMDI grafting results appear on page 19, after the polymer-capacity bounds. The initial grafting-reference checkpoint follows later in Backup. Its dated values and plot are frozen with the submission review; the
[grafting visual report](../results/grafting/index.html) continues updating as
the finite calculation queue finishes.

The main story is generated from the same audited snapshot as the comparison
figures. Numerical findings update with that feed; polymer properties remain distinct from local molecular predictions. Earlier
vinyl-monomer work is preserved as a superseded route exploration.

Refresh the presentation and maintained figures:

```bash
PYTHONPATH=src .venv-figures/bin/python scripts/pipeline/update_verified_slides.py \
  --out results/reviews/2026-10-19/research_presentation/slides --update-plots
```

Edit the main-story generator at
`scripts/plotting/candidates/research_presentation.py`. The slide wrapper is
[redox_workflow.tex](redox_workflow.tex); supporting source is in [backup/](backup/).
The earlier detailed deck is preserved in
[archive/2026-10-19_detailed/](archive/2026-10-19_detailed/README.md).

The current benchmark slide uses the [expanded retrospective panel](../results/validation/publication_benchmark/current.json):
39/39 admitted potentials across 30 exact free molecules and seven studies.
The checked [figure](../results/figures/validation/retrospective_redox_benchmark.pdf)
retains all admitted outliers; the fixed panel is complete. Earlier validation subsets
in Backup are historical snapshots.

The polymer slides include all twelve finite-baseline native G values and their
two-pose partitions. The [24 October spatial evidence](../manuscript/evidence/2026-10-24_spatial_revision/manifest.json)
adds a one-spacer example with both close and separated cores. Its 7.1 kJ/mol
within-graph G difference includes the whole chain and solvent. Broader
convergence, isolated motif binding and electronic coupling remain unqualified.
The quantitative DBU capture-term discrepancy retains its separate scope and
earlier source evidence.

The 25 October neighbor-cycle revision contains 76 pages (17 main and 59 backup). The PMDI matched finite cycle and actual separated-core geometry are included; [the build record](../results/reviews/2026-10-25/pmdi_neighbor_slides/status.json) and [separate rendered review](../results/reviews/2026-10-25/pmdi_neighbor_slides/visual_review.json) bind the checked PDF. Ensemble convergence and material eligibility remain open.
