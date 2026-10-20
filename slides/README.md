# Research presentation

All result categories: [results index](../results/README.md).

Open [redox_workflow.pdf](redox_workflow.pdf).

The first 13 slides present the background, post-grafting route, computational
workflow, method, key candidate results and synthesis implications. Detailed
methods, validation tables, supplementary molecular drawings and historical
comparisons follow the **Backup** divider.

The initial grafting-reference checkpoint is in Backup, immediately after the
divider. Its dated values and plot are frozen with the submission review; the
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
