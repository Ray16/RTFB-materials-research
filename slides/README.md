# Redox workflow presentation

Open [colleague_workflow.pdf](colleague_workflow.pdf), also maintained as
[redox_workflow.pdf](redox_workflow.pdf). The rebuilt deck has **13 main slides
and 4 supporting slides** for an approximately 15–20 minute scientific discussion.

- **Motivation:** slides 1–3 explain the material-design problem.
- **UMA pipeline:** slide 4 shows exact models → route review → RDKit/UMA
  gas-phase warm starts → DFT in acetonitrile → the three property analyses.
  A [standalone pipeline PDF](../results/figures/presentation/story/uma_pipeline.pdf)
  is available for reuse.
- **Expansion and tuning:** slides 5–6 distinguish the 837-design library from
  the evaluated cohort and show matched PMDI voltage shifts.
- **Validation:** slide 7 reports the 39-endpoint experimental benchmark.
- **Grafting ΔG:** slide 8 compares the HCl reference with the DBU-assisted route.
- **Polymer models and results:** slides 9–11 explain the occupancy cycle,
  finite PMDI neighbor effect and frozen local-motif interactions.
- **Selection and lessons:** slides 12–13 distinguish descriptor leaders from
  candidates supported on every required material axis.

The editable source is [redox_workflow.tex](redox_workflow.tex), with the main
story and supporting slides in [story/](story/). [Speaker notes](colleague_speaker_notes.md)
record each slide's purpose, result and interpretation. The earlier deck is
preserved in [archive/2026-10-27_before_story_rebuild/](archive/2026-10-27_before_story_rebuild/).

Use `scripts/plotting/presentation/plot_redox_story.py` to regenerate the three
new charts from frozen evidence; other panels retain their checked existing
figure sources. Use `scripts/presentations/build_redox_story.py --work NEW_DIR
--review NEW_DIR` with the redox interpreter, then render and visually review
with `scripts/presentations/render_colleague_workflow.py`. Cap native CPU pools
at one thread. The builder checks input hashes, compiler errors, slide count,
and explicit UMA/grafting content; visual review is recorded separately.

[Current build check](../results/reviews/2026-10-27/uma_story_v2/build_check.json)
and [visual review](../results/reviews/2026-10-27/uma_story_v2/visual_review.json)
bind the current PDF. No quantum calculation was launched for this revision.
