# Focused workflow presentation

Open [colleague_workflow.pdf](colleague_workflow.pdf), also maintained as
[redox_workflow.pdf](redox_workflow.pdf). The deck has **16 main slides and
7 supporting slides**, intended for an approximately 20-minute discussion.

The narrative connects systematic candidate expansion, matched design effects,
validation of computed properties, defined-route grafting and the finite PMDI
neighbor example. Each slide has one purpose and a take-home claim. Detailed
QC and implementation evidence appear only in backup.

- [Speaker notes](colleague_speaker_notes.md) explain each slide and likely questions.
- [Scientific brief](colleague_brief.md) records the current findings and limits.
- [Editable source](redox_workflow.tex) uses the [main story](generated/colleague_story.tex)
  and [supporting evidence](generated/colleague_backup.tex).
- [Presentation evidence](../results/reviews/2026-10-26/colleague_workflow/README.md)
  contains frozen inputs, independent numerical checks, build and separate visual review.
- The [DIPEA NMR source revision](../results/reviews/2026-10-26/colleague_workflow_dipea_nmr/build_check.json)
  and its [visual check](../results/reviews/2026-10-26/colleague_workflow_dipea_nmr/visual_review.json)
  preserve the earlier salt-only presentation. The [current build](../results/reviews/2026-10-26/colleague_workflow_dipea_uv/build_check.json)
  and [visual check](../results/reviews/2026-10-26/colleague_workflow_dipea_uv/visual_review.json)
  bind slide 12 to both checked experimental sources.
- [Plotting code](../scripts/plotting/presentation/plot_colleague_workflow.py)
  produces the maintained [presentation figures](../results/figures/presentation/).

The previous 76-page deck remains in
[the earlier-deck archive](archive/2026-10-26_before_focused_workflow/README.md).
Historical scientific records and manuscript evidence are unchanged.

The expansion denominators are deliberately separate: 837 enumerated designs,
75 evaluated-cohort candidates and 34 future comparison slots. The evaluated
cohort overlaps 33 enumerated graphs, of which 32 have both computed potentials.
The 34-slot plan has generated zero new exact molecules. The benchmark has
39/39 exact-parent MeCN endpoints; it does not validate exact-graft or polymer
accuracy. No candidate currently passes every required material-selection axis.

For reproduction, use one CPU and set OMP_NUM_THREADS, OMP_THREAD_LIMIT,
MKL_NUM_THREADS, OPENBLAS_NUM_THREADS, NUMEXPR_NUM_THREADS and XTB_THREADS to 1.
Set PYTHONPATH=src. Run the following scripts with their --help options:

1. `scripts/presentations/collect_colleague_evidence.py` with the figure environment,
   writing a new evidence directory (existing snapshots are never overwritten).
2. `scripts/plotting/presentation/plot_colleague_workflow.py` with the same figure
   environment, the checked evidence and `results/figures/presentation` output.
3. `scripts/presentations/build_colleague_workflow.py` with pypdf and pdflatex,
   a new build directory, the checked evidence and `--dipea-evidence` pointing
   to the checked DIPEA NMR summary, plus `--dipea-uv-evidence` pointing to
   the checked DIPEA UV summary.
4. `scripts/presentations/render_colleague_workflow.py` with Pillow and Ghostscript,
   followed by separate visual inspection.

The figure environment is `.venv-figures/bin/python`; the project redox interpreter
supplies pypdf and Pillow. Respect shared-node CPU allocation. No script launches
quantum calculations. Compiler checks and visual-review records remain separate.
