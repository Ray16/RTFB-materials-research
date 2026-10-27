# Focused workflow presentation

Open [colleague_workflow.pdf](colleague_workflow.pdf), also maintained as
[redox_workflow.pdf](redox_workflow.pdf). The deck has **17 main slides and
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
  contains the original frozen discussion inputs and independent numerical checks. The current build and visual review are linked below.
- The [DIPEA NMR source revision](../results/reviews/2026-10-26/colleague_workflow_dipea_nmr/build_check.json)
  and its [visual check](../results/reviews/2026-10-26/colleague_workflow_dipea_nmr/visual_review.json)
  preserve the earlier salt-only presentation. The [DIPEA revision build](../results/reviews/2026-10-26/colleague_workflow_dipea_uv/build_check.json)
  and [visual check](../results/reviews/2026-10-26/colleague_workflow_dipea_uv/visual_review.json)
  bind slide 12 to both checked experimental sources. The [AQ-proposal build](../results/reviews/2026-10-26/colleague_workflow_aq_proposal/build_check.json)
  and [changed-slide review](../results/reviews/2026-10-26/colleague_workflow_aq_proposal/visual_review.json)
  add the separate held AQ proposal to slide 6. The [latest build](../results/reviews/2026-10-26/colleague_workflow_empty_range/build_check.json)
  and [changed-slide review](../results/reviews/2026-10-26/colleague_workflow_empty_range/visual_review.json)
  add the alternative-model empty-chain range to slide 15. The [latest build](../results/reviews/2026-10-26/colleague_workflow_empty_set/build_check.json)
  and [changed-slide review](../results/reviews/2026-10-26/colleague_workflow_empty_set/visual_review.json)
  clarify that the four geometries are distinct under the checked comparison.
  The [latest build](../results/reviews/2026-10-26/colleague_workflow_ndi_challenge/build_check.json)
  and [changed-slide review](../results/reviews/2026-10-26/colleague_workflow_ndi_challenge/visual_review.json)
  note the sourced, uncomputed NDI challenge on slide 8. The earlier [TEMPO build](../results/reviews/2026-10-26/colleague_workflow_tempo_exact_support_v2/build_check.json)
  established the source-supported route example. The [current build](../results/reviews/2026-10-26/colleague_workflow_local_contact_v2/build_check.json)
  and [visual review](../results/reviews/2026-10-26/colleague_workflow_local_contact_v2/visual_review.json) add a focused local-contact result as slide 15; TEMPO is now slide 16.
- [Plotting code](../scripts/plotting/presentation/plot_colleague_workflow.py)
  produces the maintained [presentation figures](../results/figures/presentation/).

The previous 76-page deck remains in
[the earlier-deck archive](archive/2026-10-26_before_focused_workflow/README.md).
Historical scientific records and manuscript evidence are unchanged.

The expansion denominators are deliberately separate: 837 enumerated designs,
75 evaluated-cohort candidates and 34 future comparison slots. The evaluated
cohort overlaps 33 enumerated graphs, of which 32 have both computed potentials.
Partial PMDI and phenothiazine matrices now contain mapped proposals, but the 34-slot plan is not complete and their new routes are held. The benchmark has
39/39 exact-parent MeCN endpoints; it does not validate exact-graft or polymer
accuracy. No candidate currently passes every required material-selection axis.

Separately from that 34-slot plan, one [AQ-2-carboxylate graft proposal](../results/discovery/aq_carboxylate_2026-10-26/README.md)
has been mapped and route-screened. Its route is held and its redox properties
are not computed; the current deck's 837 count refers to the original library.

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
   the checked DIPEA UV summary and `--tempo-evidence` pointing to the TEMPO
   route and potential summary.
4. `scripts/presentations/render_colleague_workflow.py` with Pillow and Ghostscript,
   followed by separate visual inspection.

The figure environment is `.venv-figures/bin/python`; the project redox interpreter
supplies pypdf and Pillow. Respect shared-node CPU allocation. No script launches
quantum calculations. Compiler checks and visual-review records remain separate.
