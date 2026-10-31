# Focused redox workflow presentation

Open [redox_workflow.pdf](redox_workflow.pdf): fifteen main slides followed by
four supporting slides. The editable sources are
[redox_workflow.tex](redox_workflow.tex),
[story/motivation_results.tex](story/motivation_results.tex) and
[story/supporting.tex](story/supporting.tex).

The main sequence explains the motivation, criteria, UMA-to-DFT pipeline,
candidate registry, benchmark, exact TEMPO response, base-defined grafting ΔG,
four-state polymer cycle, finite PMDI result and the limits of descriptor ranking.
Slide 4 shows UMA structural search followed by final DFT in MeCN. Slide 8 is the
new exact-viologen reduction ladder; slide 9 gives the TEMPO directional comparison,
and slide 10 explains grafting ΔG. Supporting
slide 18 updates both empty-backbone ensemble comparisons. Each slide has one
claim, its supporting visual and a concise source/scope note.

The 904-graph registry is broader than the 837-design systematic three-core grid
shown as its subset. The paper's older frozen 878-candidate all-axis assessment
is separately scoped. No full material shortlist is claimed. Nitro's measured
force-repeatability failure is explicit; the standard-inference/V100 profile is
separate from unqualified Nitro inference. Molecular descriptor fronts remain
separate from synthesis priorities.

Build with the redox interpreter, `PYTHONPATH=src`, and native CPU thread pools
capped at one: `python slides/build.py`. The build verifies the exact TEMPO
numerical inputs and explicit local figure sources, then records compiler facts
in [build_check.json](build_check.json). Rendered-page inspection is recorded
separately in [visual_review.json](visual_review.json). Figures are collected
under `results/figures/slides/2026-10-30_tempo/`; plotting code stays in
`scripts/plotting/`. The unmodified archived source is identified in
[restoration.json](restoration.json).
