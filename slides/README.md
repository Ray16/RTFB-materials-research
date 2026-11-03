# Focused redox workflow presentation

Open [redox_workflow.pdf](redox_workflow.pdf): sixteen main slides followed by
seven supporting slides. The editable sources are
[redox_workflow.tex](redox_workflow.tex),
[story/motivation_results.tex](story/motivation_results.tex) and
[story/supporting.tex](story/supporting.tex).

The main sequence explains the motivation, criteria, UMA-to-DFT pipeline,
candidate registry, benchmark, exact TEMPO response, base-defined grafting ΔG,
four-state polymer cycle, finite PMDI result and the limits of descriptor ranking.
Slide 4 shows UMA structural search followed by final DFT in MeCN. Slide 8 is the
exact-viologen reduction ladder and −11.49 kJ/mol local quaternization; slide 9 gives the TEMPO directional comparison,
and slide 10 explains grafting ΔG. Slide 13 compares whole-chain and frozen-motif electronic definitions. Supporting
slide 19 reports the first complete occupied phthalimide-control sampling comparison. Each slide has one
claim, its supporting visual and a concise source/scope note.

The 912-graph registry is broader than the 837-design systematic three-core grid
shown as its subset. The paper's older frozen 878-candidate all-axis assessment
is separately scoped. No full material shortlist is claimed. Nitro's measured
force-repeatability failure is explicit; the standard-inference/V100 profile is
separate from unqualified Nitro inference. Molecular descriptor fronts remain
separate from synthesis priorities.

Build with the redox interpreter, `PYTHONPATH=src`, and native CPU thread pools
capped at one: `python slides/build.py`. The build verifies the exact TEMPO
and exact-viologen redox/grafting numerical inputs and explicit local figure sources, then records compiler facts
in [build_check.json](build_check.json). Rendered-page inspection is recorded
separately in [visual_review.json](visual_review.json). Figures are collected
under `results/figures/slides/2026-10-30_tempo/`; plotting code stays in
`scripts/plotting/`. The unmodified archived source is identified in
[restoration.json](restoration.json).

Supporting slide 20 reports all four checked finite ionic neighboring-graft G
cycles: −3.94 kJ/mol for all-same adjacent, +22.58 for all-same one-spacer,
+8.28 for alternating adjacent, and −5.30 for alternating one-spacer at the
primary cutoff. Twelve independently checked local minima support the cycles.
Both signs reverse with backbone configuration. Occupied ensembles remain pending. Slide 21 reports six primary chloride-contact outcomes in five duplicate groups
(1-M association G +5.33 to +22.36 kJ/mol), with the spin and particle-count
conventions explicit. Slide 22 reports the primary molecular control; slide 23
retains reproducibility limits. The current immutable manuscript
evidence epoch also supplies twelve motif triplets and 72 independently checked
SCFs for slide 13. Its electronic comparison is distinct from binding G,
electronic coupling and an experimental grafting distribution.

Supporting slide 19 now shows the failed predeclared 10-to-20 occupied sampling
check: -2.40 kJ/mol ensemble G change at the main rotor cutoff and -3.88 kJ/mol
minimum electronic-energy change, beyond the 2-kJ/mol tolerance. All three
rotor settings and numerical duplicate thresholds retain the failure. Its
[frozen evidence](evidence/2026-11-02_occupied_sampling/manifest.json),
[geometry/partition replay](check_sampling.py) and standalone
[figure](../results/figures/slides/occupied_sampling_convergence.pdf) are bound
into the build checks. This is one complete corrected phthalimide control cell;
the other eleven cells, balanced occupied ensembles and material admission
remain pending. The earlier empty-backbone comparisons remain in the archived
source/evidence history. All 23 pages were rendered against the prior checked
deck; only page 19 changed and passed visual inspection.

The subsequent phenothiazine update adds its checked local-model first oxidation
(+0.626 V versus Fc in MeCN) to slide 5. The deck remains 23 pages. Its value
uses the shared frozen-evidence checker; neither route yield nor polymer
performance is inferred. The current visual-review record supersedes the
earlier changed-page observations above.
