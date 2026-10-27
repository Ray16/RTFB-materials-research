# Colleague discussion: candidate expansion and validation

Use [the focused presentation](colleague_workflow.pdf): 17 main slides for an
approximately 20-minute discussion, followed by supporting evidence. The maintained
[workflow deck](redox_workflow.pdf) contains the same presentation.

The scientific question is which molecular redox advantages survive attachment
to a preformed polymer. Our workflow connects exact molecular identity, reviewed
attachment chemistry, solution redox calculations, defined-route grafting
thermodynamics and conditional polymer environments. A candidate must meet every
required criterion before entering the material shortlist.

The systematic library contains **837 exact local graft-site designs from
1,593 generation recipes**: 459 viologens, 243 NDI and 135 PMDI designs. The
75-candidate evaluated cohort overlaps **33** library graphs; **32** of those
have both computed potentials. Across the evaluated cohort, 62 candidates have
at least one reported native-calculation potential and 60 have both. All
**26 matched comparisons**, including 18 systematic three-core comparisons,
reconstruct from the unrounded feed. Neither this overlap nor enumeration
constitutes a completed 837-molecule quantum campaign.

The expansion is a **34-slot comparison plan**: seven cores × baseline/shielded ×
flexible/constrained linker, plus six electronic contrasts for BQ, AQ and
phenothiazine. Partial PMDI and phenothiazine matrices now contain mapped
proposals, but their routes are held and the full comparison plan is incomplete.
A separate positive attachment control is proposed. Exact structures, route
evidence and deduplication precede new expensive calculations. Shielding and
rigidification are testable tradeoffs, not assumed improvements.

The completed experimental redox benchmark contains **39 endpoints, 30 exact
free molecules and seven studies in MeCN**. Primary unfitted MAE/RMSE are
**0.247/0.279 V**, with bias **+0.190 V**. The molecule-cluster 95% MAE interval
is **0.213–0.286 V**. First-step errors have mean +0.266 V; second steps have
mean −0.064 V. Post hoc molecule-held-out offset diagnostics give MAE
0.158 V for one offset and 0.106 V for step-specific offsets. These are internal
diagnostics, not prospective validation or replacements for the primary results.
The panel has no NDI endpoint and is quinone-heavy.

Five neutral-partner SVPD reoptimization/Hessian controls change first potentials
by −11.7 to −35.1 mV, holding the existing anion and Fc fixed. The largest
25/50/100 cm⁻¹ rotor-cutoff span over the benchmark is 25.43 mV. These quantify
specific sensitivities; they do not bound missing-conformer error or supply
candidate prediction intervals. Geometry, spin, electronic-energy identity,
native-Hessian checks and independent arithmetic address numerical consistency.
Experimental accuracy and model robustness require separate evidence.

A later [DIPEAH hexafluorophosphate NMR estimate](../results/validation/dipea_experimental_anchor_2026-10-26/README.md)
reports pKa 18.7 in MeCN. A conditional DIPEA/DBU proton-transfer comparison
gives +31.80 kJ/mol computed and +32.02 kJ/mol from the pKa values. The
counterion, activities and temperature conventions are not demonstrated to
match the modeled chloride system. Replacing a single HCl capture term is a
sensitivity rather than a validated correction to the full grafting reaction.
An independent [spectrophotometric measurement](../results/validation/dipea_uv_anchor_2026-10-26/README.md)
reports log Kᴴ = 18.64(1) in MeCN at 25 °C. With the same separate DBU scale,
the implied relative transfer is +32.36 kJ/mol; indicator calibration and
activity conventions still limit the comparison.

The framework now has a finite exact-PMDI example linking redox, precursor
access, grafting and an occupancy cycle. At the primary setting, first-site
DBU/separated-ion attachment is −97.72/−96.96 kJ/mol at two positions. The
neighbor nonadditivity is **−1.86 kJ/mol** in a sampled structure with **16.30 Å**
core separation. This is a whole-chain conditional grafting quantity, not
isolated motif binding or electronic coupling. Selectivity, redox retention,
broader conformer convergence and material-level criteria remain unresolved.

The four sampled phthalimide capped-fragment poses add a distinct local
electronic comparison: close adjacent and one-spacer poses are approximately
−24.8 and −26.2 kJ/mol in the MeCN model, while separated poses are near zero.
These uncorrected frozen-pose values include cap and continuum response; they
are not association free energies, populations or electronic couplings.

Since the frozen 21 October draft, the benchmark progressed from 38/39 to 39/39;
the spacing baseline now has all 12 pose free energies and six finite two-pose
partitions; and the exact PMDI native cycle is complete. Raw-state and Hessian
replay now work after relocation. Public raw-data deposition/DOI and broader
deployment portability remain release tasks. The main manuscript and separate
SI contain the checked results and their scientific limits.

The most useful discussion decisions are the operating window and material
performance criteria, the first route-feasible exact structures for the broader
chemistry panel, and the conformer evidence needed before using spacing trends
for design. The eight molecular descriptor-front models are discussion leads;
**zero candidates currently pass every required material-selection axis**.

Detailed sources and independently reconstructed values are in
[the discussion evidence](../results/reviews/2026-10-26/colleague_workflow/evidence.json).
No new quantum jobs or candidate registrations were required to prepare this summary.
