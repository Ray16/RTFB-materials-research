> **Superseded (2026-11-18).** This per-section draft is kept for history only and is out of date:
> several numbers and claims in it were later corrected. The current text is
> `manuscript/manuscript_draft_v0.md` (journal version: `manuscript/jcim/`; SI:
> `manuscript/supporting_information.md`). Do not cite or copy from this file.

## Validation against experiment

The credibility of any downstream ranking rests on whether the pipeline reproduces
measured electrochemistry, so before applying it to the candidate pool we test it
against literature redox potentials under matched conditions. The validation is
deliberately retrospective and uncorrected: we compare the pipeline's computed
potentials of *free parent* molecules against experimental values measured in
acetonitrile and referenced to ferrocene (Fc), using the identical build → DFT+SMD →
single-point → thermochemistry route applied to every screening candidate. The
experimental panel is hand-curated from the primary literature, with each wave
carrying its source, DOI, table or page, and original reference electrode; every
value is admitted (or not) *before* any computed number is seen, on criteria set by
what the source reports rather than by the residual it produces. Only waves reported
as (quasi)reversible half-wave potentials vs Fc, either directly or through the
authors' own calibration, enter the scored panel (tier A); half-peak potentials,
added-water electrolytes, unverified conversion constants, cathodic-peak-only waves,
and compounds with unspecified stereocenters or inferred identity are tiered down or
excluded with their reason recorded and never scored. Potentials are computed as
E_vs_Fc = -ΔG/n - E_abs(Fc), where the Fc/Fc⁺ reference is evaluated *live* at the
active level of theory with the same free-energy definition (G = E_SMD + G_thermal)
as every molecule, so the reference can never drift from the molecules it calibrates.
No empirical correction, fitted offset, or outlier removal is applied, and the panel
is not held out from model development — it is a descriptive accuracy report, not a
calibrated predictor.

**Headline result.** Across the frozen tier-A panel the pipeline reproduces
experimental potentials with a mean absolute error of **0.247 V**, an RMSE of 0.279 V,
and a **systematic bias of +0.19 V** (computed potentials sit above experiment),
with a maximum absolute error of 0.547 V. The panel comprises **39 scored
endpoints (redox couples) drawn from 30 molecules and seven independent studies**
[source: `manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json`;
frozen archive `manuscript.tex`]. Per-molecule (macro) MAE is 0.263 V. Uncertainty is
quantified by cluster bootstrap (2000 replicates, seed 20261021): resampling by
molecule gives a 95% interval of [0.213, 0.286] V for MAE and [0.128, 0.257] V for
bias; the weaker seven-study clustering gives [0.178, 0.306] V and [0.053, 0.285] V,
reported in the SI. These are descriptive intervals for the observed panel, not
prediction intervals for a new molecule.

We note a discrepancy with an earlier internal brief that quoted "MAE 0.226 V across
44 couples." The frozen benchmark record is unambiguously **39 scored tier-A endpoints (30
molecules, seven studies) with MAE 0.247 V**, obtained by applying the
pre-registered admission rules to the larger raw curated table. The +0.19 V bias
is common to both figures. We therefore report the verified frozen numbers and
treat the earlier "0.226 V / 44 couples" pair as a stale or miscounted citation.

**Interpretation.** The dominant error is a *systematic* positive offset, not scatter:
the bias (+0.19 V) is a large fraction of the MAE, and grouped diagnostics localize it.
The first reduction/oxidation step carries nearly all of it (wave 1 MAE 0.266 V, bias
+0.266 V), while the second step is both smaller and centered (wave 2 MAE 0.185 V, bias
-0.064 V). A systematic offset is exactly the error structure under which *relative*
rankings are far more reliable than absolute values — a uniform shift reorders nothing —
and this is what licenses the comparative screening in the following sections, where
candidates are ranked within pools rather than reported as absolute cell potentials.
We resist attributing the offset to a single cause. The live Fc reference (E_abs(Fc) =
4.379 V) differs from other-functional literature values by an amount comparable to the
first-step bias, making a reference contribution plausible but not proven; reference-free
wave spacings (E1 - E2), which cancel any common Fc shift, are themselves systematically
too large across all nine paired ladders, so a pure reference error cannot explain the
residuals. Charge-state solvation and ion-pairing errors can coexist with it. We state
these as explicit limits: the panel validates *free parents*, not exact grafts or
polymer interactions (`exact_graft_accuracy_validated = False` in the record); external
naphthalene-diimide (NDI) accuracy is not covered, since the available NDI entry was
excluded for inferred identity and unspecified stereocenters; and causal attribution of
the bias to the Fc reference versus solvation remains unresolved (reviewer item M2).

**Figure 4.** The parity plot reports computed vs experimental E_vs_Fc for all 39
endpoints, colored by chemical family (quinone n=22, MAE 0.286 V; aromatic imide n=9,
MAE 0.158 V; phenothiazine n=5, MAE 0.296 V; viologen n=2; nitroxide n=1), with the
1:1 line and the +0.19 V offset line annotated. An accompanying residual histogram
makes the systematic positive shift visible as a displaced, roughly unimodal
distribution, annotated with MAE = 0.247 V and bias = +0.19 V. The worst outlier
(N-ethyl-3,7-dimethoxyphenothiazine, 0.547 V) is retained, not trimmed.

**Reproducibility status (reviewer item M2).** The numbers above are the *frozen
archive* benchmark, regenerated here from the immutable evidence epoch and its
content-hashed panel. The maintained, refactored code paths that produce them —
`properties/potentials.py`, `reporting/benchmark.py`, and the content-hashed panel
builder/validator in `quality/benchmark_panel.py` — are in place and preserve the
scoring, admission, and freeze logic verbatim. Re-confirming that an end-to-end run of
the *maintained* code reproduces the frozen 39-endpoint / MAE 0.247 V / +0.19 V values
requires regenerating `results/redox_potentials.csv` from the DFT energies, which is a
GPU/QM step deferred under the current no-QM constraint (reviewer item M2 pending).
