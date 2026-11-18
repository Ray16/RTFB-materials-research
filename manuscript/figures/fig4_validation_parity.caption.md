# Figure 4 — validation parity

**Figure 4.** Computed versus experimental redox potentials (V vs Fc/Fc+, in MeCN) for
the tier-A benchmark panel. Each point is one redox endpoint (a single measured wave /
couple), coloured by redox family (quinone, aromatic imide, phenothiazine, viologen,
nitroxide). The solid line is the 1:1 identity; the shaded envelope is the +/-MAE
agreement band. Over **n = 39 endpoints from 30 molecules and 7 independent studies**,
the protocol reaches **MAE = 0.247 V, RMSE = 0.279 V**, with a **systematic bias of
+0.19 V** (computed potentials run slightly high); the largest single error is 0.55 V
(an N,N'-disubstituted phenothiazine second oxidation). Level of theory: geometries at
r2SCAN-D4/def2-SVP(D)/SMD(MeCN), energies at wB97M-V/def2-TZVPD/SMD(MeCN).

Because the error is dominated by a near-constant positive offset rather than scatter
(wave-1 bias +0.27 V, wave-2 bias -0.06 V), **relative rankings and trends within a
family are more reliable than absolute potentials**; absolute values should be read with
the +0.19 V offset in mind. The companion panel **Figure 4b** shows the signed-residual
distribution, making the offset and its spread explicit.

**Scope and limits.** The panel validates **free soluble parent molecules, not capped
grafts or polymer-bound/oligomeric forms**, which are distinct models and are not
established by this comparison. Naphthalene diimide (NDI) endpoints are excluded from the
scored set. The confidence intervals are panel-resampling (bootstrap) intervals over the
observed endpoints, not predictive intervals, and no empirical correction, outlier
removal, or acceptance threshold was fit to these errors. Values are from the **frozen
2026-10-22 benchmark archive**; a re-run confirming the maintained/refactored code
reproduces them (the "M2" QM re-run) is still pending under the no-GPU constraint, so
these are not yet "validated by the maintained code."

Source: `archive/0_redox/manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json`.
Generator: `scripts/plotting/manuscript/fig4_validation_parity.py`.
