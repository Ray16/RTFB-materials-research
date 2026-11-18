# Figure 4 — validation parity

**Figure 4.** Computed versus experimental redox potentials (V vs Fc/Fc+, in MeCN) for
the tier-A benchmark panel. Each point is one redox endpoint (a single measured wave /
couple), coloured by redox family (quinone, aromatic imide, phenothiazine, viologen,
nitroxide). The solid line is the 1:1 identity; the shaded envelope is the +/-MAE
agreement band. Over **n = 39 endpoints from 30 molecules and 7 independent studies**,
the protocol reaches **MAE = 0.247 V, RMSE = 0.279 V**, with a **systematic bias of
+0.19 V** (computed potentials run high); the largest single error is 0.55 V (the
second oxidation of N-ethyl-3,7-dimethoxyphenothiazine), kept in the panel. By family:
quinone n = 22, MAE 0.286 V; aromatic imide n = 9, MAE 0.158 V; phenothiazine n = 5,
MAE 0.296 V; viologen n = 2; nitroxide n = 1. Level of theory: geometries at
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
removal, or acceptance threshold was fit to these errors. The values were computed with
the code version at commit `46b39c3b` (tag `legacy/pre-refactor-main`); the plotted
potentials use the DFT-frequency thermal tier, and the same panel at the pool's xTB
tier gives MAE 0.232 V (SI S3).

Source: `archive/0_redox/manuscript/evidence/2026-10-22_benchmark_complete/benchmark.json`.
Generator: `scripts/plotting/manuscript/fig4_validation_parity.py`.
