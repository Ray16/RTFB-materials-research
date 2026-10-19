# Phenothiazine validation: second-wave model limitation

Three exact parent structures supply four currently checked measured potentials.
EPT first wave: 0.537 vs 0.278 V; MEEPT first: 0.559 vs 0.311 V;
DMeOEPT first/second: 0.249/1.202 vs 0.056/0.655 V (all vs Fc).

The **+0.547 V** DMeOEPT second-wave error is a material limitation for two-electron
phenothiazine window/capacity predictions. The calculated charging gap is
0.95309 eV versus 0.599 eV measured. Electronic and thermal contributions are
0.91365 and 0.03944 eV, so a reference offset or thermal-only adjustment cannot
remove the discrepancy. No value, threshold or reference is changed to fit it.

First-wave matched shifts are better: MEEPT minus EPT is +0.022 V calculated vs
+0.033 V measured; DMeOEPT minus EPT is -0.288 vs -0.222 V. These two contrasts
have errors -0.011 and -0.066 V; they do not establish general family accuracy.

The frozen-geometry RKS reproduction / zero-net-spin unrestricted-seed diagnostic
is under `results/validation/ptz_dication_scf_probe_2026-10-19/`. It tests numerical
and orbital-solution possibilities before assigning the residual to model or
electrolyte effects. Its outputs cannot replace canonical records. Calculated
values remain available as model predictions; confident two-electron window or
capacity claims for grafted analogues remain provisional.
