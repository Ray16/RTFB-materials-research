# Exact TEMPO directional solution relaxation

Six fresh electronic calculations for the same capped TEMPO graft now pass
independent energy, density, spin, geometry, source and resource checks. The
result advances one exact candidate without completing its material assessment.

| Quantity | Reduction (+1 to 0) | Oxidation (0 to +1) |
| --- | --- | --- |
| Equilibrium geometry component / eV | 0.396266 | 0.552630 |
| Frozen-slow solvent increment / eV | 1.080856 | 1.086011 |
| Directional solution relaxation / eV | 1.477122 | 1.638641 |

The exact product is `Cc1ccc(COC2CC(C)(C)N([O])C(C)(C)C2)cc1`, InChIKey
`KTIGVFIKAXWERK-UHFFFAOYSA-N`. Its separately reconstructed oxidation potential
is +0.578775 V versus Fc/Fc+ in MeCN. The reviewed attachment precedent uses
4-hydroxy-TEMPO, preformed Merrifield resin and NaH in DMF or THF; reaction
conditions and operating MeCN are distinct contexts.

The calculations use wb97m-v/def2-TZVPD electronic energies at two checked
r2SCAN-D4 geometry-level solution minima. Static/optical dielectric constants are
35.688/1.80687364. Fast polarization relaxes while the initial slow solvent
response is frozen; the geometry-only SMD CDS term is included once. No nuclear
thermal, gas-inner or Born-outer term is added. The endpoints differ by 0.137 Å
heavy-atom RMSD; three bridge torsions change by at most 10.85 degrees. These are
paired-basin diagnostics, not a converged conformer ensemble or proof of active
single-point stationarity.

The two directions differ by 10.4%, exceeding the predeclared 5%
equal-curvature diagnostic. Both values remain visible. The qualified symmetric
scalar, grafting free energy, occupancy nonadditivity and electronic coupling
remain blank in [candidate.json](candidate.json). Its all-axis decision is
pending. These electronic response costs do not establish a measured rate,
empirical prediction accuracy or a feasible material.

[independent_review.json](independent_review.json) binds the six raw SCF/density
records and original native geometry/redox sources. The archived evaluator's
103 independent checks pass, supplemented by explicit spin-density and
full-graph checks in the maintained reviewer. Both own-geometry energies reproduce
their native source active single points within 2e-9 eV. The old scientific
status file is retained unchanged: it says review pending because the independent
review is a separate artifact. Use the latter for the completed QC disposition.

The V100 run used the shared GPU gate, two disjoint CPU slots and six capped
thread pools. A prospective 16-GiB host guard passed this case, with tree-RSS peak
1.49 GiB. GPU peak was 16,518 MiB: this result does **not** qualify 16-GB GPUs.
The normal full polymer Hessian resource class remains separate.

Durable raw arrays, six SCF/density pairs, settings and as-used sources are at
`/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-10-30_TEMPO_solution_pair/`.
The original archive is mounted read-only at its recorded prefix; new outputs
are written separately. No archived quantum record is rewritten.

The maintained construction, guarded worker, independent observer, exact-case
publisher and plot are versioned in `src/redox/workflows/solution_pair_namespace.py`,
`scripts/pipeline/run_checked_solution_pair.py`,
`scripts/pipeline/observe_checked_solution_pair.py`,
`scripts/pipeline/publish_checked_solution_pair.py` and
`scripts/plotting/candidates/plot_solution_pair.py`. The small numerical inputs
are independently reconstructed by `manuscript/check_tempo_solution.py`.
Fifty-eight focused tests include rejection of changed energies, reference
densities, roles, charge/spin, geometry, figures, generated numerical macros and
unsupported promotion. The final paper, SI and focused slides include the checked
result, with visual inspection recorded separately from compiler facts.

The paper's frozen 878-candidate all-axis assessment and the maintained
904-graph registry are different scopes; the presentation labels the latter as
a registry, not a fully evaluated shortlist. Expanding exact-route evidence and
completing candidate/polymer criteria remain open under M1 and M4. Public release
assembly and deposition remain open under M7.
