# Current disposition of colleague review

As of 31 October 2026, the manuscript remains a working draft. M1, M4, complete M7 and final M9 remain open. The table distinguishes scoped reviewer responses from complete submission readiness.

**Additional release requirement:** the [screening implementation assessment](../validation/screening_implementation_2026-10-31/README.md)
reproduces an upstream GPU4PySCF grid-screening defect on V100 and H100.
The isolated correction passes fixed-density and closed/open-shell MeCN
derivative controls on both devices. All 71 benchmark electronic states now
pass independent repeat review; the largest corrected-minus-installed
electronic contribution to a relative-potential change is 8.66e-11 V.
Full polymer-Hessian, reaction-cycle and production qualification remain open. Previous arithmetic,
experimental-source and saved-record replays remain reproducible; they do not
measure the defect's effect on converged chemical predictions. Stored chemical
claims are provisional for release. New uncorrected admissions are held.

| Comment | Status | Current disposition |
| --- | --- | --- |
| M1 | partial | Exact PMDI finite cycle, baseline enumerated viologen native ladder plus three-seed local quaternization and TEMPO directional MeCN descriptors reported; converged ensembles and complete material/route criteria remain open. |
| M2 | addressed in current scope | 39 endpoints and diagnostic controls reported; reference-error attribution and external NDI accuracy remain explicitly limited. |
| M3 | addressed in current scope | DBU/HCl and two DIPEA source estimates checked; full-route signs remain conditional. |
| M4 | open | Both twenty-input empty components meet the batch criterion with 20/19 distinct structures; occupied ensembles and balanced cycles remain incomplete. |
| M5 | addressed in current scope | Five controls reported; no complete basis-convergence claim. |
| M6 | addressed in current scope | Methods and twelve-pose SI table report the gates; the tested Nitro repeatability failure is explicit and production qualification is withheld. |
| M7 | partial | The current 75-state polymer cache also passes extracted-archive replay in a fresh NumPy-only environment; a ten-state exact-viologen CPU example is now versioned and passes fresh NumPy-only extracted replay; public DOI, final metadata, complete release assembly and wider deployment portability remain open. |
| M8 | addressed in current scope | The 48-reference base audit is preserved; Achazi 2024, CREST 2020, the Oh Biying 2025/Mohanty 2026 polymer studies and You 2015 nonequilibrium-solvation theory are separately verified additions. Every benchmark study is cited, with source/reference handling tabulated. |
| M9 | revised final science pending | Current article and separate SI are revised; final narrative depends on open physical results. |
| remaining_endpoint | addressed in current scope | 39 endpoints reconstruct; alternative-spin thermochemistry is still stated as a limit. |
| m1 | addressed in current scope | Current immutable epoch and separately built PDFs replace stale reporting. |
| m2 | addressed in current scope | Molecule/study intervals and their weak calibration are explicit. |
| m3 | addressed in current scope | Actual conformers and exact occupancy construction appear in main/SI figures. |
| m4 | addressed in current scope | Main text states the exact identity relationship. |
| m5 | addressed in current scope | SI reports signed fragment terms, cancellations, 57–93 steps and 0.25–0.70 Å relaxation. |
| m6 | addressed in maintained code | Central conversion now uses exact SI defining constants; 68 focused checks retain archived arithmetic without rewriting raw scientific records. |
| m7 | addressed in current scope | Build emits compiler facts and requests a separate PDF-hash-bound visual record. |
| m8 | addressed in current scope | Current chemical artwork and exact models replace status banners; no extra author-approval gate is inferred. |
| m9 | addressed in current scope | Remaining evidence/decision terms are defined in an SI table. |
| m10 | addressed in current scope | Initial seed/counts, one-configuration scope and xTB symmetry detection are explicit. |

The current article and SI are 21/29 pages; the focused deck has fifteen main and four supporting slides. Both PDFs compile without overfull boxes and the revised pages are visually checked. The self-contained upload compiles with matching text; no Overleaf project is connected.

The exact baseline enumerated viologen has checked MeCN reductions of −0.608/−1.208 V versus matched Fc and a finite local quaternization ΔG of −11.49 kJ/mol at 298.15 K. Three precursor starts pass independent geometry, derivative, nuclear, spin and resource checks. The seed span is 1.50 kJ/mol; rotor settings give −13.30 to −11.49 kJ/mol. This is a separated-ion capped-site model, with no HCl/base term. Salt association, exchange, reflux thermodynamics and polymer penalties remain unresolved. [Sources and reproduction](../validation/viologen_grafting_2026-10-31/README.md).

The TEMPO example retains separate directional MeCN relaxation costs of 1.477/1.639 eV. Their 10.4% asymmetry fails the declared equal-curvature diagnostic, so a symmetric scalar remains unqualified. Neither example is a fully admitted material.

The frozen polymer comparison contains 75 reconstructed states. Both twenty-input empty-backbone sets meet the finite-batch criterion with 20/19 distinct structures. The [working census](../validation/screening_implementation_2026-10-31/native_census.json) now contains 122 independently reconstructed input states outside that frozen comparison; balanced occupied ensembles remain incomplete. The original atom-31/z sampled-derivative failure remains preserved at 1.2138%, above the unchanged 1% gate. A separate full native repeat passes the existing protocol, with a maximum column error of 0.7044%, and passes independent physical reconstruction. Force repeatability and finite-difference step convergence remain distinct open questions; neither the repeat nor its census count qualifies the corrected screening implementation.

The 75-state extracted cache reconstructs Hessians, sampled derivative checks and nuclear terms in a fresh NumPy-only environment with original paths blocked. This does not establish complete historical executable closure, wider deployment portability or public deposition. Final archive assembly, DOI and author/deposition metadata remain open.

The older [machine-readable comment audit](../results/validation/review_comment_audit_2026-10-26/current.json) binds its historical inspected sources; it is not a certificate for the current revision. [Publication acceptance criteria](../docs/PUBLICATION_READINESS.md) remain authoritative for completion.

The exact viologen now has [twelve audited initial charged-polymer models](../validation/viologen_polymer_preparation_2026-10-31/README.md), including adjacent/spacer graphs in two configurations. SI Section 3.4 separates connectivity, global state hypotheses and PF6 dry-salt mass bookkeeping from uncomputed charged-polymer thermochemistry. This preparation advances M1 without closing its polymer assessment or M4.

An interrupted site-4 case now passes independent reconstruction. V100 host-allocation controls reproduce a NumPy huge-page-advice stall. The earlier seven-block/315-coordinate response execution probe used solute H1 and a frozen-solvent operator; its array/resource review remains scoped to instrumentation, without native equilibrium-response or G qualification. Fresh full Hessian and corrected checkpoint-reuse controls preserve the PCM/SMD H1 contribution. An unqualified restart that bypassed that term was explicitly rejected and stopped; it supplies no accepted result. See the [current assessment](../validation/screening_chemical_impact_2026-10-31/README.md).
