# Current manuscript review disposition

As of 2 November 2026, this is a 21-page JCIM working article with a separate
35-page SI. Submission remains incomplete. Historical review dispositions are
preserved in Git; the source-frozen manuscript epoch is selected by
[evidence/current.json](evidence/current.json).

The [six-outcome chloride followup](../validation/counterion_completion_2026-11-03/README.md)
is now integrated into the article, SI Table S11 and slide 21. The finite 1-M
range is +5.33 to +22.36 kJ/mol, with five numerical duplicate groups and
positive values at all tested cutoffs. This expands the declared sample while
retaining conditional reaction and incomplete-speciation limits. The earlier
four-outcome editorial snapshot remains preserved.

## Next publication decisions

Keep the paper's current identity inventory fixed while completing a bounded
comparison on exact viologen, TEMPO and PMDI. Phenothiazine electrochemistry is
an optional extension requiring qualified state/reference components; its
unfinished controls do not block the core paper. Phthalimide remains a method control. Compare molecule-only,
route-aware local-graft and polymer-context interpretations using explicitly
matched inputs. Record changed priorities, supported conclusions and missing
evidence separately; an evidence hold is not a demonstrated synthetic failure.
The existing five interpretation-transfer tests are the starting evidence,
not proof of superior screening accuracy.

The next acceptance targets are complete independently checked lead-state
cycles, convergence/uncertainty for the occupancy comparisons supporting the
main conclusions, and an externally supported comparison with simpler
screening. Representative ion and reaction-condition sensitivity should test
the claimed decisions. Additional candidate enumeration is useful only when it
fills a declared coverage gap. Exact-graft experimental data would strengthen
application claims; the retrospective parent benchmark cannot establish
polymer performance. These priorities supplement the remaining requirements
below without declaring the current article submission-ready.

| Requirement | Current evidence | Still required |
| --- | --- | --- |
| M1: linked exact-candidate assessment | Exact viologen redox ladder and three-start local quaternization; TEMPO directional MeCN response; finite PMDI first-site/neighbor cycles; twelve checked alternate-model polymer components and all four declared matched neighboring-graft G cycles | Retain finite-model and conditional-reaction scope; qualify additional components only if their quantities enter a headline claim; full material admission is not claimed |
| M2: redox accuracy | 39 sourced endpoints from 30 molecules/seven studies; MAE 0.247 V, grouped diagnostics and reference handling | Retain retrospective/domain limits; no NDI validation or exact-graft experimental accuracy is inferred |
| M3: reaction convention | DBU/HCl and DIPEA source comparisons; four of five tested DBU sign screens change relative to HCl while pairwise ranks remain unchanged; six exact viologen/chloride outcomes pass full native QC and form five duplicate groups | Test the scenarios needed for the claimed interpretation, or retain conditional ΔG without a synthesis-feasibility conclusion; ΔG is not yield or rate |
| M4: polymer sampling | Earlier forty-input empty cohorts pass their finite diagnostic; the first corrected occupied 20-input cell fails the 2-kJ/mol diagnostic across all duplicate/rotor scenarios; SI Table S17 reports the signed changes | For equilibrium claims, complete and converge every required cycle component; otherwise restrict conclusions to retained geometries. Phthalimide failure does not establish viologen convergence |
| M5: basis sensitivity | Five neutral-partner alternative-basis controls | Preserve limited scope; no complete basis-convergence claim |
| M6: numerical implementation | Native derivative/mode/partition checks; 71 benchmark electronic repeats; exact viologen/base-reaction panels; three full 119-atom Hessian controls | Other benchmark/base-reaction nuclear coverage and wider corrected production qualification |
| M7: reproducibility | Versioned exact-viologen CPU example; 75-state numerical cache; checked xTB component relocation; local reviewer collection and new control archive | Final synchronized code/data/environment release, author/license metadata and public deposition |
| M8: literature | Verified base reference audit and separately checked methodological/polymer additions | Final built bibliography and claim/source consistency |
| M9: publication narrative | Rewritten article; five interpretation-transfer tests and a failed prospective occupied sampling diagnostic; concrete contribution in the introduction, results and conclusion; separate SI and focused slides | Lead with matched motif/chain and finite-backbone comparisons; keep the common-base shift secondary; reconcile final figures, evidence and release statements |

## Numerical reliability and its limits

The [implementation assessment](../validation/screening_chemical_impact_2026-10-31/README.md)
reproduces a GPU4PySCF screening defect on V100/H100 and checks an isolated
upstream correction. All 71 benchmark electronic states pass independent repeat
review; the largest reference-cancelled electronic potential shift is 8.66e-11 V.
Three [119-atom full-Hessian controls](../validation/polymer_hessian_reproduction_2026-10-31/README.md)
pass physical reconstruction and prospective comparisons; their largest nuclear
term difference is 1.67e-8 eV. These are scoped controls rather than broad
production or experimental qualification. New uncorrected admissions remain held.

The [new component panels](../validation/jcim_completion_2026-11-02/README.md)
cover ten exact viologen states and twenty-one base-reaction inputs, sharing two
Fc executions. All 70 actual SCF molecular frames are independently checked.
Corrected-minus-installed electronic contributions are at most 1.24e-11 V for
viologen potentials, 1.35e-9 kJ/mol for three local grafting sums and 4.18e-9
kJ/mol for five HCl-reference/twenty base-assisted sums. Retained nuclear terms
are not promoted to newly qualified total G by these electronic checks.

## Lead and oligomer evidence

The exact viologen's checked molecular reductions are −0.608/−1.208 V versus
matched Fc in MeCN. Local separated-ion quaternization is −11.49 kJ/mol at
298.15 K; three precursor starts span 1.50 kJ/mol and rotor settings give
−13.30 to −11.49 kJ/mol. This route has no HCl/base term. Salt association,
exchange, reflux conditions and polymer accessibility remain unresolved.
[Source/reproduction](../validation/viologen_grafting_2026-10-31/README.md).

The [ionic polymer collection](../validation/viologen_native_2026-11-02/README.md)
now contains twelve checked stationary components and four balanced
whole-chain electronic nonadditivities: adjacent/one-spacer values are
−3.46/+18.68 kJ/mol for all-same and +1.20/−13.48 kJ/mol for alternating
configuration. They are single-conformer E terms under one corrected
r2SCAN-D4/SVP smooth-surface MeCN model. They supply neither interaction G,
isolated motif binding nor operating-state coupling.

The alternate polymer collection now has twelve independently checked components
and all four declared finite G cycles: all-same adjacent −3.94 kJ/mol, all-same
one-spacer +22.58, alternating adjacent +8.28, and alternating one-spacer −5.30
at the primary cutoff. Both finite signs reverse with backbone configuration.
Each double-graft component has a full Hessian and all 25 actual SCF frames checked.
The alternating adjacent recovery uses its original source-checked full matrix
and response, all 25 fresh center/displaced evaluations, and independent own-mode
and nuclear reconstruction. No ensemble or material inference follows.
[Current collection](../validation/ionic_cycle_completion_2026-11-02/complete_cycles_12.json).

The separate full 240-input corrected phthalimide-control campaign now has
90 checked electronic inputs and eighteen checked native inputs. The largest
119-atom pilot passes full independent minimum, derivative, wavefunction and
resource review; matching native followups continue under source-bound resource and physical gates. No occupied ensemble is complete. Six exact local
viologen/chloride geometries now pass independent stationarity and actual
SCF-energy checks. The first pose's −50.52-cm⁻¹ minimum failure is independently
reconstructed and preserved unqualified. Four checked outcomes from three original starts pass full native QC and give finite
1-M association G of +15.4 to +22.4 kJ/mol at the primary cutoff. Two other
failed minima have diagnosed signed repair branches prepared. Other poses and
actual-condition speciation remain open. Twelve frozen
motif triplets additionally pass 72 actual SCF reviews; local repulsion and the
whole-chain E cycle have opposite signs in two of four parent frames. These
motif and six chloride-contact outcomes are integrated into the current frozen manuscript and slide evidence. The N0 repairs share one numerical duplicate group; no hit-count entropy is introduced.
[Sampling progress](../validation/corrected_sampling_2026-11-02/README.md).

The [full primary molecular cohort](../validation/ionic_cycle_completion_2026-11-02/README.md)
now repeats all ten exact viologen reaction/reference states with nine new
molecular Hessians and independent analytic chloride translation. Eighty-one
native SCF frames and occupied-orbital spin values pass reconstruction; 70
additional electronic-control wavefunctions also pass independent spin review.
The six complete cycles change by at most 2.44e-8 V or 1.54e-6 kJ/mol from the
source across all cutoffs. This qualifies numerical reproduction of the tested
primary cohort, without experimental/speciation or occupied-ensemble admission.

TEMPO retains directional solution costs 1.477/1.639 eV. Their 10.4% asymmetry
fails the declared 5% equal-curvature diagnostic; a symmetric scalar remains
unqualified. PMDI retains a finite primary-model −1.86-kJ/mol four-state G
cycle and +5.3–6.1-kJ/mol first-site penalties, with one minimum per graph.
Its practical route is still on hold. The broader 912-graph registry and older
878-candidate all-axis snapshot have distinct scopes; neither is a material shortlist.

## Reviewer access and completion

The new 29-execution local control archive contains 70 SCF checkpoints. Every
archive byte is verified and an extracted CPU replay succeeds with original
project/data paths blocked and no GPU-library imports. It reconstructs cached
frames/energies/cycle shifts, not fresh quantum calculations or wavefunction
spin diagnostics. Public deposition and final metadata remain open.

The prior measured Nitro repeatability failure, alternative-spin limits and
incomplete occupied ensembles remain visible. A successful later repeat does
not erase the failure or establish force-step convergence. Compiler, visual
inspection and self-contained export records each retain separate scopes;
no Overleaf remote is connected. Author details and final data reuse terms
require the authors' input. The [publication acceptance criteria](../docs/PUBLICATION_READINESS.md)
remain authoritative; green tests or a built PDF do not establish completion.
