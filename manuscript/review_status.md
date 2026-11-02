# Current manuscript review disposition

As of 2 November 2026, this is an 18-page JCIM working article with a separate
31-page SI. Submission remains incomplete. Historical review dispositions are
preserved in Git; the source-frozen manuscript epoch is selected by
[evidence/current.json](evidence/current.json).

| Requirement | Current evidence | Still required |
| --- | --- | --- |
| M1: linked exact-candidate assessment | Exact viologen redox ladder and three-start local quaternization; TEMPO directional MeCN response; finite PMDI first-site/neighbor cycles; eight checked alternate-model polymer components | Complete the chosen lead's matched polymer G cycles, sampling and scoped decisions; no candidate currently passes every material requirement |
| M2: redox accuracy | 39 sourced endpoints from 30 molecules/seven studies; MAE 0.247 V, grouped diagnostics and reference handling | Retain retrospective/domain limits; no NDI validation or exact-graft experimental accuracy is inferred |
| M3: reaction convention | DBU/HCl and DIPEA source comparisons; four of five tested DBU sign screens change relative to HCl while pairwise ranks remain unchanged | Actual counterion/speciation and reaction-condition sensitivity; ΔG is not yield or rate |
| M4: polymer sampling | Forty alternate-model empty inputs meet the finite 10-to-20 criterion; primary finite PMDI/phthalimide cycles remain separately scoped | Occupied graph/configuration batches, duplicate review and balanced occupied ensemble cycles |
| M5: basis sensitivity | Five neutral-partner alternative-basis controls | Preserve limited scope; no complete basis-convergence claim |
| M6: numerical implementation | Native derivative/mode/partition checks; 71 benchmark electronic repeats; exact viologen/base-reaction panels; three full 119-atom Hessian controls | Other benchmark/base-reaction nuclear coverage and wider corrected production qualification |
| M7: reproducibility | Versioned exact-viologen CPU example; 75-state numerical cache; checked xTB component relocation; local reviewer collection and new control archive | Final synchronized code/data/environment release, author/license metadata and public deposition |
| M8: literature | Verified base reference audit and separately checked methodological/polymer additions | Final built bibliography and claim/source consistency |
| M9: publication narrative | Rewritten article; five quantitative interpretation-transfer tests; separate SI and focused slides | Integrate final lead results, demonstrate transferable scientific utility and synchronize final figures/data statements |

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

The alternate polymer collection now has eight independently checked components:
six charge-2 singlet single grafts and two matched neutral empty references.
Minimum internal frequencies range from 7.90 to 12.74 cm⁻¹. All retain their own
Hessian, nuclear, checkpoint, source and resource reviews. No complete ionic
occupied G cycle or converged ensemble is reported. The initial 137-atom pilot
stopped at an RKS H1 allocation limit, retaining its converged SCF, exact grids
and completed solute partial tensor. A separately declared H100 recovery reuses
only that frame/tensor, expands the RKS workspace to 16 GiB and recomputes full
solvent H1 and response under unchanged physical gates. Three corresponding
larger followups remain held for its independent review. The alternate model is
not mixed with primary molecular/polymer G values.

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
Its practical route is still on hold. The broader 904-graph registry and older
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
