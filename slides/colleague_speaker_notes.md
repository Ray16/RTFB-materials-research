# Speaker notes for the focused workflow presentation

Suggested length: 18–22 minutes, with backup used only for questions. Slides 3–6
cover expansion; slides 7–10 cover validation. Keep those sections if time is short.

| Slide | Purpose | What to say |
| --- | --- | --- |
| 1 | Establish the scientific question | “We are selecting groups for an existing polymer. A good free molecule is a starting point; we need its advantages to survive grafting and nearby occupied sites.” |
| 2 | Show how the components work together | “We preserve exact precursor/product identity and the attachment bond. We review the route before an expensive new redox campaign. Different calculations then contribute evidence to the same design decision.” |
| 3 | Define the expansion denominator | “837 is the enumerated graph library, not a count of completed quantum campaigns. There are 1,593 recipes because several recipes can produce the same graph. The evaluated cohort overlaps 33 graphs; 32 have both potentials.” |
| 4 | Make a design axis concrete | “These are actual calculated PMDI structures, not generic sketches. The ethyl cap and graft-site model are held fixed while the declared ring substitution changes. Each has its own precursor graph.” |
| 5 | Show a result from systematic expansion | “The circles and squares are first and second potential shifts relative to a matched baseline. Ring electronics gives a larger tuning range than these cap changes. More positive or negative is not automatically better; it depends on the electrode and operating window.” |
| 6 | Explain the next expansion scientifically | “A large three-core library does not solve core diversity. The next plan deliberately compares seven cores and two polymer-relevant design factors. The 34 entries are slots, not new completed molecules. Shielding and constraints can also hinder grafting or favor persistent contacts.” |
| 7 | Define what ‘validated’ means | “We separate arithmetic/numerical consistency, comparison against experiment, and sensitivity to model choices. A correct sum or converged SCF is necessary but cannot establish chemical accuracy.” |
| 8 | State the measured accuracy honestly | “All 39 endpoints now finish. The primary MAE is about a quarter volt, with no fitted offset and all outliers retained. This is a retrospective free-parent MeCN benchmark. It does not calibrate exact grafts or polymer association, and it has no NDI endpoint.” |
| 9 | Explain the bias diagnosis | “Every first-step error is positive, but the second-step mean is slightly negative. Held-out offset diagnostics show structure in the error. We retain the unfitted benchmark; fitting a constant does not prove Fc is the cause, and a shared shift cannot fix step spacings.” |
| 10 | Quantify specific robustness checks | “We tested the geometry/Hessian basis imbalance on five neutral partners. We also varied the rotor cutoff. These shifts are tens of millivolts, but they are settings sensitivities, not confidence intervals or bounds on missing conformers.” |
| 11 | Prevent the Pareto plot becoming a synthesis claim | “The rings mark eight leaders in three molecular descriptors. The color shows the second potential, which is also part of the dominance comparison. None is currently all-axis eligible. Nominal molecular capacity omits real loading and accessibility.” |
| 12 | Explain the base-assisted result | “A fixed base model adds the same term to every eligible reference reaction, so ranks do not change. BQ should not be rejected because a different HCl-reference reaction is positive. The capture-term pKa discrepancy and uncertain speciation keep the absolute signs conditional, especially for DIPEA.” |
| 13 | Demonstrate an exact end-to-end case | “The same PMDI identity now connects redox, exact precursor access, a defined base-assisted reaction and polymer occupancy. This is a concrete worked example. It still does not establish selective attachment, redox-core retention or synthesis eligibility.” |
| 14 | Define the polymer observable precisely | “The four-state cycle asks how the modeled cost of attaching one group changes when another is present. The sampled PMDI pair is separated by 16.3 Å. A small cycle value in this structure says nothing by itself about isolated binding, charge-state coupling or all conformers.” |
| 15 | Make progress since the reviewed draft visible | “We closed the missing endpoint, completed the finite baseline thermochemistry and added the PMDI cycle and portable raw-state replay. The important remaining calculation is a defensible broader ensemble, rather than more arithmetic on one structure.” |
| 16 | Turn the discussion into concrete next decisions | “We should agree on application criteria, prioritize exact route-feasible designs that fill chemistry gaps, and establish sufficient conformer evidence before drawing polymer spacing rules. Those steps turn descriptor leaders into a defensible material shortlist.” |

Numerical answers to likely questions:

- Four checked empty-chain free energies span 7.19 kJ/mol at the 50 cm⁻¹
  cutoff in the separate smooth-solvent model. This is not a grafting-cycle
  uncertainty or an ensemble, and should not be mixed with the primary model.

- A separate AQ-2-carboxylate graft has a source-linked precursor and mapped
  attachment, but remains on route hold with no computed redox properties. It
  is not counted among the 837 earlier designs or the 34 unfilled comparison slots.

- The evaluation feed is dated **24 October**, while the exact PMDI result and
  reproducibility addendum were checked on **25–26 October**. The deck does not
  pretend they are one simultaneous live campaign snapshot.
- The 879-graph decision ledger contains **878 pending candidate designs and
  one comparison control**. The control is already represented in the graph
  library; do not add it again to the union count.
- The 75-candidate cohort has **62 candidates with at least one native-calculation
  potential** and **60 with both**. This does not mean every candidate has every
  material property.
- Unfitted benchmark MAE/RMSE/bias: **0.24697/0.27929/+0.18959 V**.
  Molecule-held-out one-offset MAE/RMSE: **0.15833/0.21084 V**; step-offset:
  **0.10605/0.16605 V**. These are internal diagnostic estimates.
- DBU/HCl experimental MeCN pKa values in the checked source analysis are
  **24.31/10.30**. The computed capture term differs by **22.81 kJ/mol**.
  HCl cancels from the full reaction cycle, so this discrepancy cannot simply
  be added to every complete grafting reaction as a validated correction.
- Cook et al. report an NMR pKa estimate of **18.7** in MeCN for the DIPEAH
  hexafluorophosphate salt (DOI 10.1039/D3SC02611F). The conditional cross-source
  DIPEA/DBU transfer comparison is +31.80 kJ/mol computed versus +32.02 kJ/mol
  inferred from the pKa values. The salt and chloride/activity conditions are
  not matched; this does not establish a correction to the complete grafting G.
- Barišić et al. independently report **log Kᴴ = 18.64(1)** for DIPEA in
  MeCN at 25 °C using bromocresol green (DOI 10.1039/D1RA04738H, SI Fig. S51).
  Combining it with the separate DBU scale gives +32.36 kJ/mol versus +31.80
  computed. Indicator, salt and activity conventions remain important.
- PMDI rotor-cutoff neighbor values at 25/50/100 cm⁻¹ are
  **−1.07/−1.86/−1.82 kJ/mol**. The sampled phthalimide two-pose adjacent/spacer
  results are **−5.17/−11.39 kJ/mol**, under a different exact model. Neither
  panel establishes conformer convergence.
- Native thermochemistry uses DFT Hessian modes. xTB serves as the qRRHO
  integrator; it does not supply an xTB Hessian to those free energies.
- Nitro is installed and preferred by the public entry point where supported.
  Actual redox-panel Nitro parity and matched-device speedup remain unqualified.
  Standard V100 inference is separately tested. Do not claim a Nitro speedup.

For detailed questions, use the seven backup slides. The full historical deck
is preserved in [the earlier-deck archive](archive/2026-10-26_before_focused_workflow/README.md).

Sources: [checked discussion evidence](../results/reviews/2026-10-26/colleague_workflow/evidence.json),
[review response ledger](../manuscript/revision_plan.md),
[candidate expansion plan](../docs/TARGETED_CANDIDATE_EXPANSION.md),
[benchmark](../results/validation/reference_bias_2026-10-22/README.md),
[basis controls](../results/validation/neutral_basis_2026-10-22/README.md),
[PMDI cycle](../results/reviews/2026-10-25/pmdi_neighbor_cycle_checked/README.md),
[raw replay](../results/releases/neighbor_cycle/README.md).
