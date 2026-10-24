# Exact PMDI site-3 native free-energy component

The [independent reconstruction](summary.json) verifies one 115-atom, five-repeat
PMDI-grafted polymer component with site 3 occupied. The open backbone has ten
carbons, nine backbone bonds and two ends; the graph and declared stereochemistry
pass the XYZ audit. This is a specified local occupancy scenario, not a measured
polymer loading or conformer ensemble.

| Quantity | Checked value |
| --- | ---: |
| Active SMD(MeCN) electronic energy | −94100.82627312276 eV |
| Electronic plus native nuclear/spin free energy, 50 cm⁻¹ cutoff | −94077.69081516036 eV |
| Internal modes | 339 |
| Lowest internal frequency | +5.439987 cm⁻¹ |
| Maximum gradient | 3.1115e-5 Ha/bohr |
| Maximum Hessian asymmetry | 5.0300e-11 |
| Maximum sampled Hessian-column relative error | 1.8356e-5 |

The normal native-state admission is rechecked. The independent calculation
reconstructs the projected Hessian/nuclear partition from stored records and
matches the admitted total G within 1e-8 eV. It also retains the 25/50/100 cm⁻¹
rotor-cutoff values. No xTB Hessian or substituted electronic-only result enters
this free energy.

**The matched PMDI cycle remains incomplete.** This single component does not
give an interaction free energy, neighbor penalty, conditional grafting G or
all-axis candidate pass. The site-2 and two-graft states still require completion;
broader conformer and configuration convergence is separate. Absolute component
energies are provided for reproduction, not as comparisons between molecules
with different compositions.

The [review implementation](producer_source.py) and all original state, Hessian,
geometry and calculation hashes are recorded in the summary. Raw state:
`calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_3_native_switch_curvature/thermal/baseline__design_pmdi_ethyl__occ_3__pose0/rest/`.
