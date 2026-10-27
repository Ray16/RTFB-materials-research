# Speaker notes: molecular design to grafted polymers

1. Motivation: introduce the molecular-to-material question. The displayed PMDI chain is an actual computed structure, not an assumed bead population.

2. Motivation: larger carriers motivate size-selective separation. The cartoon is conceptual; the two cited demonstrations are context, not validation of this proposed resin.

3. Question: require useful redox behavior, a feasible attachment and weak unwanted intersite interactions. Avoid a compensating weighted score.

4. Pipeline: define exact precursor/product and polymer graphs, screen new routes, generate RDKit conformers and use UMA for gas-phase warm starts. DFT in MeCN supplies the reported redox, grafting and neighbor quantities. Nitro qualification remains open; do not describe the standard V100 fallback as Nitro.

5. Expansion: 837 designs are enumerated, not all computed. The evaluated 75-candidate cohort is distinct; 60 have both potentials. Coverage is strongest in three reduction-active cores.

6. Molecular result: the shown PMDI ring-F and ring-CN shifts are +0.158 and +0.177 V. Cap and ether-linker shifts are −0.022 and −0.045 V. These are paired model differences, not externally validated improvements.

7. Validation result: primary unfitted MAE is 0.247 V and RMSE 0.279 V for 39 endpoints. The panel is retrospective and quinone-heavy, with no NDI endpoint. Do not interpret MAE as a per-candidate confidence interval.

8. Grafting result: BQ-OMe changes from +57.6 to −45.1 kJ/mol and BQ-Me3 from +56.0 to −46.7 under the computed DBU separated-ion route. The shared base offset leaves ranks unchanged. The experimental MeCN pKa comparison shows an approximately 23 kJ/mol discrepancy; thermodynamic signs and synthetic feasibility remain conditional.

9. Polymer method: explain G00, G10, G01 and G11 on matched open backbones before showing their balanced difference. Cl marks a reactive benzyl-chloride handle, not ordinary styrene. No experimental loading distribution is inferred.

10. PMDI result: the finite cycle is −1.86 kJ/mol with 16.30 Å core separation. First-site polymer penalties are +5.34 and +6.11 kJ/mol relative to the local capped model. This is whole-chain nonadditivity from one minimum per graph, not pure motif binding or global ensemble convergence.

11. Interaction result: close phthalimide control poses are approximately −25 to −26 kJ/mol in MeCN; separated poses are near zero. The plot uses frozen capped fragments. The cap-length variation is not a statistical confidence interval. No BSSE correction, association G or electronic coupling is established.

12. Selection: eight models lie on the molecular three-objective descriptor front, whose projection is shown. None has demonstrated every material requirement. Nominal molecular capacity is not usable polymer capacity.

13. Closing: tune redox chemistry systematically, evaluate a defined reaction and test local polymer contexts. Broader chemistry and complete occupied-state evidence are next scientific requirements.

14. Backup: first/second-step errors differ. Offset statistics are molecule-held-out internal diagnostics, not replacements for unfitted primary statistics.

15. Backup: distinguish experimental property accuracy, internal numerical consistency and model/sampling robustness.

16. Backup: the alternate-model empty-backbone 10-to-20 diagnostic passes for one configuration. This does not complete all occupied graph/configuration comparisons or validate global sampling.

17. Backup: raw-state packages and scoped relocation replay exist; public deposition and full material assessment remain open. UMA installation and numerical qualification are separate.
