# Empty-chain finite-set sampling sensitivity

[Results](summary.json) · [Independent arithmetic and negative check](independent_check.json) · [Tests](tests.xml) · [Four-state geometry review](../empty_chain_finite_set_2026-10-26/README.md)

For the four independently reviewed empty five-repeat structures, we exhaustively evaluated every nonempty proper subset. All four use one backbone configuration and the alternative smooth-surface acetonitrile model. At the 50 cm⁻¹ rotor cutoff, the six two-structure subsets give finite-set partition free energies **0.60–3.82 kJ/mol above** the full four-structure partition. The worst two-structure subset misses the full-set minimum by **2.34 kJ/mol**. Three-structure subsets differ by 0.07–1.63 kJ/mol in their partitions. At the 25 and 100 cm⁻¹ cutoffs, the worst three-structure differences are 2.49 and 2.06 kJ/mol, respectively.

These are comparisons within a known four-structure set, not statistical confidence intervals or upper bounds on unsampled structures. The check says nothing about convergence of occupied graphs, a different backbone configuration, the primary solvent model, inter-motif interactions or a polymer cycle. The conditional finite partition also carries the qRRHO conformer/mode overlap limitation described in the [four-state review](../empty_chain_finite_set_2026-10-26/README.md). No new quantum calculation or physical QC was performed.

The producer checks the source review hashes, shared graph/configuration/method, distinct geometry groups, and full-partition consistency before evaluating subsets. The independent check reconstructs 42 subset differences from conditional weights and rejects a tampered duplicate grouping.

Reproduce with `PYTHONPATH=src python scripts/validation/empty_conformer_subset_sensitivity.py --manifest results/validation/empty_chain_native_range_2026-10-26/manifest.json --finite results/validation/empty_chain_finite_set_2026-10-26/summary.json --out NEW_OUTPUT.json`.
