# Finite empty-chain conformer set

[Checked result](summary.json) · [Frozen input reviews](../empty_chain_native_range_2026-10-26/manifest.json) ·
[Partition tests](tests.xml) · [Independent check](independent_check.json) ·
[Earlier range diagnostic](../empty_chain_native_range_2026-10-26/README.md)

Four previously reviewed native free-energy records describe the same empty
five-repeat backbone and one configuration in the alternative smooth-surface
acetonitrile model. Each geometry was checked against its declared connectivity
and stereochemistry and its source hashes. Symmetry-aware heavy-atom pairwise
RMSDs span **2.01–4.27 Å**. The four structures remain separate under each
declared duplicate threshold of 0.05, 0.10 and 0.20 Å, together with the
0.20 kJ/mol electronic-energy criterion.

At 298.15 K and the 50 cm⁻¹ rotor cutoff, a conditional partition over these
four structures is **1.81 kJ/mol below the lowest individual G**. At 25 and
100 cm⁻¹ the corresponding lowerings are 1.13 and 1.42 kJ/mol. Each distinct
minimum is counted once; search-hit frequency is not used as a statistical
weight. The value is an approximate finite-set model, not an experimentally
measured population. Soft backbone modes and conformer partitioning may overlap
in the approximate qRRHO treatment.

These four structures do not establish a converged ensemble, quantify a polymer
grafting cycle, compare tacticities, or measure motif association or electronic
coupling. This alternate model cannot be mixed with the primary-model reaction
cycle. The analysis checks the geometry, plan and electronic-review hashes bound
by the independent native review; it does not repeat the quantum calculations
or all the reviews' raw-record checks. A separate arithmetic reconstruction
reproduced the 1.81 kJ/mol lowering, and a deliberately altered review hash was
rejected before output. The partition and input-gating tests passed 22 cases.
