# Exact PMDI first-site grafting: checked sampled minima

[Independent reconstruction](summary.json) · [Generic analysis](../../../../scripts/validation/analyze_first_site_grafting.py) · [Main manuscript](../../../../manuscript/manuscript.pdf)

The exact mono-N-ethyl/NH PMDI precursor is attached at site 2 or 3 of a five-repeat open polymer model. Sites 2–4 are initially reactive benzyl chlorides; the two end repeats are plain styrene. Every reaction is independently balanced for atoms and charge, and each nuclear partition is reconstructed from its stored modes.

| Sampled site | DBU, separated ions | DBU, contact pair | Polymer minus capped local model |
| --- | ---: | ---: | ---: |
| 2 | −97.72 | −102.52 | +5.34 |
| 3 | −96.96 | −101.75 | +6.11 |

Values are kJ/mol at 298.15 K in SMD acetonitrile, primary 50 cm⁻¹ rotor setting. The contact-pair reaction includes the −7.92587 kJ/mol, Δn=−1 standard-state term exactly once. The new calculation uses the SI-derived eV conversion; frozen legacy arithmetic remains unchanged.

Common precursor, base and byproduct terms cancel in the polymer-minus-local double-mutation cycle. Its positive penalties persist at the prescribed 25/50/100 cm⁻¹ settings. The small difference between sites changes sign, so these results do not establish a preferred site. The shared empty-chain reference is the lower-electronic-energy member of the two inherited exact empty graphs, not a globally converged ensemble.

The independently checked paired-state electronic cycle is −2.30582 kJ/mol. **Paired-state G and neighboring-graft ΔΔG remain blank.** Electronic energy is never substituted for free energy. Proton-transfer/speciation uncertainty still qualifies absolute DBU reaction signs. Barriers, selectivity, post-grafting redox retention and material utilization remain unverified; no candidate becomes synthesis-eligible.

Raw single-site records: `calcs/polymer_models/2026-10-24_streamed_native/` and `calcs/polymer_models/2026-10-23_crest_refinement/`. Exact state paths, input/output hashes, empty/capped/base sources and all cutoff terms are stored in the JSON. The immutable manuscript copy and analysis producer are in `manuscript/evidence/2026-10-25_first_graft_checked/`.
