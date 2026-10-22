# Broader conformer revision: active CPU pilots

Two empty-backbone CREST 3.0.2 iMTD-GC searches use GFN2/ALPB(acetonitrile), an 8 kcal/mol energy window and two disjoint CPU slots each. The all-same and alternating mapped stereochemical configurations are distinct constructed scenarios; neither is asserted to represent an atactic ensemble. [Launch](launch.json) and [actual process verification](start_verification.json) record ownership, affinity and six native caps=2.

Raw durable trajectories/logs: `calcs/polymer_models/2026-10-22_crest_all_same/` and `calcs/polymer_models/2026-10-22_crest_alternating/`. Follow their status.json for actual progress. CREST binary and inputs are hash-pinned in each status. The installed 3.0.2 CLI has no documented iMTD RNG-seed override; the program default and realized trajectories are retained. Initial alternate-configuration embedding uses seed 20261022. [environment_explicit.txt](environment_explicit.txt) locks this sampling environment separately from the primary DFT/qRRHO environment.

Only the two empty-backbone pilots have started. Other occupancy graphs, DFT re-ranking, minimum-convergence evidence and ensemble G are pending. A low-level ensemble is never promoted directly to a DFT result. Producer: `scripts/pipeline/sample_polymer_crest.py`; it checks exact topology/stereochemistry and rejects truncated XYZ frames.
