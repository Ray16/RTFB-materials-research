# Validation datasets (external, git-ignored)

Re-clone (see docs/DATASETS.md for details):

    git clone https://github.com/Liu-group/SI_data_redox_paper   # OROP/ROP313 — experimental MeCN potentials
    git clone https://github.com/grynova-ccc/ReSolvedDB          # computed MeCN reduction potentials (RP_ACN) + SMILES

- OROP: `implicit_solvation_results.csv` (exp. + computed, multi-functional benchmark:
  B3LYP-D3, CAM-B3LYP-D3, ωB97X-D3, PBE0-D3, ωPBEh-D3). Molecules indexed by number;
  structures in `features/` + `implicit_optimized_xyz/`.
- ReSolvedDB: `ReSolvedData.csv` — SMILES, RP_ACN, HOMO/LUMO, EA, stability, per-solvent
  solvation free energies. Use for computed cross-check + candidate expansion.
