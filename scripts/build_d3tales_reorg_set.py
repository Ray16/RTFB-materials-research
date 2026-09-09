#!/usr/bin/env python
"""Extract all quinone + imide molecules (our focus families that D3TaLES actually contains;
viologens are absent from D3TaLES) with reorganization-energy values, to validate our machinery
against D3TaLES under their protocol. Writes results/d3tales_reorg_validation/molecules.csv.

Only neutral-ground-state molecules are kept (groundState_charge == 0) so the electron couple is
neutral<->anion and the hole couple is neutral<->cation, matching our 4-point definition.
"""
from __future__ import annotations
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
DUMP = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
OUT = ROOT / "results" / "d3tales_reorg_validation"

FAM = {
    "quinone": ["O=C1C=CC(=O)C=C1", "O=C1c2ccccc2C(=O)c2ccccc21", "O=C1c2ccccc2C(=O)C=C1"],
    "imide":   ["O=C1[#6]~[#6]C(=O)[#7]1", "O=C1[#6]=[#6]C(=O)[#7]1", "O=C1c2ccccc2C(=O)[#7]1"],
}
PATS = {f: [Chem.MolFromSmarts(s) for s in ss] for f, ss in FAM.items()}


def family_of(m):
    for fam, ps in PATS.items():
        if any(m.HasSubstructMatch(p) for p in ps if p is not None):
            return fam
    return None


def main():
    cols = ["_id", "smiles", "groundState_charge", "sa_score",
            "hole_reorganization_energy", "electron_reorganization_energy",
            "relaxation_groundState_anion1", "relaxation_anion1_groundState",
            "relaxation_groundState_cation1", "relaxation_cation1_groundState"]
    df = pd.read_csv(DUMP, usecols=cols)
    rows = []
    for _, r in df.iterrows():
        s = r["smiles"]
        if not isinstance(s, str):
            continue
        if pd.notna(r["groundState_charge"]) and int(r["groundState_charge"]) != 0:
            continue  # keep neutral-ground-state only
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        fam = family_of(m)
        if fam is None:
            continue
        has_e = pd.notna(r["electron_reorganization_energy"])
        has_h = pd.notna(r["hole_reorganization_energy"])
        if not (has_e or has_h):
            continue
        rows.append(dict(id=r["_id"], smiles=s, family=fam,
                         d3_electron=r["electron_reorganization_energy"],
                         d3_hole=r["hole_reorganization_energy"],
                         d3_relax_gs_anion=r["relaxation_groundState_anion1"],
                         d3_relax_anion_gs=r["relaxation_anion1_groundState"],
                         d3_relax_gs_cation=r["relaxation_groundState_cation1"],
                         d3_relax_cation_gs=r["relaxation_cation1_groundState"],
                         n_atoms=m.GetNumAtoms(), n_rot=rdMolDescriptors.CalcNumRotatableBonds(m)))
    out = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    fp = OUT / "molecules.csv"
    out.to_csv(fp, index=False)
    print(f"wrote {fp}  ({len(out)} molecules)")
    print("\nby family:")
    print(out.groupby("family").agg(n=("id", "size"),
                                    e_reorg=("d3_electron", lambda x: x.notna().sum()),
                                    h_reorg=("d3_hole", lambda x: x.notna().sum())).to_string())
    print(f"\nsize distribution (heavy atoms): median={out.n_atoms.median():.0f} "
          f"max={out.n_atoms.max():.0f}  (>40 atoms: {(out.n_atoms>40).sum()})")
    print(f"rotatable bonds: median={out.n_rot.median():.0f} max={out.n_rot.max():.0f}")


if __name__ == "__main__":
    main()
