"""Decorate the Merrifield-monomer benzylic site with redox groups and build 3D structures.

Reads the scaffold + groups from config/redox_groups.py, zips each group onto the
scaffold, enumerates redox states (charge/spin), embeds a 3D conformer, and writes:

  library/<group_id>/<group_id>.smiles     canonical decorated SMILES (state-independent)
  library/<group_id>/<state>.xyz           3D geometry (same geom seed per state)
  library/manifest.csv                     id, name, family, state, charge, mult, n_e, smiles

Run:  python -m redox.build          (from src/, with the `redox` env active)
"""
from __future__ import annotations
import csv
import importlib.util
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "redox_groups.py"
LIBRARY = ROOT / "library"


def _load_config():
    spec = importlib.util.spec_from_file_location("redox_groups", CONFIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.SCAFFOLD, mod.GROUPS


def decorate(scaffold_smi: str, frag_smi: str) -> Chem.Mol:
    """Zip a fragment onto the scaffold at matching [*:1] dummies."""
    scaffold = Chem.MolFromSmiles(scaffold_smi)
    frag = Chem.MolFromSmiles(frag_smi)
    if scaffold is None or frag is None:
        raise ValueError(f"bad SMILES: {scaffold_smi!r} / {frag_smi!r}")
    combined = Chem.CombineMols(scaffold, frag)
    mol = Chem.molzip(combined)          # joins the two [*:1] dummies into a bond
    Chem.SanitizeMol(mol)
    return mol


def unassigned_stereo(mol: Chem.Mol) -> int:
    """Count stereo elements (centers/bonds) left UNSPECIFIED. Non-zero means the 3D
    embedding would pick an arbitrary isomer — a hard guard for library expansion."""
    from rdkit.Chem import FindPotentialStereo, StereoSpecified
    return sum(1 for e in FindPotentialStereo(mol)
               if e.specified != StereoSpecified.Specified)


def n_conformers(mol: Chem.Mol) -> int:
    """Heuristic ensemble size from rotatable-bond count."""
    from rdkit.Chem import Descriptors
    n_rot = Descriptors.NumRotatableBonds(mol)
    if n_rot <= 3:
        return 25
    if n_rot <= 6:
        return 75
    return 150


def radius_of_gyration(molH: Chem.Mol) -> float:
    """Unweighted radius of gyration (Å) — a cheap collapse indicator."""
    import numpy as np
    conf = molH.GetConformer()
    xyz = np.array([list(conf.GetAtomPosition(i)) for i in range(molH.GetNumAtoms())])
    return float(np.sqrt(((xyz - xyz.mean(0)) ** 2).sum(1).mean()))


def best_conformer(mol: Chem.Mol, net_charge: int = 0,
                   seed: int = 0xC0FFEE, prune_rms: float = 0.5):
    """Conformer search: embed an ETKDG ensemble, FF-rank, RMSD-prune, pick lowest.

    ETKDG (distance geometry) yields extended seeds. For charged species we deliberately
    do a LIGHT FF cleanup only (few iters), to avoid gas-phase over-folding into
    artificial salt bridges / intramolecular H-bonds; those states are re-optimized in
    implicit solvent (xtb-ALPB / DFT-SMD) downstream. Returns
    (mol_best_conf, ff_used, energy, n_kept, rg)."""
    molH = Chem.AddHs(mol)
    n = n_conformers(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = prune_rms
    params.useRandomCoords = False
    cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n, params=params))
    if not cids:  # retry with random coords for hard cases
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n, params=params))
    if not cids:
        raise RuntimeError("conformer embedding failed")

    # Gentle for charged species (avoid gas-phase collapse); fuller for neutral.
    maxiters = 200 if net_charge != 0 else 1000
    ff_used = "MMFF"
    props = AllChem.MMFFGetMoleculeProperties(molH)
    if props is not None:
        res = AllChem.MMFFOptimizeMoleculeConfs(molH, maxIters=maxiters)
    else:
        ff_used = "UFF"
        res = AllChem.UFFOptimizeMoleculeConfs(molH, maxIters=maxiters)
    energies = [e for _, e in res]

    best = min(range(len(cids)), key=lambda i: energies[i])
    out = Chem.Mol(molH)
    out.RemoveAllConformers()
    out.AddConformer(molH.GetConformer(cids[best]), assignId=True)
    return out, ff_used, energies[best], len(cids), radius_of_gyration(out)


def to_xyz(molH: Chem.Mol, comment: str) -> str:
    conf = molH.GetConformer()
    lines = [str(molH.GetNumAtoms()), comment]
    for atom in molH.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return "\n".join(lines) + "\n"


def main():
    scaffold, groups = _load_config()
    LIBRARY.mkdir(exist_ok=True)
    rows = []
    for g in groups:
        gdir = LIBRARY / g["id"]
        gdir.mkdir(exist_ok=True)
        mol = decorate(scaffold, g["frag"])
        smi = Chem.MolToSmiles(mol)
        n_unspec = unassigned_stereo(mol)
        if n_unspec:
            print(f"  WARNING {g['id']}: {n_unspec} UNSPECIFIED stereo element(s) — "
                  f"embedding picks an arbitrary isomer; specify stereo in the SMILES "
                  f"or enumerate. SMILES={smi}")
        (gdir / f"{g['id']}.smiles").write_text(smi + "\n")
        # One conformer search per molecular framework -> shared seed geometry.
        net_q = Chem.GetFormalCharge(mol)
        molH, ff, e, nkept, rg = best_conformer(mol, net_charge=net_q)
        seed = to_xyz(molH, f"conformer seed ff={ff} E={e:.3f} nconf={nkept} Rg={rg:.2f} smiles={smi}")
        (gdir / "conformer.xyz").write_text(seed)
        print(f"  {g['id']:16s} seed ff={ff:4s} E={e:9.2f} nconf={nkept:3d} Rg={rg:.2f} q={net_q:+d}")
        # Per-state seeds start from the same conformer; UMA/DFT relax each state.
        # Charged states get solvated pre-opt (xtb-ALPB) to avoid gas-phase artifacts.
        for label, charge, mult, n_e in g["states"]:
            (gdir / f"{label}.xyz").write_text(
                to_xyz(molH, f"charge={charge} mult={mult} n_e={n_e} smiles={smi}"))
            rows.append(dict(id=g["id"], name=g["name"], family=g["family"],
                             state=label, charge=charge, mult=mult, n_e=n_e,
                             solv_preopt=int(charge != 0), smiles=smi))
            print(f"      {label:5s} q={charge:+d} m={mult} n_e={n_e:+d}"
                  f"{'  [solv-preopt]' if charge != 0 else ''}")
    manifest = LIBRARY / "manifest.csv"
    with manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "name", "family", "state",
                                          "charge", "mult", "n_e", "solv_preopt", "smiles"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} states across {len(groups)} groups -> {manifest}")


if __name__ == "__main__":
    sys.exit(main())
