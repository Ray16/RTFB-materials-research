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
    """Count GENUINE stereo elements (centers/bonds) left UNSPECIFIED. Non-zero means the 3D
    embedding would pick an arbitrary isomer — a hard guard for library expansion.

    Excludes symmetric hypervalent centers that have no real stereoisomers: RDKit's
    FindPotentialStereo over-eagerly flags e.g. octahedral PF6- (six identical F) as an
    unspecified stereocenter, but with all ligands symmetry-equivalent there is only one
    isomer. We drop any atom-centered element whose neighbors are all equivalent (identical
    canonical ranks); genuine centers have distinguishable substituents and are still counted.
    """
    from rdkit.Chem import FindPotentialStereo, StereoSpecified, CanonicalRankAtoms
    ranks = list(CanonicalRankAtoms(mol, breakTies=False))
    n = 0
    for e in FindPotentialStereo(mol):
        if e.specified == StereoSpecified.Specified:
            continue
        if e.type.name.startswith("Atom"):
            nbr_ranks = [ranks[nb.GetIdx()]
                         for nb in mol.GetAtomWithIdx(e.centeredOn).GetNeighbors()]
            if nbr_ranks and len(set(nbr_ranks)) == 1:
                continue  # all ligands identical -> not a real stereocenter
        n += 1
    return n


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


def conformer_ensemble(mol: Chem.Mol, net_charge: int = 0, top_k: int = 10,
                       seed: int = 0xC0FFEE, prune_rms: float = 0.1):
    """Conformer search: embed an ETKDG ensemble, FF-rank, RMSD-prune, return the TOP-K
    lowest-FF-energy distinct conformers (each re-ranked per redox state downstream with
    UMA, which is charge/spin-aware; FF ranking here is only a coarse pre-filter).

    Tighter pruning (0.1 A) keeps benzylic/linker rotamers the old 0.5 A collapsed. For
    charged species we do a LIGHT FF cleanup only, to avoid gas-phase over-folding into
    artificial salt bridges; final geometries come from UMA + DFT-SMD. Returns
    (molH_with_topk_confs, conf_ids_sorted, ff_used, energies_sorted)."""
    molH = Chem.AddHs(mol)
    n = n_conformers(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = prune_rms
    params.useRandomCoords = False
    cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=n, params=params))
    if not cids:
        raise RuntimeError("conformer embedding failed")

    maxiters = 200 if net_charge != 0 else 1000
    ff_used = "MMFF"
    props = AllChem.MMFFGetMoleculeProperties(molH)
    if props is not None:
        res = AllChem.MMFFOptimizeMoleculeConfs(molH, maxIters=maxiters)
    else:
        ff_used = "UFF"
        res = AllChem.UFFOptimizeMoleculeConfs(molH, maxIters=maxiters)
    energies = [e for _, e in res]

    order = sorted(range(len(cids)), key=lambda i: energies[i])[:top_k]
    return molH, [cids[i] for i in order], ff_used, [energies[i] for i in order]


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
        # Conformer ENSEMBLE per molecular framework -> conf_00.xyz .. conf_KK.xyz.
        # These are charge/spin-independent geometries; UMA re-ranks them per redox state
        # (charge/spin-aware) and scans spin multiplicities downstream (see uma.py).
        net_q = Chem.GetFormalCharge(mol)
        confdir = gdir / "conformers"; confdir.mkdir(exist_ok=True)
        molH, conf_ids, ff, energies = conformer_ensemble(mol, net_charge=net_q)
        for f in confdir.glob("conf_*.xyz"):
            f.unlink()  # clear stale ensemble
        for i, (cid, e) in enumerate(zip(conf_ids, energies)):
            single = Chem.Mol(molH); single.RemoveAllConformers()
            single.AddConformer(molH.GetConformer(cid), assignId=True)
            rg = radius_of_gyration(single)
            (confdir / f"conf_{i:02d}.xyz").write_text(
                to_xyz(single, f"conf {i} ff={ff} E={e:.3f} Rg={rg:.2f} smiles={smi}"))
        print(f"  {g['id']:16s} {len(conf_ids)} conformers (ff={ff}, q={net_q:+d})")
        for label, charge, mult, n_e in g["states"]:
            rows.append(dict(id=g["id"], name=g["name"], family=g["family"],
                             state=label, charge=charge, mult_hint=mult, n_e=n_e,
                             n_conf=len(conf_ids), solv_preopt=int(charge != 0), smiles=smi))
            print(f"      {label:5s} q={charge:+d} mult_hint={mult} n_e={n_e:+d}")
    manifest = LIBRARY / "manifest.csv"
    with manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "name", "family", "state", "charge",
                                          "mult_hint", "n_e", "n_conf", "solv_preopt", "smiles"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} states across {len(groups)} groups -> {manifest}")


if __name__ == "__main__":
    sys.exit(main())
