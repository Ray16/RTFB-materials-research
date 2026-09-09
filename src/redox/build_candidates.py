"""Build the starting-candidate set (config/starting_candidates.py) into the SAME library +
manifest as the decorated monomers, so they run the identical UMA -> DFT+SMD -> redox pipeline.

Mirrors build_validation.py's MERGE behaviour: existing manifest rows for ids we are NOT
(re)building are preserved; freshly built candidate rows are appended. This never clobbers
the validation/ferrocene rows the way the from-scratch build.py would.

Reuses the grafting (decorate) + conformer ensemble from build.py, so candidates are built
byte-for-byte the same way as redox_groups.py monomers.

  python -m redox.build_candidates              # build all starting candidates
  python -m redox.build_candidates --only pmdi  # one candidate
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
from pathlib import Path

from rdkit import Chem

from redox.build import (decorate, conformer_ensemble, to_xyz, unassigned_stereo,
                         radius_of_gyration)

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "config" / "starting_candidates.py"
LIBRARY = ROOT / "library"
MANIFEST = LIBRARY / "manifest.csv"

COLS = ["id", "name", "family", "state", "charge", "mult_hint",
        "n_e", "n_conf", "solv_preopt", "smiles"]


def _load_config():
    spec = importlib.util.spec_from_file_location("starting_candidates", CONFIG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.SCAFFOLD, mod.GROUPS


def read_manifest_rows():
    if not MANIFEST.exists():
        return []
    with MANIFEST.open() as f:
        return list(csv.DictReader(f))


def build_group(scaffold: str, g: dict) -> tuple[str, int]:
    """Graft the candidate onto the scaffold, embed a conformer ensemble, write xyz files.
    Returns (decorated_smiles, n_conf). Raises on unassigned stereo (hard guard)."""
    gid = g["id"]
    mol = decorate(scaffold, g["frag"])
    smi = Chem.MolToSmiles(mol)
    n_unspec = unassigned_stereo(mol)
    if n_unspec:
        raise ValueError(f"{gid}: {n_unspec} unassigned stereo element(s) — fix SMILES "
                         f"or enumerate before building. SMILES={smi}")
    (LIBRARY / gid).mkdir(parents=True, exist_ok=True)
    (LIBRARY / gid / f"{gid}.smiles").write_text(smi + "\n")
    net_q = Chem.GetFormalCharge(mol)
    confdir = LIBRARY / gid / "conformers"
    confdir.mkdir(exist_ok=True)
    molH, conf_ids, ff, energies = conformer_ensemble(mol, net_charge=net_q)
    for f in confdir.glob("conf_*.xyz"):
        f.unlink()  # clear stale ensemble
    for i, (cid, e) in enumerate(zip(conf_ids, energies)):
        single = Chem.Mol(molH); single.RemoveAllConformers()
        single.AddConformer(molH.GetConformer(cid), assignId=True)
        rg = radius_of_gyration(single)
        (confdir / f"conf_{i:02d}.xyz").write_text(
            to_xyz(single, f"conf {i} ff={ff} E={e:.3f} Rg={rg:.2f} smiles={smi}"))
    print(f"  {gid:16s} {len(conf_ids)} conformers (ff={ff}, q={net_q:+d})")
    return smi, len(conf_ids)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="restrict to one candidate id")
    args = ap.parse_args()

    scaffold, groups = _load_config()
    if args.only:
        groups = [g for g in groups if g["id"] == args.only]
        if not groups:
            raise SystemExit(f"no candidate matches --only {args.only!r}")

    existing = read_manifest_rows()
    new_rows = []
    for g in groups:
        smi, n_conf = build_group(scaffold, g)
        for (label, charge, mult, n_e) in g["states"]:
            new_rows.append(dict(id=g["id"], name=g["name"], family=g["family"],
                                 state=label, charge=charge, mult_hint=mult, n_e=n_e,
                                 n_conf=n_conf, solv_preopt=int(charge != 0), smiles=smi))
            print(f"      {label:5s} q={charge:+d} mult_hint={mult} n_e={n_e:+d}")

    rebuilt_ids = {g["id"] for g in groups}
    kept = [r for r in existing if r["id"] not in rebuilt_ids]
    all_rows = kept + new_rows
    with MANIFEST.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in all_rows:
            w.writerow({c: r.get(c, "") for c in COLS})
    print(f"\n[manifest] {len(new_rows)} candidate states written; "
          f"manifest now {len(all_rows)} rows -> {MANIFEST}")


if __name__ == "__main__":
    main()
