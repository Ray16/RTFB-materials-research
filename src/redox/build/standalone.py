"""Build the standalone (un-tethered) candidate forms (config/standalone.py) into the shared
library + manifest, so they run the identical UMA -> DFT+SMD -> reorg pipeline. Direct SMILES
(no scaffold decoration). MERGE-safe like build_validation.py.

  python -m redox.build.standalone
  python -m redox.build.standalone --only aq_benzyloxy_sa,nq_benzyloxy_sa
"""
from __future__ import annotations
import csv
import importlib.util
from pathlib import Path

from rdkit import Chem

from redox.build.structures import conformer_ensemble, to_xyz, unassigned_stereo, radius_of_gyration

ROOT = Path(__file__).resolve().parents[3]
CONFIG = ROOT / "config" / "standalone.py"
LIBRARY = ROOT / "library"
MANIFEST = LIBRARY / "manifest.csv"
COLS = ["id", "name", "family", "state", "charge", "mult_hint", "n_e", "n_conf",
        "solv_preopt", "smiles"]


def _load():
    spec = importlib.util.spec_from_file_location("standalone", CONFIG)
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m.STANDALONE


def _rows():
    if not MANIFEST.exists():
        return []
    with MANIFEST.open() as f:
        return list(csv.DictReader(f))


def build_one(mol_def):
    gid = mol_def["id"]
    mol = Chem.MolFromSmiles(mol_def["smiles"])
    if mol is None:
        raise ValueError(f"{gid}: bad SMILES {mol_def['smiles']!r}")
    smi = Chem.MolToSmiles(mol)
    n_unspec = unassigned_stereo(mol)
    if n_unspec:
        raise ValueError(f"{gid}: {n_unspec} unassigned stereo element(s) — fix before building")
    (LIBRARY / gid).mkdir(parents=True, exist_ok=True)
    (LIBRARY / gid / f"{gid}.smiles").write_text(smi + "\n")
    net_q = Chem.GetFormalCharge(mol)
    confdir = LIBRARY / gid / "conformers"; confdir.mkdir(exist_ok=True)
    molH, conf_ids, ff, energies = conformer_ensemble(mol, net_charge=net_q)
    for f in confdir.glob("conf_*.xyz"):
        f.unlink()
    for i, (cid, e) in enumerate(zip(conf_ids, energies)):
        single = Chem.Mol(molH); single.RemoveAllConformers()
        single.AddConformer(molH.GetConformer(cid), assignId=True)
        (confdir / f"conf_{i:02d}.xyz").write_text(
            to_xyz(single, f"conf {i} ff={ff} E={e:.3f} Rg={radius_of_gyration(single):.2f} smiles={smi}"))
    print(f"  {gid:16s} {len(conf_ids)} conformers (q={net_q:+d})")
    return smi, len(conf_ids)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None,
                    help="comma-separated standalone ids to (re)build; default: all. "
                         "Use this to add new entries without regenerating the conformer "
                         "ensembles of standalones that are already computed.")
    args = ap.parse_args()

    mols = _load()
    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        mols = [m for m in mols if m["id"] in want]
        missing = want - {m["id"] for m in mols}
        if missing:
            raise SystemExit(f"no standalone entry for: {sorted(missing)}")
    existing = _rows()
    new_rows = []
    for md in mols:
        smi, n_conf = build_one(md)
        for (label, charge, mult, n_e) in md["states"]:
            new_rows.append(dict(id=md["id"], name=md["name"], family=md["family"],
                                 state=label, charge=charge, mult_hint=mult, n_e=n_e,
                                 n_conf=n_conf, solv_preopt=int(charge != 0), smiles=smi))
            print(f"      {label:5s} q={charge:+d} mult_hint={mult} n_e={n_e:+d}")
    rebuilt = {md["id"] for md in mols}
    kept = [r for r in existing if r["id"] not in rebuilt]
    with MANIFEST.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
        for r in kept + new_rows:
            w.writerow({c: r.get(c, "") for c in COLS})
    print(f"\n[manifest] {len(new_rows)} standalone states written; now {len(kept)+len(new_rows)} rows")


if __name__ == "__main__":
    main()
