"""Build the §V validation set (parent redox cores) into the SAME library + manifest as the
decorated monomers, so they run through the identical UMA -> DFT+SMD -> redox pipeline.

Most cores embed with RDKit like the monomers. Ferrocene is a metallocene: RDKit cannot
embed the sandwich, so we synthesize an eclipsed-D5h starting geometry analytically and let
UMA (charge/spin-aware) + DFT+SMD relax it. Ferrocene is also the internal reference — its
computed Fc/Fc+ absolute potential (redox.py) becomes electrolyte.FC_ABS_COMPUTED_V, and is
itself a validation check (should land near the OROP MeCN B3LYP value ~4.63 V).

  python -m redox.build_validation            # build all validation cores
  python -m redox.build_validation --only ferrocene
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import math
from pathlib import Path

from rdkit import Chem

from redox.build import conformer_ensemble, to_xyz, unassigned_stereo

ROOT = Path(__file__).resolve().parents[2]
LIBRARY = ROOT / "library"
MANIFEST = LIBRARY / "manifest.csv"


def _load_validation():
    spec = importlib.util.spec_from_file_location(
        "validation", ROOT / "config" / "validation.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.VALIDATION


def ferrocene_xyz() -> str:
    """Eclipsed-D5h ferrocene starting geometry (Fe + two Cp rings).

    Standard metrics: ring C-C 1.43 A -> C radius 1.216 A; ring centroid-Fe 1.65 A
    (=> Fe-C 2.05 A); C-H 1.08 A radially outward, in-plane. UMA/DFT relax from here.
    """
    r_c = 1.43 / (2.0 * math.sin(math.pi / 5))   # 1.216 A
    r_h = r_c + 1.08                              # H radially outward, in ring plane
    z = 1.65                                      # ring centroid to Fe along axis
    atoms = [("Fe", 0.0, 0.0, 0.0)]
    for sign in (+1, -1):                         # top and bottom rings, eclipsed
        for k in range(5):
            ang = 2.0 * math.pi * k / 5.0
            cx, cy = math.cos(ang), math.sin(ang)
            atoms.append(("C", r_c * cx, r_c * cy, sign * z))
            atoms.append(("H", r_h * cx, r_h * cy, sign * z))
    lines = [str(len(atoms)), "ferrocene eclipsed-D5h starting geometry"]
    for s, x, y, zz in atoms:
        lines.append(f"{s:2s} {x:14.8f} {y:14.8f} {zz:14.8f}")
    return "\n".join(lines) + "\n"


def build_core(core: dict) -> tuple[list[str], int]:
    """Write conformer xyz(s) for one validation core. Returns (xyz_texts, n_conf)."""
    gid = core["id"]
    confdir = LIBRARY / gid / "conformers"
    confdir.mkdir(parents=True, exist_ok=True)

    if core.get("special_geometry"):
        (confdir / "conf_00.xyz").write_text(ferrocene_xyz())
        return (["conf_00.xyz"], 1)

    mol = Chem.MolFromSmiles(core["smiles"])
    if mol is None:
        raise ValueError(f"{gid}: bad SMILES {core['smiles']!r}")
    n_unspec = unassigned_stereo(mol)
    if n_unspec:
        raise ValueError(f"{gid}: {n_unspec} unassigned stereocenter(s) — fix before building")
    net_q = Chem.GetFormalCharge(mol)
    molH, conf_ids, ff, energies = conformer_ensemble(mol, net_charge=net_q)
    names = []
    for i, cid in enumerate(conf_ids):
        # write each conformer by selecting it in to_xyz via a per-conf mol
        txt = _conf_to_xyz(molH, cid, f"{gid} conf {i} ff={ff} E={energies[i]:.3f}")
        name = f"conf_{i:02d}.xyz"
        (confdir / name).write_text(txt)
        names.append(name)
    return (names, len(names))


def _conf_to_xyz(molH, conf_id, comment):
    conf = molH.GetConformer(conf_id)
    lines = [str(molH.GetNumAtoms()), comment]
    for atom in molH.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {p.x:14.8f} {p.y:14.8f} {p.z:14.8f}")
    return "\n".join(lines) + "\n"


def read_manifest_rows():
    if not MANIFEST.exists():
        return []
    with MANIFEST.open() as f:
        return list(csv.DictReader(f))


COLS = ["id", "name", "family", "state", "charge", "mult_hint",
        "n_e", "n_conf", "solv_preopt", "smiles"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="restrict to one core id")
    args = ap.parse_args()

    cores = _load_validation()
    if args.only:
        cores = [c for c in cores if c["id"] == args.only]
        if not cores:
            raise SystemExit(f"no validation core matches --only {args.only!r}")

    existing = read_manifest_rows()
    have_ids = {r["id"] for r in existing}
    new_rows = []
    for core in cores:
        gid = core["id"]
        names, n_conf = build_core(core)
        print(f"[built] {gid}: {n_conf} conformer(s)")
        for (label, charge, mult, n_e) in core["states"]:
            new_rows.append({
                "id": gid, "name": core["name"], "family": "validation",
                "state": label, "charge": charge, "mult_hint": mult,
                "n_e": n_e, "n_conf": n_conf, "solv_preopt": "",
                "smiles": core["smiles"],
            })

    # merge: keep all existing rows that are NOT for the cores we just (re)built,
    # then append the freshly built validation rows.
    rebuilt_ids = {c["id"] for c in cores}
    kept = [r for r in existing if r["id"] not in rebuilt_ids]
    all_rows = kept + new_rows
    with MANIFEST.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader()
        for r in all_rows:
            w.writerow({c: r.get(c, "") for c in COLS})
    print(f"[manifest] {len(new_rows)} validation states written; "
          f"manifest now {len(all_rows)} rows -> {MANIFEST}")


if __name__ == "__main__":
    main()
