"""Build the starting-candidate set (config/starting_candidates.py) into the SAME library +
manifest as the decorated monomers, so they run the identical UMA -> DFT+SMD -> redox pipeline.

Mirrors build_validation.py's MERGE behaviour: existing manifest rows for ids we are NOT
(re)building are preserved; freshly built candidate rows are appended. This never clobbers
the validation/ferrocene rows the way the from-scratch build.py would.

Reuses the grafting (decorate) + conformer ensemble from build.py, so candidates are built
byte-for-byte the same way as redox_groups.py monomers.

  python -m redox.build.candidates              # build all starting candidates
  python -m redox.build.candidates --only pmdi  # one candidate
  python -m redox.build.candidates --config config/merrifield_multielectron.py
"""
from __future__ import annotations
import argparse
import importlib.util
from pathlib import Path

from rdkit import Chem
from redox.build.contracts import validate_group
from redox.build.manifest import COLS, locked_manifest, read_rows, write_rows

from redox.build.structures import (decorate, conformer_ensemble, to_xyz, unassigned_stereo,
                         radius_of_gyration)

ROOT = Path(__file__).resolve().parents[3]
CONFIG = ROOT / "config" / "starting_candidates.py"   # default; override with --config
LIBRARY = ROOT / "library"
MANIFEST = LIBRARY / "manifest.csv"

def _load_config(path=None):
    path = Path(path) if path else CONFIG
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.SCAFFOLD, mod.GROUPS


def read_manifest_rows():
    return read_rows(MANIFEST)


def build_group(scaffold: str, g: dict) -> tuple[str, int]:
    """Graft the candidate onto the scaffold, embed a conformer ensemble, write xyz files.
    Returns (decorated_smiles, n_conf). Raises on unassigned stereo (hard guard)."""
    gid = g["id"]
    mol = validate_group(scaffold, g)
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
    global LIBRARY, MANIFEST
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="comma-separated candidate IDs")
    ap.add_argument("--config", default=None,
                    help="candidate config .py to build (default: config/starting_candidates.py)")
    ap.add_argument("--library", type=Path, default=LIBRARY, help="output library (for isolated smoke tests)")
    ap.add_argument("--dry-run", action="store_true", help="validate all definitions without writing")
    args = ap.parse_args()
    LIBRARY = args.library.resolve()
    MANIFEST = LIBRARY / "manifest.csv"

    scaffold, groups = _load_config(args.config)
    if args.only is not None:
        selected = {value.strip() for value in args.only.split(",") if value.strip()}
        unknown = selected - {group["id"] for group in groups}
        if not selected or unknown:
            raise SystemExit(f"invalid --only selection: {sorted(unknown) if unknown else 'empty'}")
        groups = [group for group in groups if group["id"] in selected]

    if args.dry_run:
        _build_selected(scaffold, groups, dry_run=True)
        return
    # Different motif writers share one manifest. Hold its lock across validation,
    # building and the final merge so a second process cannot lose the first's rows
    # or create an unregistered directory for the same ID.
    LIBRARY.mkdir(parents=True, exist_ok=True)
    with locked_manifest(MANIFEST):
        _build_selected(scaffold, groups, dry_run=False)


def _build_selected(scaffold, groups, dry_run=False):
    # Re-read only after acquiring the global writer lock in main().
    existing = read_manifest_rows()
    if len({g["id"] for g in groups}) != len(groups):
        raise ValueError("duplicate candidate IDs in config")
    # Validate the whole batch before writing anything; never reuse an ID for changed chemistry.
    pending = []
    for g in groups:
        mol = validate_group(scaffold, g)
        smi = Chem.MolToSmiles(mol)
        old = [r for r in existing if r["id"] == g["id"]]
        if old:
            wanted = {(s, q, m, ne) for s, q, m, ne in g["states"]}
            found = {(r["state"], int(r["charge"]), int(r["mult_hint"]), int(r["n_e"])) for r in old}
            if wanted != found or any(Chem.MolToSmiles(Chem.MolFromSmiles(r["smiles"])) != smi for r in old):
                raise ValueError(f"{g['id']}: ID already belongs to different chemistry/states; use a new ID")
            confs = list((LIBRARY / g["id"] / "conformers").glob("conf_*.xyz"))
            if not confs:
                raise ValueError(f"{g['id']}: existing library is incomplete; inspect before rebuilding")
            print(f"[skip] {g['id']}: matching existing library")
        else:
            if (LIBRARY / g["id"]).exists():
                raise ValueError(f"{g['id']}: directory exists without a manifest entry; inspect before rebuilding")
            pending.append(g)
        print(f"[valid] {g['id']}: {smi}")
    if dry_run or not pending:
        return
    groups = pending
    LIBRARY.mkdir(parents=True, exist_ok=True)
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
    write_rows(MANIFEST, all_rows)
    print(f"\n[manifest] {len(new_rows)} candidate states written; "
          f"manifest now {len(all_rows)} rows -> {MANIFEST}")


if __name__ == "__main__":
    main()
