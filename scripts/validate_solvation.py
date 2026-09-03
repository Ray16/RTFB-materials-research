#!/usr/bin/env python
"""Validate the pipeline's computed ΔG_solv against EXPERIMENT (Ray's Sept-2 action item).

We report ΔG_solv = e_smd - e_gas (full SMD solvation free energy — gpu4pyscf's SMD includes
the non-electrostatic CDS term). This checks it is physically reasonable by comparing to
measured solvation free energies from multiple experimental databases:

  FreeSolv  (WATER, 642 exp hydration free energies)  -> machinery check, gold standard, open.
  MNSol     (ACETONITRILE + others; the set SMD was parameterized against)  -> on-target,
            gold standard. MeCN has 7 neutrals + 69 ions (39 cations, 30 anions).

Protocol matches SMD's definition: ΔG_solv = e_smd - e_gas at a FIXED geometry, both single
points at ωB97M-V/def2-TZVP(+diffuse for ions)+VV10. FreeSolv rows are RDKit-embedded then
optimized; MNSol rows use MNSol's own M06-2X/MG3S gas geometry as-is (do_opt=False) — the
standard SMD solvation protocol. Resumable (per-mol result cache). Every row carries an
explicit charge + spin multiplicity; ion spin comes from MNSol's geometry header.

  python scripts/validate_solvation.py --build-freesolv 24
  python scripts/validate_solvation.py --build-mnsol acetonitrile
  CUDA_VISIBLE_DEVICES=3 python scripts/validate_solvation.py --run <refset.csv> --shard 4:0
  python scripts/validate_solvation.py --report
"""
from __future__ import annotations
import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
REFDIR = ROOT / "data" / "raw" / "validation"
CALC = ROOT / "calcs" / "solvation"
RESULTS = ROOT / "results"
EV2KCAL = 23.060541945329334
MNSOL = REFDIR / "MNSol" / "extracted" / "MNSolDatabase_v2012"

REFCOLS = ["id", "source", "solvent", "charge", "mult", "exp_kcal", "geom", "smiles", "name"]


# ---------------------------------------------------------------- reference sets
def build_freesolv(n: int):
    """Stratified sample of n small neutral molecules spanning FreeSolv's exp range (water)."""
    from rdkit import Chem
    rows = []
    for line in (REFDIR / "FreeSolv" / "database.txt").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = [x.strip() for x in line.split(";")]
        cid, smi, name, exp = f[0], f[1], f[2], f[3]
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        if m.GetNumHeavyAtoms() > 10 or not {a.GetSymbol() for a in m.GetAtoms()} <= {
                "C", "N", "O", "F", "S", "Cl", "H", "P", "Br"}:
            continue
        rows.append(dict(id=cid, source="FreeSolv", solvent="water", charge=0, mult=1,
                         exp_kcal=float(exp), geom="", smiles=smi, name=name))
    rows.sort(key=lambda r: r["exp_kcal"])
    idx = sorted(set(round(i * (len(rows) - 1) / (n - 1)) for i in range(n)))
    _write_refset("solvation_refset_water.csv", [rows[i] for i in idx])


def _mnsol_geom(handle: str):
    """Parse an MNSol all_solutes/<handle>.xyz (Z x y z; header line has 'charge mult').
    Returns (charge, mult, xyz_text_with_symbols)."""
    from ase.data import chemical_symbols
    lines = (MNSOL / "all_solutes" / f"{handle}.xyz").read_text().splitlines()
    charge = mult = None
    atoms = []
    for ln in lines:
        t = ln.split()
        if len(t) == 2:
            try:
                charge, mult = int(t[0]), int(t[1]); continue
            except ValueError:
                pass
        if len(t) == 4:
            try:
                z = int(float(t[0])); sym = chemical_symbols[z]
                atoms.append(f"{sym} {t[1]} {t[2]} {t[3]}")
            except (ValueError, IndexError):
                pass
    xyz = f"{len(atoms)}\nMNSol {handle} q={charge} m={mult}\n" + "\n".join(atoms) + "\n"
    return charge, mult, xyz


def build_mnsol(solvent: str):
    """All MNSol rows for a solvent (type=abs), using MNSol's own gas geometries."""
    tbl = MNSOL / "MNSol_alldata.txt"
    out = []
    for r in csv.DictReader(tbl.open(), delimiter="\t"):
        if r["Solvent"].lower() != solvent.lower() or r["type"] != "abs":
            continue
        h = r["FileHandle"]
        if not (MNSOL / "all_solutes" / f"{h}.xyz").exists():
            continue
        charge, mult, _ = _mnsol_geom(h)
        out.append(dict(id=h, source="MNSol", solvent=solvent,
                        charge=charge if charge is not None else r["Charge"],
                        mult=mult if mult is not None else 1,
                        exp_kcal=float(r["DeltaGsolv"]), geom=f"mnsol:{h}",
                        smiles="", name=r["SoluteName"]))
    _write_refset(f"solvation_refset_{solvent.lower()}.csv", out)


def _write_refset(fname, rows):
    out = REFDIR / fname
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=REFCOLS); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in REFCOLS})
    exps = [float(r["exp_kcal"]) for r in rows]
    print(f"wrote {out}  ({len(rows)} molecules, exp {min(exps):.2f} .. {max(exps):.2f} kcal/mol)")


# ---------------------------------------------------------------- run one molecule
def _embed_xyz(smiles: str, out_xyz: Path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, randomSeed=0xf00d) != 0:
        AllChem.EmbedMolecule(m, useRandomCoords=True, randomSeed=7)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    lines = [str(m.GetNumAtoms()), smiles]
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"{a.GetSymbol():2s} {p.x:.6f} {p.y:.6f} {p.z:.6f}")
    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    out_xyz.write_text("\n".join(lines) + "\n")


def run(refset: str, shard: str | None, backend: str, force: bool):
    from redox.dft import dft_smd
    rows = list(csv.DictReader((REFDIR / refset).open()))
    if shard:
        nsh, ish = (int(x) for x in shard.split(":"))
        rows = [r for k, r in enumerate(rows) if k % nsh == ish]
    for r in rows:
        cid, solv = r["id"], r["solvent"]
        charge, mult = int(r["charge"]), int(r["mult"])
        d = CALC / solv / cid
        rj = d / "result.json"
        if rj.exists() and not force:
            print(f"[skip] {cid}", flush=True); continue
        d.mkdir(parents=True, exist_ok=True)
        xyz = d / "seed.xyz"
        try:
            if r["geom"].startswith("mnsol:"):
                _, _, xyztext = _mnsol_geom(r["geom"].split(":", 1)[1])
                xyz.write_text(xyztext)
                do_opt = False          # use MNSol's gas geometry as-is (SMD protocol)
            else:
                _embed_xyz(r["smiles"], xyz)
                do_opt = True           # relax RDKit seed
            res = dft_smd(xyz, charge=charge, mult=mult, solvent=solv,
                          do_opt=do_opt, do_freq=False, do_gas=True, do_smd=True,
                          opt_out=(d / "opt.xyz") if do_opt else None, backend=backend)
        except Exception as e:
            print(f"[fail] {cid}: {type(e).__name__}: {str(e)[:100]}", flush=True); continue
        dg = res.get("dG_solv_eV")
        res.update(id=cid, source=r["source"], solvent=solv, charge=charge, mult=mult,
                   exp_kcal=r["exp_kcal"], name=r.get("name", ""),
                   dG_solv_kcal=(dg * EV2KCAL) if dg is not None else None)
        rj.write_text(json.dumps(res, indent=2))
        got = res["dG_solv_kcal"]
        print(f"[done] {cid} q={charge} calc={got:+.2f} exp={float(r['exp_kcal']):+.2f} "
              f"kcal/mol ({solv}/{r['source']})" if got is not None
              else f"[warn] {cid} no dG_solv", flush=True)


# ---------------------------------------------------------------- report
def _charge_class(q):
    return "neutral" if int(q) == 0 else ("cation" if int(q) > 0 else "anion")


def report():
    import numpy as np
    rows = []
    for rj in sorted(CALC.glob("*/*/result.json")):
        d = json.loads(rj.read_text())
        if d.get("dG_solv_kcal") is None:
            continue
        q = int(d.get("charge", 0))
        rows.append(dict(solvent=d["solvent"], source=d.get("source", "?"),
                         charge_class=_charge_class(q), id=d["id"], name=d.get("name", ""),
                         charge=q, calc_kcal=round(d["dG_solv_kcal"], 2),
                         exp_kcal=round(float(d["exp_kcal"]), 2),
                         err=round(d["dG_solv_kcal"] - float(d["exp_kcal"]), 2)))
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "solvation_validation.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["solvent", "source", "charge_class", "charge",
                                           "id", "name", "calc_kcal", "exp_kcal", "err"])
        w.writeheader()
        w.writerows(sorted(rows, key=lambda x: (x["solvent"], x["charge_class"], x["exp_kcal"])))
    print(f"wrote {out}  ({len(rows)} molecules)\n")

    def stats(sub, label):
        e = np.array([r["err"] for r in sub]); c = np.array([r["calc_kcal"] for r in sub])
        x = np.array([r["exp_kcal"] for r in sub])
        r2 = np.corrcoef(c, x)[0, 1] ** 2 if len(e) > 1 else float("nan")
        print(f"  {label:34s} n={len(e):3d}  MAE={np.abs(e).mean():5.2f}  "
              f"RMSE={np.sqrt((e**2).mean()):5.2f}  bias={e.mean():+5.2f}  R²={r2:5.3f}")

    print("=== by solvent × charge-class (kcal/mol) ===")
    for solv in sorted({r["solvent"] for r in rows}):
        for cc in ["neutral", "cation", "anion"]:
            sub = [r for r in rows if r["solvent"] == solv and r["charge_class"] == cc]
            if sub:
                src = sorted({r["source"] for r in sub})[0]
                stats(sub, f"{solv}/{cc} [{src}]")
    print("=== overall neutrals (the screening-relevant property) ===")
    stats([r for r in rows if r["charge_class"] == "neutral"], "ALL neutrals")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-freesolv", type=int, metavar="N")
    ap.add_argument("--build-mnsol", metavar="SOLVENT")
    ap.add_argument("--run", metavar="REFSET.csv")
    ap.add_argument("--shard", default=None, help="'n:i'")
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.build_freesolv:
        build_freesolv(a.build_freesolv)
    if a.build_mnsol:
        build_mnsol(a.build_mnsol)
    if a.run:
        run(a.run, a.shard, a.backend, a.force)
    if a.report:
        report()


if __name__ == "__main__":
    main()
