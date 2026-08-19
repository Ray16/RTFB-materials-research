"""DFT + SMD(acetonitrile) single-point energies (the solvation step).

Takes a geometry (e.g. a UMA-relaxed xyz) plus explicit charge + spin multiplicity and
returns the SMD-solvated DFT energy (and the gas-phase energy for the solvation term).
This supplies the solution-phase energies the redox-potential calculation needs.

Level of theory is PROVISIONAL (see docs/PLAN.md §9 [DECIDE] and the validation gate §V):
default B3LYP/def2-SVP is a fast first pass; anions want diffuse functions (def2-SVPD)
and redox generally wants a dispersion-corrected range-separated hybrid (wB97X-D).

  python -m redox.dft --xyz calcs/uma/pyridinium/ox/relaxed.xyz --charge 1 --mult 1
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HARTREE_EV = 27.211386245988


def _electrolyte():
    spec = importlib.util.spec_from_file_location(
        "electrolyte", ROOT / "config" / "electrolyte.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dft_smd(xyz_path: Path, charge: int, mult: int,
            xc: str = "b3lyp", basis: str = "def2-svp",
            solvent: str | None = None, do_gas: bool = True) -> dict:
    """Return DFT energies (Ha) in SMD solvent and (optionally) gas phase."""
    from pyscf import gto, dft
    from ase.io import read

    atoms = read(str(xyz_path))
    atom_str = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                         for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    nunpaired = int(mult) - 1  # pyscf mol.spin = 2S = (mult-1)
    mol = gto.M(atom=atom_str, basis=basis, charge=int(charge), spin=nunpaired,
                verbose=0)

    KS = dft.RKS if nunpaired == 0 else dft.UKS
    solvent = solvent or _electrolyte().SOLVENT["name"]

    # Solvated
    mfs = KS(mol); mfs.xc = xc
    mfs = mfs.SMD()
    mfs.with_solvent.solvent = solvent
    e_smd = float(mfs.kernel())
    conv_smd = bool(mfs.converged)

    out = dict(xyz=str(xyz_path), charge=int(charge), mult=int(mult),
               xc=xc, basis=basis, solvent=solvent,
               e_smd_Ha=e_smd, e_smd_eV=e_smd * HARTREE_EV, converged_smd=conv_smd)

    if do_gas:
        mfg = KS(mol); mfg.xc = xc
        e_gas = float(mfg.kernel())
        out.update(e_gas_Ha=e_gas, e_gas_eV=e_gas * HARTREE_EV,
                   converged_gas=bool(mfg.converged),
                   dG_solv_eV=(e_smd - e_gas) * HARTREE_EV)
    return out


def read_manifest():
    import csv
    with (ROOT / "library" / "manifest.csv").open() as f:
        return list(csv.DictReader(f))


def geom_for(gid: str, state: str) -> Path:
    """Prefer the UMA-relaxed geometry; fall back to the conformer seed."""
    relaxed = ROOT / "calcs" / "uma" / gid / state / "relaxed.xyz"
    return relaxed if relaxed.exists() else ROOT / "library" / gid / f"{state}.xyz"


def run_batch(rows, xc, basis, force):
    outroot = ROOT / "calcs" / "dft"
    for r in rows:
        gid, st = r["id"], r["state"]
        outdir = outroot / gid / st
        rj = outdir / "result.json"
        if rj.exists() and not force:
            print(f"[skip] {gid}/{st}", flush=True); continue
        geom = geom_for(gid, st)
        if not geom.exists():
            print(f"[miss] {gid}/{st} no geometry ({geom})", flush=True); continue
        print(f"[run ] {gid}/{st} q={r['charge']} m={r['mult']} ...", flush=True)
        try:
            res = dft_smd(geom, int(r["charge"]), int(r["mult"]), xc, basis)
        except Exception as e:
            print(f"[fail] {gid}/{st}: {type(e).__name__}: {str(e)[:120]}", flush=True)
            continue
        res.update(id=gid, state=st, n_e=int(r["n_e"]))
        outdir.mkdir(parents=True, exist_ok=True)
        rj.write_text(json.dumps(res, indent=2))
        print(f"[done] {gid}/{st} E_smd={res['e_smd_eV']:.3f} eV "
              f"dGsolv={res.get('dG_solv_eV', float('nan')):.3f} eV "
              f"conv={res['converged_smd']}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", default=None, help="single-geometry mode")
    ap.add_argument("--charge", type=int)
    ap.add_argument("--mult", type=int)
    ap.add_argument("--all", action="store_true", help="batch over manifest states")
    ap.add_argument("--only", default=None, help="'group' or 'group:state'")
    ap.add_argument("--shard", default=None, help="'n:i' for CPU fan-out")
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--out", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--nthreads", type=int, default=0)
    args = ap.parse_args()

    if args.nthreads:
        import pyscf.lib
        pyscf.lib.num_threads(args.nthreads)

    if args.all or args.only or args.shard:
        rows = read_manifest()
        if args.only:
            gid, _, st = args.only.partition(":")
            rows = [r for r in rows if r["id"] == gid and (not st or r["state"] == st)]
        if args.shard:
            n, i = (int(x) for x in args.shard.split(":"))
            rows = [r for k, r in enumerate(rows) if k % n == i]
        run_batch(rows, args.xc, args.basis, args.force)
        return

    res = dft_smd(Path(args.xyz), args.charge, args.mult, args.xc, args.basis)
    print(json.dumps(res, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
