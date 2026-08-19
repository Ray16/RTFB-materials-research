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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", required=True)
    ap.add_argument("--charge", type=int, required=True)
    ap.add_argument("--mult", type=int, required=True)
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--out", default=None, help="write result JSON here")
    ap.add_argument("--nthreads", type=int, default=0, help="OMP threads (0=leave as is)")
    args = ap.parse_args()

    if args.nthreads:
        import pyscf.lib
        pyscf.lib.num_threads(args.nthreads)

    res = dft_smd(Path(args.xyz), args.charge, args.mult, args.xc, args.basis)
    print(json.dumps(res, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
