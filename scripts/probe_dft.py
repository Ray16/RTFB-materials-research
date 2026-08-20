"""Verify the composite-protocol level of theory actually works in this env BEFORE a
production submit: gpu4pyscf import, r2SCAN-D4 and wB97X-D3 functionals, D3/D4 dispersion,
diffuse basis on an anion, and one SMD geometry step. Prints a PASS/FAIL per capability.

  CUDA_VISIBLE_DEVICES=1 python scripts/probe_dft.py --backend gpu
  python scripts/probe_dft.py --backend cpu     # fallback check
"""
import argparse
import sys


def dft_module(backend):
    if backend == "gpu":
        from gpu4pyscf import dft
        return dft
    from pyscf import dft
    return dft


def try_sp(backend, xc, basis, charge, spin, disp=None, smd=True, label=""):
    from pyscf import gto
    dft = dft_module(backend)
    # tiny test system: formate anion (has an anion + is small) or water
    atom = "O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467" if charge == 0 else \
           "C 0 0 0; O 0 0 1.25; O 1.1 0 -0.6; H -1.0 0 -0.4"
    mol = gto.M(atom=atom, basis=basis, charge=charge, spin=spin, verbose=0)
    KS = dft.RKS if spin == 0 else dft.UKS
    mf = KS(mol); mf.xc = xc
    if disp:
        mf.disp = disp
    if smd:
        mf = mf.SMD(); mf.with_solvent.solvent = "acetonitrile"
    e = float(mf.kernel())
    print(f"  PASS  {label:34s} E={e:.6f} Ha  conv={bool(mf.converged)}")
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    args = ap.parse_args()

    print(f"== probe DFT capabilities (backend={args.backend}) ==")
    try:
        dft_module(args.backend)
        print(f"  PASS  import {args.backend} dft module")
    except Exception as e:
        print(f"  FAIL  import {args.backend} dft module: {e}")
        sys.exit(1)

    checks = [
        ("r2scan", "def2-svp", 0, 0, "d4",  "r2SCAN-D4 / def2-svp  (opt, neutral)"),
        ("r2scan", "def2-svpd", -1, 1, "d4", "r2SCAN-D4 / def2-svpd (opt, anion)"),
        ("wb97x-d3", "def2-tzvp", 0, 0, None, "wB97X-D3 / def2-tzvp  (energy, neutral)"),
        ("wb97x-d3", "def2-tzvpd", -1, 1, None, "wB97X-D3 / def2-tzvpd (energy, anion)"),
    ]
    npass = 0
    for xc, basis, q, s, disp, label in checks:
        try:
            try_sp(args.backend, xc, basis, q, s, disp=disp, label=label)
            npass += 1
        except Exception as e:
            print(f"  FAIL  {label:34s} {type(e).__name__}: {str(e)[:90]}")

    # one SMD geometry-opt step (proves gradients work end to end)
    try:
        from pyscf import gto
        dft = dft_module(args.backend)
        if args.backend == "gpu":
            from gpu4pyscf.geomopt.geometric_solver import optimize
        else:
            from pyscf.geomopt.geometric_solver import optimize
        mol = gto.M(atom="O 0 0 0.117; H 0 0.757 -0.467; H 0 -0.757 -0.467",
                    basis="def2-svp", verbose=0)
        mf = dft.RKS(mol); mf.xc = "r2scan"; mf.disp = "d4"
        mf = mf.SMD(); mf.with_solvent.solvent = "acetonitrile"
        optimize(mf, maxsteps=2)
        print("  PASS  SMD geometry-opt step (r2SCAN-D4)")
        npass += 1
    except Exception as e:
        print(f"  FAIL  SMD geometry-opt step: {type(e).__name__}: {str(e)[:90]}")

    print(f"== {npass}/{len(checks)+1} capabilities OK ==")
    sys.exit(0 if npass == len(checks) + 1 else 2)


if __name__ == "__main__":
    main()
