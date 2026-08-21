"""DFT + SMD(acetonitrile) energies (the solvation step).

Takes a geometry (e.g. a UMA-relaxed xyz) plus explicit charge + spin multiplicity and
returns the SMD-solvated DFT energy (and the gas-phase energy for the solvation term).
This supplies the solution-phase energies the redox-potential calculation needs.

By default the geometry is OPTIMIZED IN SOLVENT (SMD gradients, geomeTRIC): UMA is only
a gas-phase pre-optimizer and cannot see solvation, so the final structure/energy must be
relaxed on the SMD potential-energy surface. Pass do_opt=False for a single-point only.

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


def _dft_module(backend: str):
    """Return the dft module for the chosen backend ('gpu' -> gpu4pyscf, else pyscf)."""
    if backend == "gpu":
        from gpu4pyscf import dft as gdft
        return gdft
    from pyscf import dft as cdft
    return cdft


def _geomopt_fn(backend: str):
    """geomeTRIC driver. pyscf's is backend-agnostic and drives a gpu4pyscf mean-field too
    (gpu4pyscf ships no geometric_solver of its own, only an ase_solver)."""
    from pyscf.geomopt.geometric_solver import optimize as fn
    return fn


# --- Decided level of theory (composite protocol; see docs/PLAN.md §9) ---
# Geometry optimized IN SMD at a cheap-but-accurate meta-GGA (r2SCAN-D4); energy scored at
# a range-separated hybrid (wB97M-V: wB97M + VV10 nonlocal dispersion) that controls the
# self-interaction error dominating redox potentials. (wB97X-D3 was the plan but is not
# implemented in gpu4pyscf; wB97M-V is the top-benchmark range-separated substitute.)
# Anionic/reduced states get diffuse functions (physics: a loosely bound extra electron
# needs them). Strings verified against the env by scripts/probe_dft.py before a submit.
OPT_XC = "r2scan"
OPT_DISP = "d4"           # empirical D4 dispersion
OPT_NLC = None
OPT_BASIS = "def2-svp"
OPT_BASIS_ANION = "def2-svpd"
SP_XC = "wb97m-v"
SP_DISP = None
SP_NLC = "vv10"           # VV10 nonlocal correlation (built into wB97M-V)
SP_BASIS = "def2-tzvp"
SP_BASIS_ANION = "def2-tzvpd"


def _basis_for(charge: int, base: str, diffuse: str) -> str:
    """Anions need diffuse functions; neutrals/cations do not."""
    return diffuse if int(charge) < 0 else base


def _build_smd_mf(mol, xc, solvent, backend="cpu", disp=None, nlc=None):
    dft = _dft_module(backend)
    KS = dft.RKS if mol.spin == 0 else dft.UKS
    mf = KS(mol); mf.xc = xc
    if disp:
        mf.disp = disp          # empirical dispersion (e.g. r2SCAN-D4)
    if nlc:
        mf.nlc = nlc            # nonlocal correlation (VV10, for wB97M-V)
    mf = mf.SMD()
    mf.with_solvent.solvent = solvent
    return mf


def _kernel_robust(mf):
    """Run SCF; if it doesn't converge, retry with escalating stabilizers.

    Only the fallbacks fire on hard cases (e.g. open-shell metal cations in SMD);
    a system that converges on the first pass is scored exactly as before, so this
    never perturbs results that were already fine. Same functional/basis/solvent —
    we only make the SCF find its minimum, we do not change the physics.
    """
    mf.max_cycle = 200
    mf.conv_tol = 1e-9
    e = float(mf.kernel())
    if bool(mf.converged):
        return e, True
    # (1) level shift + damping, restarted from the current density
    dm = mf.make_rdm1()
    mf.level_shift = 0.5
    mf.damp = 0.3
    e = float(mf.kernel(dm0=dm))
    if bool(mf.converged):
        return e, True
    # (2) second-order (Newton/SOSCF) restart — robust for stubborn metal SCFs
    try:
        dm = mf.make_rdm1()
        mf.level_shift = 0.0
        mf.damp = 0.0
        mf2 = mf.newton()
        e = float(mf2.kernel(dm0=dm))
        return e, bool(mf2.converged)
    except Exception as exc:  # pragma: no cover - backend-dependent
        print(f"[warn] SOSCF fallback failed: {exc}", flush=True)
        return e, bool(mf.converged)


def _thermal_correction(atoms, charge, mult, xc=None, basis=None, disp=None,
                        backend=None, temperature=298.15):
    """Gibbs thermal free-energy correction (eV) from a GFN2-xTB Hessian (RRHO, 298.15 K).

    Returns dict(g_thermal_eV, zpe_eV, n_imag, freq_min_cm, thermal_method) so redox uses
    G = E_smd + g_thermal instead of the bare electronic energy. g_thermal is xtb's
    'G(RRHO) contrib.' = ZPE + H_thermal - T*S (the part of G beyond the electronic energy).

    Why semi-empirical for the thermal term: (1) RRHO thermal corrections are nearly
    independent of the electronic-structure method, so a cheap Hessian is standard practice
    (DFT electronic energy + GFN2 thermal is a well-established composite); (2) gpu4pyscf's
    UKS analytic Hessian is numerically BROKEN in this version (it inflates open-shell
    frequencies ~2x -> corrupt ZPE), so DFT Hessians are not a safe option for the radical
    states that dominate this dataset. The DFT+SMD wb97m-v energy is untouched; only the
    (ZPE + thermal - TS) term is xTB. xc/basis/disp/backend are accepted for signature
    compatibility and IGNORED (the correction is method-transferable).

    xtb runs GFN2 with the FIXED input geometry (no re-optimization), so the correction is
    evaluated at our DFT+SMD-optimized structure. It uses xtb's free-/rigid-rotor treatment
    of low/imaginary modes, which is robust for floppy species (e.g. ion pairs).
    """
    import subprocess
    import tempfile
    import os
    import re

    q = int(charge)
    uhf = int(mult) - 1
    xtb_bin = os.environ.get("XTB_BIN", "xtb")
    with tempfile.TemporaryDirectory(prefix="xtbhess_") as d:
        _write_xyz(atoms, Path(d) / "mol.xyz")
        cmd = [xtb_bin, "mol.xyz", "--gfn", "2", "--chrg", str(q), "--uhf", str(uhf),
               "--hess", "--acc", "1.0"]
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("XTB_THREADS", "8"))
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env)
        out = r.stdout + "\n" + r.stderr

    def _grab(pat, cast=float):
        m = re.search(pat, out)
        return cast(m.group(1)) if m else None

    g_rrho = _grab(r"G\(RRHO\) contrib\.\s+(-?\d+\.\d+)\s+Eh")   # ZPE + thermal - T*S, Hartree
    zpe = _grab(r"zero point energy\s+(-?\d+\.\d+)\s+Eh")
    n_imag = _grab(r"#\s*imaginary freq\.\s+(\d+)", int)
    if g_rrho is None:
        raise RuntimeError(f"xtb thermal parse failed (rc={r.returncode}): ...{out[-600:]}")
    return dict(g_thermal_eV=g_rrho * HARTREE_EV,
                zpe_eV=zpe * HARTREE_EV if zpe is not None else None,
                n_imag=n_imag,
                freq_min_cm=None,
                thermal_method="gfn2-xtb")


def dft_smd(xyz_path: Path, charge: int, mult: int,
            opt_xc: str = OPT_XC, opt_basis: str = OPT_BASIS,
            opt_basis_anion: str = OPT_BASIS_ANION, opt_disp: str | None = OPT_DISP,
            opt_nlc: str | None = OPT_NLC,
            sp_xc: str = SP_XC, sp_basis: str = SP_BASIS,
            sp_basis_anion: str = SP_BASIS_ANION, sp_disp: str | None = SP_DISP,
            sp_nlc: str | None = SP_NLC,
            solvent: str | None = None, do_gas: bool = True, do_smd: bool = True,
            do_opt: bool = True, do_freq: bool = True, opt_out: Path | None = None,
            max_opt_steps: int = 100, backend: str = "cpu") -> dict:
    """Composite DFT+SMD: optimize geometry in solvent at (opt_xc/opt_basis), then score
    the energy at (sp_xc/sp_basis). Anions use the diffuse basis variant. Returns energies
    (Ha) in SMD solvent and (optionally) gas phase. backend='gpu' uses gpu4pyscf.
    """
    from pyscf import gto
    from ase.io import read
    dft = _dft_module(backend)

    q = int(charge)
    nunpaired = int(mult) - 1  # pyscf mol.spin = 2S = (mult-1)
    obas = _basis_for(q, opt_basis, opt_basis_anion)
    sbas = _basis_for(q, sp_basis, sp_basis_anion)

    atoms = read(str(xyz_path))
    def _mol(basis):
        atom_str = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                             for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
        return gto.M(atom=atom_str, basis=basis, charge=q, spin=nunpaired, verbose=0)

    solvent = solvent or _electrolyte().SOLVENT["name"]

    out = dict(xyz=str(xyz_path), charge=q, mult=int(mult), solvent=solvent,
               optimized=bool(do_opt), backend=backend,
               opt_xc=opt_xc, opt_basis=obas, opt_disp=opt_disp, opt_nlc=opt_nlc,
               sp_xc=sp_xc, sp_basis=sbas, sp_disp=sp_disp, sp_nlc=sp_nlc,
               # xc/basis mirror the ENERGY level (what redox referencing keys on)
               xc=sp_xc, basis=sbas)

    # Optimize the geometry IN SOLVENT (SMD gradients) at the opt level.
    mol = _mol(obas if do_opt else sbas)
    if do_opt:
        geomopt = _geomopt_fn(backend)
        mf_opt = _build_smd_mf(mol, opt_xc, solvent, backend, opt_disp, opt_nlc)
        mol = geomopt(mf_opt, maxsteps=max_opt_steps)  # relaxed Mole (at opt basis)
        opt_atoms = _mol_to_atoms(mol)
        if opt_out is not None:
            _write_xyz(opt_atoms, Path(opt_out),
                       comment=f"SMD({solvent})-opt {opt_xc}/{obas} q={q} m={mult}")
            out["opt_out"] = str(opt_out)
        # rebuild at the (larger) single-point basis on the optimized coordinates
        atoms = opt_atoms

    # Energy: solvated single point at the sp level on the optimized geometry.
    # do_smd=False skips the (expensive) SMD SCF for callers that only need the gas energy
    # (e.g. inner-sphere reorganization cross-points) — a ~2x saving there.
    mol_sp = _mol(sbas)
    e_smd = None
    if do_smd:
        mfs = _build_smd_mf(mol_sp, sp_xc, solvent, backend, sp_disp, sp_nlc)
        e_smd, conv_smd = _kernel_robust(mfs)
        out.update(e_smd_Ha=e_smd, e_smd_eV=e_smd * HARTREE_EV, converged_smd=conv_smd)

    if do_gas:
        KS = dft.RKS if nunpaired == 0 else dft.UKS
        mfg = KS(mol_sp); mfg.xc = sp_xc
        if sp_disp:
            mfg.disp = sp_disp
        if sp_nlc:
            mfg.nlc = sp_nlc
        e_gas, conv_gas = _kernel_robust(mfg)
        out.update(e_gas_Ha=e_gas, e_gas_eV=e_gas * HARTREE_EV, converged_gas=conv_gas)
        if e_smd is not None:
            out["dG_solv_eV"] = (e_smd - e_gas) * HARTREE_EV

    # Thermal free-energy correction (GFN2-xTB RRHO on the DFT+SMD-optimized geometry) so
    # redox uses G = E_smd + g_thermal, not the bare electronic energy. See _thermal_correction.
    if do_freq:
        try:
            th = _thermal_correction(atoms, q, int(mult))
            out.update(**th, freq_level="gfn2-xtb (RRHO, 298.15K)")
        except Exception as exc:  # a failed Hessian must not lose the (expensive) energies
            print(f"[warn] freq/thermal failed: {type(exc).__name__}: {str(exc)[:120]}",
                  flush=True)
            out["g_thermal_eV"] = None
    return out


def _mol_to_atoms(mol):
    from ase import Atoms
    BOHR = 0.52917721067
    syms = [mol.atom_symbol(i) for i in range(mol.natm)]
    coords = mol.atom_coords() * BOHR  # pyscf stores coords in Bohr
    return Atoms(symbols=syms, positions=coords)


def _write_xyz(atoms, path: Path, comment: str = ""):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(atoms)), comment]
    for s, p in zip(atoms.get_chemical_symbols(), atoms.positions):
        lines.append(f"{s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def read_manifest():
    import csv
    with (ROOT / "library" / "manifest.csv").open() as f:
        return list(csv.DictReader(f))


def geom_for(gid: str, state: str) -> Path:
    """Prefer the UMA-relaxed geometry; fall back to the conformer seed."""
    relaxed = ROOT / "calcs" / "uma" / gid / state / "relaxed.xyz"
    return relaxed if relaxed.exists() else ROOT / "library" / gid / f"{state}.xyz"


def _uma_mult(gid: str, state: str, fallback: int) -> int:
    """Use the spin multiplicity UMA actually chose (spin scan), not the manifest hint."""
    uj = ROOT / "calcs" / "uma" / gid / state / "result.json"
    if uj.exists():
        try:
            return int(json.loads(uj.read_text())["chosen_mult"])
        except Exception:
            pass
    return fallback


def run_batch(rows, force, do_opt=True, do_freq=True, backend="cpu"):
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
        fallback_mult = int(r.get("mult_hint") or r.get("mult") or 1)
        mult = _uma_mult(gid, st, fallback_mult)
        print(f"[run ] {gid}/{st} q={r['charge']} m={mult} opt={do_opt} "
              f"backend={backend} ...", flush=True)
        try:
            outdir.mkdir(parents=True, exist_ok=True)
            res = dft_smd(geom, int(r["charge"]), mult,
                          do_opt=do_opt, do_freq=do_freq,
                          opt_out=outdir / "opt.xyz", backend=backend)
        except Exception as e:
            print(f"[fail] {gid}/{st}: {type(e).__name__}: {str(e)[:120]}", flush=True)
            continue
        res.update(id=gid, state=st, n_e=int(r["n_e"]))
        rj.write_text(json.dumps(res, indent=2))
        gth = res.get("g_thermal_eV")
        gth_s = f"{gth:+.3f}" if isinstance(gth, (int, float)) else "n/a"
        print(f"[done] {gid}/{st} E_smd={res['e_smd_eV']:.3f} eV "
              f"dGsolv={res.get('dG_solv_eV', float('nan')):.3f} eV "
              f"Gtherm={gth_s} eV (imag={res.get('n_imag', '-')}) "
              f"conv={res['converged_smd']} ({res['sp_xc']}/{res['sp_basis']})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xyz", default=None, help="single-geometry mode")
    ap.add_argument("--charge", type=int)
    ap.add_argument("--mult", type=int)
    ap.add_argument("--all", action="store_true", help="batch over manifest states")
    ap.add_argument("--only", default=None, help="'group' or 'group:state'")
    ap.add_argument("--shard", default=None, help="'n:i' for CPU fan-out")
    ap.add_argument("--out", default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-opt", dest="no_opt", action="store_true",
                    help="single-point only; skip in-solvent geometry optimization")
    ap.add_argument("--no-freq", dest="no_freq", action="store_true",
                    help="skip the gas-phase Hessian / thermal free-energy correction")
    ap.add_argument("--backend", default="cpu", choices=["cpu", "gpu"],
                    help="'gpu' uses gpu4pyscf (set CUDA_VISIBLE_DEVICES to pin a GPU)")
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
        run_batch(rows, args.force, do_opt=not args.no_opt,
                  do_freq=not args.no_freq, backend=args.backend)
        return

    res = dft_smd(Path(args.xyz), args.charge, args.mult,
                  do_opt=not args.no_opt, do_freq=not args.no_freq, backend=args.backend,
                  opt_out=(Path(args.out).parent / "opt.xyz") if args.out else None)
    print(json.dumps(res, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
