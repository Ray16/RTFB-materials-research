"""DFT + SMD(acetonitrile) energies (the solvation step).

Takes a geometry (e.g. a UMA-relaxed xyz) plus explicit charge + spin multiplicity and
returns the SMD-solvated DFT energy (and the gas-phase energy for the solvation term).
This supplies the solution-phase energies the redox-potential calculation needs.

By default the geometry is OPTIMIZED IN SOLVENT (SMD gradients, geomeTRIC): UMA is only
a gas-phase pre-optimizer and cannot see solvation, so the final structure/energy must be
relaxed on the SMD potential-energy surface. Pass do_opt=False for a single-point only.

Level of theory (decided; see `redox.core.protocol.ACTIVE_SP`): geometry optimized IN SMD at
r2SCAN-D4/def2-SVP(D); production energies at wB97M-V/def2-TZVPD (VV10 nonlocal dispersion)
for every charge state and Fc — a dispersion-corrected range-separated hybrid with a uniform
diffuse basis.

  python -m redox.qm.dft --xyz calcs/uma/pyridinium/ox/relaxed.xyz --charge 1 --mult 1
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from redox.core.common import HARTREE_EV, ROOT, load_config, read_manifest
from redox.core.common import write_xyz as _write_xyz


def _electrolyte():
    return load_config("electrolyte")


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
# GEOMETRY: anionic states are optimized with diffuse functions (def2-SVPD), others def2-SVP.
# ENERGY: ONE diffuse basis (def2-TZVPD) for EVERY charge state and for Fc/Fc+ — a ladder
# scored with TZVP on the neutral and TZVPD on the anion has the diffuse freedom on one side
# only, so neither E° nor dG_disp gets clean basis-set cancellation. The energy protocol is
# defined once in redox.core.protocol.ACTIVE_SP; the SP_* constants below mirror it.
OPT_XC = "r2scan"
OPT_DISP = "d4"           # empirical D4 dispersion
OPT_NLC = None
OPT_BASIS = "def2-svp"
OPT_BASIS_ANION = "def2-svpd"
from redox.core.protocol import ACTIVE_SP as _SP
SP_XC = _SP["xc"]
SP_DISP = _SP["disp"]
SP_NLC = _SP["nlc"]       # VV10 nonlocal correlation (built into wB97M-V)
SP_BASIS = _SP["basis"]   # uniform diffuse basis for all charge states
SP_BASIS_ANION = _SP["basis"]


def _basis_for(charge: int, base: str, diffuse: str) -> str:
    """Anions need diffuse functions; neutrals/cations do not."""
    return diffuse if int(charge) < 0 else base


# def2 replaces the core with an effective core potential from Rb (Z=37) onward, so for those
# elements the def2 basis is a VALENCE-ONLY set. Building such an atom without also loading the
# matching ECP silently produces an all-electron calculation in a basis that has no core
# functions (e.g. W: 74 electrons crammed into the 40 AOs meant for 14) — it converges and
# returns a number, so nothing downstream catches it. Attach the ECP whenever a heavy element
# is present; light-only systems pass ecp=None so their (correct) all-electron treatment and
# existing results are bit-for-bit unchanged.
_ECP_MIN_Z = 37


def _ecp_for(symbols, basis: str):
    """def2 ECP name if any atom needs one, else None."""
    from pyscf.data.elements import charge as _z
    return basis if any(_z(s) >= _ECP_MIN_Z for s in symbols) else None


def _build_smd_mf(mol, xc, solvent, backend="cpu", disp=None, nlc=None):
    dft = _dft_module(backend)
    KS = dft.RKS if mol.spin == 0 else dft.UKS
    # RI-J (density fitting) by default: ~3-5x faster; the DF error on the absolute energy is
    # ~1 meV and cancels to sub-meV in the energy DIFFERENCES this pipeline reports (redox
    # Delta-G, reorg lambda). density_fit() MUST be applied BEFORE .SMD() (solvent wraps the DF
    # SCF; the reverse order raises in gpu4pyscf).
    mf = KS(mol).density_fit().SMD()
    mf.xc = xc
    if disp:
        mf.disp = disp          # empirical dispersion (e.g. r2SCAN-D4)
    if nlc:
        mf.nlc = nlc            # nonlocal correlation (VV10, for wB97M-V)
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


# Open-shell SCF can converge to a higher solution from the default (minao) guess: ferrocenium
# gas-phase landed 2.5 eV high, an OROP radical anion 0.11 eV high in SMD
# (scripts/validation/reference/scf_level_check.py). Open shells are therefore solved from
# every guess here and the lowest converged energy is kept; closed shells use minao only.
OPEN_SHELL_GUESSES = ("minao", "atom", "huckel")


def _free_gpu_pool():
    try:
        import cupy
        cupy.get_default_memory_pool().free_all_blocks()
    except Exception:          # CPU backend / no cupy: nothing to free
        pass


def _lowest_scf(make_mf, guesses, post=None):
    """Run `make_mf()` from each initial guess; return (e, converged, post_result, per_guess)
    for the lowest CONVERGED energy (falls back to the first guess if none converged).
    `post(mf)` (e.g. orbital diagnostics) is evaluated right after each solve so no SCF object
    outlives its guess — holding two UKS objects at once ran a V100-32GB out of memory."""
    best, per = None, []
    for g in guesses:
        mf = make_mf()
        mf.init_guess = g
        e, conv = _kernel_robust(mf)
        per.append(dict(guess=g, e_Ha=e, converged=conv))
        if best is None or (conv and (not best[1] or e < best[0])):
            best = (e, conv, post(mf) if post else None)
        del mf
        _free_gpu_pool()
    return best[0], best[1], best[2], per


def _orbital_diag(mf, nunpaired):
    """Cheap QC from a converged mean-field: HOMO/LUMO (eV) and <S^2>.

    For an ANION, homo_eV > 0 means the extra electron is UNBOUND (the gas-phase state is
    unphysical, so a gas single point there is meaningless — score it in SMD instead). For an
    open shell, spin_contam = <S^2> - S(S+1) flags a spin-contaminated / wrong-state SCF.
    """
    import numpy as np
    def _cpu(a):
        try:
            return np.asarray(a.get())
        except Exception:
            return np.asarray(a)
    out = {}
    try:
        mo_e, mo_o = _cpu(mf.mo_energy), _cpu(mf.mo_occ)
        occ = mo_e[mo_o > 0]; vir = mo_e[mo_o == 0]
        if occ.size:
            out["homo_eV"] = float(occ.max()) * HARTREE_EV
        if vir.size:
            out["lumo_eV"] = float(vir.min()) * HARTREE_EV
    except Exception:
        pass
    try:
        ss = mf.spin_square()
        s2 = float(ss[0] if isinstance(ss, (tuple, list)) else ss)
        S = nunpaired / 2.0
        out["s_squared"] = round(s2, 4)
        out["spin_contam"] = round(s2 - S * (S + 1), 4)
    except Exception:
        pass
    return out


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


def _xtb_alpb_opt(atoms, charge: int, mult: int, solvent: str = "acetonitrile",
                  accuracy: float = 1.0):
    """Fast SOLVATED pre-optimization with GFN2-xTB + ALPB(<solvent>). Returns
    (ASE Atoms of the ALPB-optimized geometry, ok_flag); on any failure returns the input
    geometry unchanged with ok=False (never loses the seed).

    Purpose: warm-start the DFT+SMD optimization from a geometry that already feels the
    solvent, removing the gas-phase-collapse basin-trapping risk a gas-only (UMA) seed carries
    for charged/floppy species. ALPB also binds otherwise gas-unbound dianions, so the pre-opt
    is well-defined for the -2 states. Cheap enough to run on every seed.
    """
    import subprocess
    import tempfile
    import os
    from ase.io import read as ase_read
    q = int(charge)
    uhf = int(mult) - 1
    xtb_bin = os.environ.get("XTB_BIN", "xtb")
    with tempfile.TemporaryDirectory(prefix="xtbalpb_") as d:
        _write_xyz(atoms, Path(d) / "mol.xyz")
        cmd = [xtb_bin, "mol.xyz", "--gfn", "2", "--chrg", str(q), "--uhf", str(uhf),
               "--opt", "--alpb", solvent, "--acc", str(accuracy)]
        env = dict(os.environ, OMP_NUM_THREADS=os.environ.get("XTB_THREADS", "2"))
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True, env=env)
        optxyz = Path(d) / "xtbopt.xyz"
        if r.returncode == 0 and optxyz.exists():
            return ase_read(str(optxyz)), True
    return atoms, False


def conformer_seeds(gid: str, state: str, n: int) -> list:
    """Starting geometries for a state's DFT+SMD optimization.

    n<=1 (DEFAULT): the single BEST seed = the UMA-relaxed geometry (charge/spin-aware, the
    lowest-energy conformer UMA already selected for this state), falling back to the FF-lowest
    conformer, then the per-state library seed. For the rigid aromatic redox cores here the
    conformers collapse to one solvated minimum (measured DFT seed spread <5 meV in 13/14
    candidate states), so one seed is enough.

    n>1 (opt-in, --nconf): the top-n FF-ranked conformers, multi-conformer screened (keep the
    lowest DFT+SMD energy). Reserve for genuinely floppy species / reduced dianions, where the
    seed choice can matter (~0.05 eV in the one case it did, pmdi/red2)."""
    if n <= 1:
        g = geom_for(gid, state)          # UMA-relaxed best seed (charge/spin-aware)
        if g.exists():
            return [g]
    confdir = ROOT / "library" / gid / "conformers"
    confs = sorted(confdir.glob("conf_*.xyz"))[:n]
    if confs:
        return confs
    g = geom_for(gid, state)
    return [g] if g.exists() else []


def dft_smd(xyz_path: Path, charge: int, mult: int,
            opt_xc: str = OPT_XC, opt_basis: str = OPT_BASIS,
            opt_basis_anion: str = OPT_BASIS_ANION, opt_disp: str | None = OPT_DISP,
            opt_nlc: str | None = OPT_NLC,
            sp_xc: str = SP_XC, sp_basis: str = SP_BASIS,
            sp_basis_anion: str = SP_BASIS_ANION, sp_disp: str | None = SP_DISP,
            sp_nlc: str | None = SP_NLC,
            solvent: str | None = None, do_gas: bool = True, do_smd: bool = True,
            do_opt: bool = True, do_freq: bool = True, opt_out: Path | None = None,
            max_opt_steps: int = 100, constraints: str | None = None,
            conv_params: dict | None = None, assert_convergence: bool = True,
            backend: str = "cpu") -> dict:
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
        syms = atoms.get_chemical_symbols()
        atom_str = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                             for s, p in zip(syms, atoms.positions))
        return gto.M(atom=atom_str, basis=basis, ecp=_ecp_for(syms, basis),
                     charge=q, spin=nunpaired, verbose=0)

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
        # `constraints` (a geomeTRIC constraints file, e.g. frozen dihedrals) keeps a flexible
        # anion in the neutral's conformer -> a conformer-matched inner-sphere lambda.
        # `conv_params` (convergence_gmax/grms/energy...) + assert_convergence=False let a
        # constrained opt stop at a loose-but-energy-converged geometry instead of oscillating.
        _optkw = dict(maxsteps=max_opt_steps, assert_convergence=assert_convergence)
        if constraints:
            _optkw["constraints"] = constraints
        if conv_params:
            _optkw.update(conv_params)
        mol = geomopt(mf_opt, **_optkw)  # relaxed Mole (at opt basis)
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
    guesses = OPEN_SHELL_GUESSES if nunpaired else ("minao",)
    # provenance of the energy: what an adopter must match (redox.qm.sp._adoptable)
    from redox.core.protocol import software_versions
    out.update(density_fit=True, conv_tol=1e-9, scf_guesses=list(guesses),
               software=software_versions())
    e_smd = None
    if do_smd:
        e_smd, conv_smd, _, per = _lowest_scf(
            lambda: _build_smd_mf(mol_sp, sp_xc, solvent, backend, sp_disp, sp_nlc), guesses)
        out.update(e_smd_Ha=e_smd, e_smd_eV=e_smd * HARTREE_EV, converged_smd=conv_smd,
                   scf_guess_energies_smd=per)

    if do_gas:
        KS = dft.RKS if nunpaired == 0 else dft.UKS

        def _gas_mf():
            m = KS(mol_sp).density_fit(); m.xc = sp_xc   # RI-J: ~3-5x faster (see _build_smd_mf)
            if sp_disp:
                m.disp = sp_disp
            if sp_nlc:
                m.nlc = sp_nlc
            return m
        e_gas, conv_gas, diag, per = _lowest_scf(_gas_mf, guesses,
                                                 post=lambda m: _orbital_diag(m, nunpaired))
        out.update(e_gas_Ha=e_gas, e_gas_eV=e_gas * HARTREE_EV, converged_gas=conv_gas,
                   scf_guess_energies_gas=per)
        # QC: anion-bound check (gas HOMO > 0 => unbound) + spin contamination.
        for k, v in (diag or {}).items():
            out[f"gas_{k}"] = v
        if out.get("gas_homo_eV") is not None and q < 0:
            out["anion_unbound"] = bool(out["gas_homo_eV"] > 0.0)
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


def _opt_protocol_mismatch(res: dict, charge: int) -> str:
    """'' if an existing result.json was optimized at the current GEOMETRY protocol, else a
    description of the mismatch. (Energies are protocol-addressed separately: redox.qm.sp.)"""
    want = dict(opt_xc=OPT_XC, opt_basis=_basis_for(charge, OPT_BASIS, OPT_BASIS_ANION),
                opt_disp=OPT_DISP, solvent=_electrolyte().SOLVENT["name"])
    bad = [f"{k}={res.get(k)!r}!={v!r}" for k, v in want.items()
           if (str(res.get(k)).lower() if res.get(k) is not None else None)
           != (str(v).lower() if v is not None else None)]
    if not res.get("optimized"):
        bad.append("not optimized")
    return "; ".join(bad)


def run_batch(rows, force, do_opt=True, do_freq=True, backend="cpu", nconf=1, preopt="alpb",
              torsion_scan=False):
    """For each state: multi-conformer DFT+SMD screening. Take up to `nconf` conformer seeds,
    optionally pre-optimize each in solvent with xtb-ALPB (`preopt`), DFT+SMD-optimize every
    seed, and keep the LOWEST-E_smd result (the winning solvated conformer). Thermal (RRHO) is
    computed only on the winner (RRHO is ~conformer-independent → saves N-1 xtb Hessians).
    Resumable: a state whose result.json exists is skipped.

    torsion_scan (opt-in, default OFF): for floppy molecules, replace the FF/UMA seeds with
    geometries rotated around the primary ring-ring torsion (>= max(nconf, 8) angles), so the
    lowest-E_smd selection finds the global-min conformer. Rigid molecules (no such torsion)
    fall back to normal seeding automatically. See src/redox/qm/torsion_scan.py."""
    import shutil
    from ase.io import read as ase_read
    outroot = ROOT / "calcs" / "dft"
    for r in rows:
        gid, st = r["id"], r["state"]
        outdir = outroot / gid / st
        rj = outdir / "result.json"
        if rj.exists() and not force:
            stale = _opt_protocol_mismatch(json.loads(rj.read_text()), int(r["charge"]))
            if stale:
                # never silently accept a geometry from a different optimization protocol,
                # and never silently overwrite it either: needs an explicit --force.
                print(f"[STALE] {gid}/{st}: result.json optimized at a different protocol "
                      f"({stale}); keeping it — rerun with --force to replace", flush=True)
            else:
                print(f"[skip] {gid}/{st}", flush=True)
            continue
        q = int(r["charge"])
        seeds = conformer_seeds(gid, st, nconf) if do_opt else [geom_for(gid, st)]
        seeds = [s for s in seeds if s and s.exists()]
        if not seeds:
            print(f"[miss] {gid}/{st} no geometry", flush=True); continue
        if do_opt and torsion_scan:
            outdir.mkdir(parents=True, exist_ok=True)
            from redox.qm.torsion_scan import torsion_rotated_seeds
            rot = torsion_rotated_seeds(seeds[0], max(nconf, 8), outdir)
            if len(rot) > 1:                       # a ring-ring torsion was found; else no-op
                seeds = rot
                print(f"[tscan] {gid}/{st}: {len(seeds)} torsion-rotated seeds", flush=True)
        fallback_mult = int(r.get("mult_hint") or r.get("mult") or 1)
        mult = _uma_mult(gid, st, fallback_mult)
        print(f"[run ] {gid}/{st} q={q} m={mult} seeds={len(seeds)} preopt={preopt} "
              f"opt={do_opt} backend={backend} ...", flush=True)
        outdir.mkdir(parents=True, exist_ok=True)
        best = None            # (e_smd_eV, seed_idx, res_dict, opt_path)
        seed_es = []
        for si, seed in enumerate(seeds):
            try:
                atoms0 = ase_read(str(seed))
                if do_opt and preopt == "alpb":
                    atoms0, ok = _xtb_alpb_opt(atoms0, q, mult)  # solvated warm start
                seed_tmp = outdir / f"_seed_{si:02d}.xyz"
                _write_xyz(atoms0, seed_tmp, comment=f"seed {si} preopt={preopt}")
                opt_tmp = outdir / f"_opt_{si:02d}.xyz"
                res = dft_smd(seed_tmp, q, mult, do_opt=do_opt, do_freq=False,
                              opt_out=opt_tmp, backend=backend)
                e = res.get("e_smd_eV")
                seed_es.append(e)
                if e is not None and (best is None or e < best[0]):
                    best = (e, si, res, opt_tmp if opt_tmp.exists() else seed_tmp)
            except Exception as e:
                print(f"[warn] {gid}/{st} seed {si} failed: {type(e).__name__}: "
                      f"{str(e)[:100]}", flush=True)
                seed_es.append(None)
                continue
        if best is None:
            print(f"[fail] {gid}/{st}: all {len(seeds)} seeds failed", flush=True); continue
        e_best, si_best, res, optpath = best
        # thermal (RRHO) only on the winning conformer
        if do_freq:
            try:
                th = _thermal_correction(ase_read(str(optpath)), q, mult)
                res.update(**th, freq_level="gfn2-xtb (RRHO, 298.15K)")
            except Exception as exc:
                print(f"[warn] {gid}/{st} thermal failed: {str(exc)[:100]}", flush=True)
                res["g_thermal_eV"] = None
        shutil.copy(optpath, outdir / "opt.xyz"); res["opt_out"] = str(outdir / "opt.xyz")
        valid_es = [x for x in seed_es if x is not None]
        res.update(id=gid, state=st, n_e=int(r["n_e"]),
                   n_seeds=len(seeds), seed_preopt=preopt, winning_seed=si_best,
                   seed_e_smd_eV=[round(x, 4) if x is not None else None for x in seed_es],
                   seed_spread_eV=(round(max(valid_es) - min(valid_es), 4)
                                   if len(valid_es) > 1 else 0.0))
        rj.write_text(json.dumps(res, indent=2))
        for f in outdir.glob("_seed_*.xyz"):
            f.unlink()
        for f in outdir.glob("_opt_*.xyz"):
            f.unlink()
        gth = res.get("g_thermal_eV")
        gth_s = f"{gth:+.3f}" if isinstance(gth, (int, float)) else "n/a"
        print(f"[done] {gid}/{st} E_smd={e_best:.3f} eV (winner {si_best}/{len(seeds)}, "
              f"spread {res['seed_spread_eV']:.3f} eV) Gtherm={gth_s} "
              f"conv={res.get('converged_smd')} ({res['sp_xc']}/{res['sp_basis']})", flush=True)


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
    ap.add_argument("--nconf", type=int, default=1,
                    help="# conformer seeds to DFT+SMD-optimize per state, keep lowest. "
                         "Default 1 = single best UMA-relaxed seed (rigid cores need no more). "
                         "Set >1 to multi-conformer screen floppy species / dianions.")
    ap.add_argument("--preopt", default="alpb", choices=["alpb", "none"],
                    help="solvated pre-opt of each seed before DFT (alpb=xtb-ALPB(MeCN), none=off)")
    ap.add_argument("--torsion-scan", dest="torsion_scan", action="store_true",
                    help="opt-in: seed the multi-conformer search with geometries rotated around "
                         "the primary ring-ring torsion (>= max(nconf,8) angles) to find the "
                         "global-min conformer. For floppy species (viologens, biaryl quinones); "
                         "no-op on rigid molecules. OFF by default.")
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
                  do_freq=not args.no_freq, backend=args.backend,
                  nconf=args.nconf, preopt=args.preopt, torsion_scan=args.torsion_scan)
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
