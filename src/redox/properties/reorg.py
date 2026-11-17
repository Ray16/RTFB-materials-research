"""Inner-sphere reorganization energy lambda_i (Nelsen 4-point) for each redox couple.

For a couple O + e- -> R (O = higher charge, R = one less), on the two adiabatic SOLUTION
surfaces:

    lambda_i = [E_O(q_R) - E_O(q_O)] + [E_R(q_O) - E_R(q_R)]

where E_X(q_Y) is the SMD (solution) energy of species X (its own charge+spin) evaluated at the
SMD-OPTIMIZED geometry of species Y. Two of the four points are each state's OWN SMD energy
(already stored); the two CROSS points E_O(q_R) and E_R(q_O) are extra SMD single points at the
partner geometry.

Convention (SOLUTION-PHASE inner-sphere lambda_i)
------------------------------------------------
Every point is scored on ONE self-consistent surface: the SMD(MeCN) potential-energy surface
at the active protocol (wb97m-v / def2-TZVPD). Each own-geometry point is therefore the MINIMUM
of its own SMD surface, so the two directional distortion penalties are nonnegative BY
CONSTRUCTION (lambda_i >= 0 is structural, not a QC hope). This replaces the previous mixed
surface, which scored GAS energies on SMD geometries: those own-geometry gas points were not
stationary (distortion could go negative) and, for the n-type anion pool that dominates this
dataset, the gas anion is frequently UNBOUND so the gas energy was ill-defined. The gas value is
kept SIDE BY SIDE (lambda_i_gas_eV) for one transition only, as a diagnostic of the shift.

Orthogonality to the outer-sphere lambda_o: lambda_i is the energy to move the NUCLEAR
GEOMETRY at FIXED charge between the two equilibrium structures on the two solution surfaces;
the nonequilibrium continuum lambda_o (redox.properties.lambda_outer) is the slow-polarization
energy to move the SOLVENT CHARGE distribution at FIXED geometry. The two coordinates
(geometry-at-fixed-charge vs charge-at-fixed-geometry) are distinct, so the terms do not
double-count. Do NOT add a solvent term to lambda_i, and do NOT add lambda_i's relaxation to a
frozen-slow solution total a second time.

Matched-conformer averaging
---------------------------
When a per-state DFT+SMD conformer ensemble is available on disk
(calcs/dft/<id>/<state>/conformers/conf_*.xyz, optional sidecar conf_*.json with e_smd_eV),
lambda_i is computed per MATCHED conformer family (same family index for O and R) and
Boltzmann-weighted by the couple's per-family solution energy. If O and R optimized to
different conformers, the families are re-paired (lowest-to-lowest) and the cross point is
evaluated from the PARTNER's matching conformer rather than discarding the couple. If no
ensemble is on disk the single matched minimum (opt.xyz) is used and conformer_averaged is
recorded False.

A QC-clean lambda_i does not establish a solution rate, reversibility or decomposition
stability. Those require electronic-coupling and reaction-path evidence.

  gpu_reserve run 0 --label redox-reorg -- env PYTHONPATH=src python -m redox.properties.reorg --only viologen --backend gpu
  PYTHONPATH=src python -m redox.properties.reorg --aggregate
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import inspect
import json
import math

from redox.core.common import DFT, EV_KJ, EV_MEV, KT_EV, RESULTS, ROOT, state_names, group_ids
from redox.core.common import read_result as _res
from redox.core.protocol import ACTIVE_SP, protocol_hash
from redox.core.io import atomic_write_json

# --- cross-point cache validity -------------------------------------------------------------
# v4: the STORED QUANTITY changed from a GAS single point to an SMD (solution) single point.
# A cache entry is reusable only if it records the same SCHEMA, the same level of theory, the
# same charge/multiplicity, the same geometry (content hash of the xyz it was evaluated on),
# a CONVERGED SMD SCF, and every required diagnostic field. A v3 (gas) entry is never reused
# for the SMD surface — its gas energy is a different physical quantity — but its gas energy is
# still READ, when present, for the side-by-side lambda_i_gas_eV diagnostic.
CACHE_SCHEMA = 4          # v4: SMD cross-point energy (was v3 gas)
REQUIRED_DIAG = ("smd_s_squared", "anion_unbound")


def _geom_hash(path):
    """Content hash of a geometry file — ties a cache entry to the exact coordinates used."""
    try:
        return hashlib.sha1(path.read_bytes()).hexdigest()[:16]
    except OSError:
        return None


def _cache_ok(d, q, m, geom_hash, proto=ACTIVE_SP):
    """True if a cached cross-point entry is reusable under the current protocol: same schema,
    level of theory (xc + basis + NLC/dispersion), charge/mult, geometry, diagnostics, and a
    converged SMD SCF. The primary stored quantity is now e_smd_eV."""
    if not isinstance(d, dict) or not geom_hash or d.get("e_smd_eV") is None:
        return False
    if not isinstance(d["e_smd_eV"], (int, float)) or not math.isfinite(d["e_smd_eV"]):
        return False
    if d.get("protocol") != proto:
        return False
    if (d.get("sp_xc") or "").lower() != proto["xc"] or (d.get("sp_basis") or "").lower() != proto["basis"]:
        return False
    if (d.get("sp_nlc") or None) != proto["nlc"] or (d.get("sp_disp") or None) != proto["disp"]:
        return False
    if d.get("converged_smd") is not True:      # missing flag = unverified = not usable
        return False
    if d.get("cache_schema") != CACHE_SCHEMA:
        return False
    if any(k not in d for k in REQUIRED_DIAG):
        return False
    if d.get("charge") != q or d.get("mult") != m:
        return False
    if d.get("geom_sha1") != geom_hash:
        return False
    return True


def _couples(gid):
    """Adjacent-charge couples [(O_state, qO, mO), (R_state, qR, mR)] for a group."""
    states = []
    for st in state_names(gid):
        r = _res(gid, st, raw=True)
        states.append((st, int(r["charge"]), int(r["mult"])))
    states.sort(key=lambda x: -x[1])
    out = []
    for (sO, qO, mO), (sR, qR, mR) in zip(states, states[1:]):
        if qO - qR == 1:
            out.append(((sO, qO, mO), (sR, qR, mR)))
    return out


# ------------------------------------------------------------------ conformer ensemble access
def _own_smd_energy(gid, state):
    """The state's OWN converged SMD energy (eV) at the active protocol, or None."""
    r = _res(gid, state)
    e = r.get("e_smd_eV") if r else None
    return float(e) if isinstance(e, (int, float)) and math.isfinite(e) else None


def _own_gas_energy(gid, state):
    """The state's OWN converged gas energy (eV) at the active protocol, or None.
    Kept only for the side-by-side gas diagnostic; may be ill-defined for unbound anions."""
    r = _res(gid, state)
    e = r.get("e_gas_eV") if r else None
    return float(e) if isinstance(e, (int, float)) and math.isfinite(e) else None


def _state_conformers(gid, state):
    """Ordered conformer ensemble for a state.

    Returns a list of dicts {key, geom (Path), e_smd_eV, e_gas_eV} ordered by increasing
    solution energy. If a per-state DFT+SMD ensemble is on disk
    (calcs/dft/<id>/<state>/conformers/conf_*.xyz, optional conf_*.json with e_smd_eV/e_gas_eV)
    it is used; otherwise the single matched minimum (opt.xyz + the stored own energy) is
    returned as the one family 'min' (conformer averaging falls back, recorded by the caller).
    """
    cdir = DFT / gid / state / "conformers"
    confs = []
    if cdir.is_dir():
        for xyz in sorted(cdir.glob("conf_*.xyz")):
            key = xyz.stem
            e_smd = e_gas = None
            side = xyz.with_suffix(".json")
            if side.exists():
                try:
                    sd = json.loads(side.read_text())
                    e_smd = sd.get("e_smd_eV")
                    e_gas = sd.get("e_gas_eV")
                except Exception:
                    pass
            confs.append(dict(key=key, geom=xyz,
                              e_smd_eV=e_smd if isinstance(e_smd, (int, float)) else None,
                              e_gas_eV=e_gas if isinstance(e_gas, (int, float)) else None))
    if confs:
        confs.sort(key=lambda c: (c["e_smd_eV"] is None, c["e_smd_eV"] if c["e_smd_eV"] is not None else 0.0))
        return confs, True
    # fallback: single matched minimum
    return [dict(key="min", geom=DFT / gid / state / "opt.xyz",
                 e_smd_eV=_own_smd_energy(gid, state),
                 e_gas_eV=_own_gas_energy(gid, state))], False


def match_conformers(confs_O, confs_R):
    """Pair O and R conformer families for matched 4-point lambda_i.

    Pure, list-in / list-out so it is unit-testable without any QM. Each conf is a dict with a
    'key'. Rule:
      * families sharing a key (same conformer family on both surfaces) are paired directly;
      * otherwise the two ensembles are RE-PAIRED lowest-to-lowest (both lists are energy
        ordered), so a couple whose O and R minima fell in different families is averaged over
        re-paired families rather than discarded.
    Returns (pairs, repaired) where pairs is a list of (idx_O, idx_R, label)."""
    keys_O = {c["key"]: i for i, c in enumerate(confs_O)}
    keys_R = {c["key"]: i for i, c in enumerate(confs_R)}
    shared = [k for k in keys_O if k in keys_R]
    if shared:
        return [(keys_O[k], keys_R[k], k) for k in shared], False
    n = min(len(confs_O), len(confs_R))
    return [(i, i, f"repair_{i}") for i in range(n)], True


def _boltzmann_weights(energies_eV, kt_eV=KT_EV):
    """Boltzmann weights for a list of conformer energies (eV); None entries drop out.
    A single (or all-None) ensemble returns uniform weights that sum to 1."""
    vals = [e for e in energies_eV if isinstance(e, (int, float)) and math.isfinite(e)]
    if len(vals) <= 1 or not vals:
        n = len(energies_eV)
        return [1.0 / n] * n if n else []
    emin = min(vals)
    raw = [math.exp(-(e - emin) / kt_eV) if isinstance(e, (int, float)) and math.isfinite(e)
           else 0.0 for e in energies_eV]
    z = sum(raw)
    return [r / z for r in raw] if z > 0 else [1.0 / len(raw)] * len(raw)


def _four_point_lambda(E_O_at_O, E_O_at_R, E_R_at_R, E_R_at_O):
    """Nelsen 4-point inner-sphere lambda on one surface (all four energies in eV, same level).
    Returns (lambda_i, relax_ox, relax_red). On the self-consistent SMD surface each own point
    is a minimum, so relax_ox, relax_red >= 0 up to SCF noise."""
    relax_ox = E_O_at_R - E_O_at_O     # distortion of O onto R's geometry
    relax_red = E_R_at_O - E_R_at_R    # distortion of R onto O's geometry
    return relax_ox + relax_red, relax_ox, relax_red


# ---------------------------------------------------------------- cross-point single points
def _cache_path(gid, species_state, geom_path, q, m, tag=None):
    """Protocol-addressed cross-point cache. The hash covers the geometry CONTENT it is
    evaluated on, charge, mult and the active level, so a protocol change or a different
    conformer geometry addresses a NEW file instead of overwriting the old one. `tag` (the
    partner geometry's state/conformer key) only makes the filename human-readable."""
    h = protocol_hash(_geom_hash(geom_path) or "-", q, m)
    label = tag or geom_path.parent.name
    return DFT / gid / "reorg" / f"{species_state}_at_{label}__{h}.json"


def _cross_energy(gid, species_state, q, m, geom_path, backend, tag=None):
    """SMD (and gas, side-by-side) single point of (q, m) at `geom_path`. Cached.

    Returns the parsed cache dict (with e_smd_eV primary, e_gas_eV diagnostic). Hardened:
      * enforces the requested charge/multiplicity on the cross species,
      * requests SCF stability analysis (stable='opt') and own-density seeding IF the dft_smd
        signature supports them (feature-detected; recorded in the cache as what was actually
        done), and records <S^2> at the cross geometry,
      * keeps an anion-bound diagnostic, now INFORMATIONAL (anions are bound in SMD).
    """
    ghash = _geom_hash(geom_path)
    cache = _cache_path(gid, species_state, geom_path, q, m, tag=tag)
    if cache.exists():
        try:
            d = json.loads(cache.read_text())
            if _cache_ok(d, q, m, ghash):
                return d
        except Exception:
            pass
    import redox.backends.dft as D
    sig = inspect.signature(D.dft_smd).parameters
    extra = {}
    hardening = []
    if "stable" in sig:
        extra["stable"] = "opt"
        hardening.append("stability_analysis")
    if "dm0" in sig or "init_dm" in sig:
        hardening.append("density_seed_requested")   # seed supplied by dft_smd if it accepts one
    # SMD + gas single point on the partner geometry (do_opt=False: this geometry is fixed).
    # do_smd=True gives the solution energy used for lambda_i; do_gas=True keeps the old gas
    # energy side-by-side for the one-transition diagnostic. No geometry optimization, so the
    # charge/mult we pass is the charge/mult that is scored.
    res = D.dft_smd(geom_path, q, m, do_opt=False, do_gas=True, do_smd=True, do_freq=False,
                    backend=backend, **extra)
    if int(res.get("charge")) != int(q) or int(res.get("mult")) != int(m):
        raise RuntimeError(f"{gid}/{species_state}@{geom_path.name}: cross point scored wrong "
                           f"charge/mult ({res.get('charge')}/{res.get('mult')} != {q}/{m})")
    if _geom_hash(geom_path) != ghash:
        raise RuntimeError(f"{gid}/{geom_path}: geometry changed during cross-point evaluation")
    entry = {"cache_schema": CACHE_SCHEMA,
             "e_smd_eV": res.get("e_smd_eV"),
             "e_gas_eV": res.get("e_gas_eV"),
             "smd_homo_eV": res.get("smd_homo_eV"),
             "smd_s_squared": res.get("smd_s_squared"),
             "smd_spin_contam": res.get("smd_spin_contam"),
             "gas_homo_eV": res.get("gas_homo_eV"),
             "gas_s_squared": res.get("gas_s_squared"),
             "anion_unbound": res.get("anion_unbound"),
             "species_state": species_state, "at": tag or geom_path.parent.name,
             "charge": q, "mult": m,
             "geom_sha1": ghash,
             "sp_xc": res.get("sp_xc"), "sp_basis": res.get("sp_basis"),
             "sp_nlc": res.get("sp_nlc"), "sp_disp": res.get("sp_disp"),
             "converged_smd": res.get("converged_smd"),
             "converged_gas": res.get("converged_gas"),
             "scf_hardening": hardening,
             "protocol": ACTIVE_SP}
    cache.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(cache, entry)
    return entry


def _couple_cross_jobs(gid, O, R):
    """Every (species_state, q, m, geom_path, tag) cross point a couple needs, over matched
    conformer families. Used by both compute_group and stale_cross_points so they agree."""
    (sO, qO, mO), (sR, qR, mR) = O, R
    confs_O, _ = _state_conformers(gid, sO)
    confs_R, _ = _state_conformers(gid, sR)
    pairs, _ = match_conformers(confs_O, confs_R)
    jobs = []
    for iO, iR, label in pairs:
        gO, gR = confs_O[iO]["geom"], confs_R[iR]["geom"]
        # E_O at R's conformer geometry, and E_R at O's conformer geometry
        jobs.append((sO, qO, mO, gR, f"{sR}_{label}"))
        jobs.append((sR, qR, mR, gO, f"{sO}_{label}"))
    return jobs


def stale_cross_points(gid):
    """Cross-point cache entries for `gid` that the current protocol will NOT reuse.
    Returns [(path, reason)] — used by --audit to size a regeneration before running it."""
    out = []
    for O, R in _couples(gid):
        for (st, q, m, geom, tag) in _couple_cross_jobs(gid, O, R):
            c = _cache_path(gid, st, geom, q, m, tag=tag)
            if not c.exists():
                out.append((c, "missing")); continue
            try:
                d = json.loads(c.read_text())
            except Exception:
                out.append((c, "unreadable")); continue
            if _cache_ok(d, q, m, _geom_hash(geom)):
                continue
            if int(d.get("cache_schema", 1)) != CACHE_SCHEMA:
                out.append((c, f"schema v{d.get('cache_schema', 1)} (was gas, now SMD)"))
            elif d.get("e_smd_eV") is None:
                out.append((c, "no SMD energy (gas-only entry)"))
            elif any(k not in d for k in REQUIRED_DIAG):
                out.append((c, "missing diagnostics"))
            else:
                out.append((c, "geometry/level mismatch"))
    return out


def _heavy_rmsd(p1, p2):
    """Kabsch heavy-atom RMSD (A) between two opt.xyz geometries; None if unavailable."""
    try:
        import numpy as np
        from ase.io import read
        a1, a2 = read(str(p1)), read(str(p2))
        m = np.array(a1.get_chemical_symbols()) != "H"
        A = a1.positions[m] - a1.positions[m].mean(0)
        B = a2.positions[m] - a2.positions[m].mean(0)
        H = A.T @ B; U, S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
        return float(np.sqrt(((A @ R.T - B) ** 2).sum(1).mean()))
    except Exception:
        return None


def _sym_heavy_rmsd(p1, p2, charge=0):
    """Heavy-atom RMSD minimised over molecular AUTOMORPHISMS (RDKit GetBestRMS).

    A 180 deg rotation of a 2-fold symmetric aryl ring maps the molecule onto itself: the same
    physical conformer, but an index-wise RMSD sees every atom move. Returns None if bond
    perception fails, so callers can fall back to the Kabsch value. (Informational only now:
    the conformer-jump RMSD no longer disqualifies a couple; matched-conformer averaging and
    the self-consistent SMD surface handle different minima directly.)"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdDetermineBonds
        ms = []
        for p in (p1, p2):
            m = Chem.MolFromXYZFile(str(p))
            if m is None:
                return None
            rdDetermineBonds.DetermineConnectivity(m, charge=int(charge))
            ms.append(Chem.RemoveHs(m, sanitize=False))
        return float(AllChem.GetBestRMS(Chem.Mol(ms[1]), Chem.Mol(ms[0])))
    except Exception:
        return None


def compute_group(gid, backend="gpu"):
    for O, R in _couples(gid):
        (sO, qO, mO), (sR, qR, mR) = O, R
        print(f"[reorg] {gid}: couple {sO}(q{qO:+d})/{sR}(q{qR:+d})", flush=True)
        for (st, q, m, geom, tag) in _couple_cross_jobs(gid, O, R):
            if not geom.exists():
                raise RuntimeError(f"{gid}: cross-point geometry missing: {geom}")
            _cross_energy(gid, st, q, m, geom, backend, tag=tag)
    pending = stale_cross_points(gid)
    if pending:
        raise RuntimeError(f"{gid}: {len(pending)} cross points remain missing or invalid")
    print(f"[reorg] {gid} cross points done", flush=True)


def _read_cross(gid, species_state, geom_path, q, m, tag):
    """Reusable SMD cross-point cache entry, or None."""
    c = _cache_path(gid, species_state, geom_path, q, m, tag=tag)
    if not c.exists():
        return None
    try:
        d = json.loads(c.read_text())
    except Exception:
        return None
    return d if _cache_ok(d, q, m, _geom_hash(geom_path)) else None


def lambda_for_couple(gid, O, R):
    (sO, qO, mO), (sR, qR, mR) = O, R
    confs_O, ens_O = _state_conformers(gid, sO)
    confs_R, ens_R = _state_conformers(gid, sR)
    pairs, repaired = match_conformers(confs_O, confs_R)
    conformer_averaged = bool(ens_O and ens_R and len(pairs) > 1)

    per_family = []   # (couple_energy_eV, lam_smd, relax_ox, relax_red, lam_gas, s2_max)
    for iO, iR, label in pairs:
        cO, cR = confs_O[iO], confs_R[iR]
        E_O_at_O = cO["e_smd_eV"]
        E_R_at_R = cR["e_smd_eV"]
        if E_O_at_O is None or E_R_at_R is None:
            return None
        cross_O = _read_cross(gid, sO, cR["geom"], qO, mO, f"{sR}_{label}")   # E_O at geom_R
        cross_R = _read_cross(gid, sR, cO["geom"], qR, mR, f"{sO}_{label}")   # E_R at geom_O
        if cross_O is None or cross_R is None:
            return None
        lam, relax_ox, relax_red = _four_point_lambda(
            E_O_at_O, cross_O["e_smd_eV"], E_R_at_R, cross_R["e_smd_eV"])
        # gas side-by-side (old quantity): only when all four gas energies are finite.
        g_O_at_O, g_R_at_R = cO.get("e_gas_eV"), cR.get("e_gas_eV")
        g_O_at_R, g_R_at_O = cross_O.get("e_gas_eV"), cross_R.get("e_gas_eV")
        gas_vals = [g_O_at_O, g_O_at_R, g_R_at_R, g_R_at_O]
        lam_gas = (_four_point_lambda(g_O_at_O, g_O_at_R, g_R_at_R, g_R_at_O)[0]
                   if all(isinstance(x, (int, float)) and math.isfinite(x) for x in gas_vals)
                   else None)
        s2s = [cross_O.get("smd_s_squared"), cross_R.get("smd_s_squared")]
        s2_max = max((abs(float(s)) for s in s2s if isinstance(s, (int, float))), default=None)
        family_energy = E_O_at_O + E_R_at_R
        per_family.append((family_energy, lam, relax_ox, relax_red, lam_gas, s2_max,
                           cross_O.get("anion_unbound"), cross_R.get("anion_unbound")))

    if not per_family:
        return None
    weights = _boltzmann_weights([f[0] for f in per_family])
    lam = sum(w * f[1] for w, f in zip(weights, per_family))
    relax_ox = sum(w * f[2] for w, f in zip(weights, per_family))
    relax_red = sum(w * f[3] for w, f in zip(weights, per_family))
    gas_terms = [(w, f[4]) for w, f in zip(weights, per_family) if f[4] is not None]
    lam_gas = (sum(w * g for w, g in gas_terms) / sum(w for w, _ in gas_terms)
               if gas_terms else None)
    s2_max = max((f[5] for f in per_family if f[5] is not None), default=None)
    any_anion_unbound = any(f[6] is True or f[7] is True for f in per_family)

    # conformer diagnostic RMSD between the two state minima (informational only)
    gO, gR = DFT / gid / sO / "opt.xyz", DFT / gid / sR / "opt.xyz"
    rmsd_naive = _heavy_rmsd(gO, gR)
    rmsd_sym = _sym_heavy_rmsd(gO, gR, charge=qO)
    rmsd = rmsd_sym if rmsd_sym is not None else rmsd_naive

    own_O, own_R = _res(gid, sO), _res(gid, sR)
    # QC flags: nonnegativity is now STRUCTURAL (each own point is a minimum of its own SMD
    # surface). Flag only REAL failures: a negative beyond SCF noise (-> surface/SCF problem),
    # a spin-contaminated cross SCF, an implausibly large half, or an unresolved spin ensemble.
    reasons = []
    if lam < -0.02:
        reasons.append("negative_lambda_scf")           # structural >=0 violated -> SCF/geometry fault
    if min(relax_ox, relax_red) < -0.02:
        reasons.append("negative_half_scf")
    if s2_max is not None and s2_max > 1.0:
        reasons.append("cross_spin_contaminated")       # |<S^2>| implausible at cross geometry
    if max(relax_ox, relax_red) > 1.0:
        reasons.append("large_half>1eV")
    if any(r.get('spin_ensemble_status') == 'applied' for r in (own_O, own_R) if r):
        reasons.append('spin_ensemble_kinetics_unvalidated')
    flag = "+".join(reasons)

    return dict(id=gid, couple=f"{sO}->{sR}", q_ox=qO, q_red=qR,
                lambda_i_eV=round(lam, 4), lambda_i_meV=round(lam * EV_MEV, 1),
                lambda_i_kJmol=round(lam * EV_KJ, 2),
                lambda_i_smd_eV=round(lam, 4),
                lambda_i_gas_eV=round(lam_gas, 4) if lam_gas is not None else "",
                lambda_i_shift_eV=round(lam - lam_gas, 4) if lam_gas is not None else "",
                relax_ox_meV=round(relax_ox * EV_MEV, 3),
                relax_red_meV=round(relax_red * EV_MEV, 3),
                rmsd_A=round(rmsd, 3) if rmsd is not None else "",
                rmsd_naive_A=round(rmsd_naive, 3) if rmsd_naive is not None else "",
                rmsd_is_sym_aware=bool(rmsd_sym is not None),
                conformer_averaged=bool(conformer_averaged),
                n_conf_matched=len(per_family),
                conformer_repaired=bool(repaired),
                cross_s_squared_max=round(s2_max, 3) if s2_max is not None else "",
                anion_unbound_gas_info=bool(any_anion_unbound),
                flag=flag)


def pd_ok(x):
    """True if x is a finite number (not NaN/None)."""
    try:
        return x is not None and float(x) == float(x) and abs(float(x)) < 1e3
    except Exception:
        return False


def _canon(smi):
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None


def _our_smiles():
    """gid -> SMILES from the library manifest (neutral parent SMILES)."""
    import csv as _csv
    out = {}
    mf = ROOT / "library" / "manifest.csv"
    if mf.exists():
        with mf.open() as f:
            for r in _csv.DictReader(f):
                if r.get("smiles"):
                    out.setdefault(r["id"], r["smiles"])
    return out


def _d3tales_reorg_map():
    """canonical SMILES -> (hole_reorg_eV, electron_reorg_eV) from the D3TaLES dump.
    Returns {} if the (git-ignored) dataset isn't present."""
    import pandas as pd
    p = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, usecols=["smiles", "hole_reorganization_energy",
                                 "electron_reorganization_energy"])
    m = {}
    for _, r in df.iterrows():
        c = _canon(r["smiles"])
        if c and c not in m:
            m[c] = (r["hole_reorganization_energy"], r["electron_reorganization_energy"])
    return m


def aggregate():
    gids = group_ids()
    rows = []
    for gid in gids:
        for O, R in _couples(gid):
            r = lambda_for_couple(gid, O, R)
            if r:
                rows.append(r)
    if not rows:
        print("no lambda values yet (run --only <gid> to compute SMD cross points)")
        _write_table(rows); return

    # --- D3TaLES computed-lambda cross-check: a GAS-PHASE LITERATURE COMPARATOR, not a
    # validation target. D3TaLES uses the same 4-point FORMULA but a GAS surface at a different
    # level; our primary value is now the SMD (solution) surface. The comparison is therefore
    # made against our side-by-side GAS value (lambda_i_gas_eV) on the SAME (gas) surface, and
    # only as a sanity comparator across different functionals/basis. ---
    smi = _our_smiles()
    d3 = _d3tales_reorg_map()
    for r in rows:
        r["d3tales_lambda_eV"] = ""
        r["d3tales_type"] = ""
        c = _canon(smi.get(r["id"], "")) if smi.get(r["id"]) else None
        if c and c in d3:
            hole, elec = d3[c]
            if r["q_ox"] == 1 and r["q_red"] == 0 and pd_ok(hole):
                r["d3tales_lambda_eV"] = round(float(hole), 4); r["d3tales_type"] = "hole"
            elif r["q_ox"] == 0 and r["q_red"] == -1 and pd_ok(elec):
                r["d3tales_lambda_eV"] = round(float(elec), 4); r["d3tales_type"] = "electron"
    hdr = (f"{'id':22s} {'couple':12s} {'lam_smd(eV)':>11s} {'lam_gas(eV)':>11s} "
           f"{'D3TaLES(gas)':>12s} {'type':>8s} {'conf':>5s} {'ok?':>4s}")
    print(hdr); print("-" * len(hdr))
    diffs = []
    for r in sorted(rows, key=lambda x: x["lambda_i_meV"]):
        ok = "yes" if r["lambda_i_eV"] >= 0 else "NEG!"
        lam_gas = r.get("lambda_i_gas_eV", "")
        d3v = r.get("d3tales_lambda_eV", "")
        gas_s = f"{lam_gas:11.3f}" if lam_gas != "" else f"{'-':>11s}"
        if d3v != "" and lam_gas != "":
            diff = abs(lam_gas - d3v); diffs.append(diff)
            d3s = f"{d3v:12.3f}"
        else:
            d3s = f"{'-':>12s}"
        conf = f"{r.get('n_conf_matched', 1):>5d}"
        print(f"{r['id']:22s} {r['couple']:12s} {r['lambda_i_eV']:11.3f} {gas_s} "
              f"{d3s} {r.get('d3tales_type',''):>8s} {conf} {ok:>4s}")
    print("\nlambda_i (primary) = SOLUTION-phase inner-sphere reorganization energy (SMD, 4-point")
    print("Nelsen on one self-consistent surface; nonnegativity is structural). lambda_i_gas is the")
    print("OLD gas-on-SMD-geometry diagnostic, kept for one transition. Solution rates and")
    print("reversibility require separate validation.")
    flagged = [r for r in rows if r.get("flag")]
    if flagged:
        print(f"\n[!] {len(flagged)} couple(s) flagged for a REAL failure (SCF/surface/spin) "
              f"— inspect the cross-point SCF:")
        for r in flagged:
            print(f"    {r['id']:22s} {r['couple']:12s} lam={r['lambda_i_eV']:.3f} "
                  f"relax_ox={r['relax_ox_meV']:.0f}meV relax_red={r['relax_red_meV']:.0f}meV "
                  f"[{r['flag']}]")
    if diffs:
        import statistics
        print(f"\nD3TaLES GAS-phase comparator: n_matched={len(diffs)}  "
              f"MAD(gas vs D3TaLES)={statistics.mean(diffs):.3f} eV (same 4-point FORMULA, "
              f"different LEVEL and both GAS: ours=wB97M-V/def2-TZVPD vs D3TaLES=IP-tuned "
              f"LC-wHPBE/def2-svp — a functional+basis difference, NOT a validation of the "
              f"solution-phase lambda_i).")
    else:
        print("\nD3TaLES GAS-phase comparator: no exact-SMILES matches among current molecules.")
    _write_table(rows)


REORG_COLUMNS = ["id", "couple", "q_ox", "q_red", "lambda_i_eV", "lambda_i_meV",
                 "lambda_i_kJmol", "lambda_i_smd_eV", "lambda_i_gas_eV", "lambda_i_shift_eV",
                 "relax_ox_meV", "relax_red_meV", "rmsd_A",
                 "rmsd_naive_A", "rmsd_is_sym_aware",
                 "conformer_averaged", "n_conf_matched", "conformer_repaired",
                 "cross_s_squared_max", "anion_unbound_gas_info",
                 "flag", "d3tales_lambda_eV", "d3tales_type"]


def _write_table(rows):
    """Always (re)write the table, header-only when no couple has cross points, so a stale
    reorganization.csv from an earlier state of the tree can never survive into the scorecard."""
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "reorganization.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=REORG_COLUMNS); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in REORG_COLUMNS})
    print(f"wrote {out} ({len(rows)} couples)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="compute cross points for one group id")
    ap.add_argument("--all-cross", dest="all_cross", action="store_true",
                    help="compute/refresh cross points for EVERY group (skips entries the "
                         "current cache schema already accepts, so it is cheap to re-run)")
    ap.add_argument("--audit", action="store_true",
                    help="report which cross-point caches the current protocol will not "
                         "reuse, without computing anything")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--backend", default="gpu")
    ap.add_argument("--shard", default=None, help="'n:i' — split --all-cross groups across workers")
    a = ap.parse_args()
    if a.audit:
        total = 0
        for gid in group_ids():
            st = stale_cross_points(gid)
            if st:
                total += len(st)
                reasons = {}
                for _, r in st:
                    reasons[r] = reasons.get(r, 0) + 1
                print(f"  {gid:24s} {len(st):3d} stale  {reasons}")
        print(f"TOTAL stale cross points: {total}")
    elif a.only:
        compute_group(a.only, a.backend)
    elif a.all_cross:
        gids = group_ids()
        if a.shard:
            n, i = (int(x) for x in a.shard.split(":"))
            gids = [g for k, g in enumerate(gids) if k % n == i]
        print(f"[reorg] refreshing cross points for {len(gids)} groups")
        failures = 0
        for gid in gids:
            if not stale_cross_points(gid):
                print(f"[skip] {gid} (cache current)")
                continue
            try:
                compute_group(gid, a.backend)
            except Exception as exc:
                failures += 1
                print(f"[fail] {gid}: {exc}", flush=True)
        if failures:
            raise SystemExit(1)
    elif a.aggregate:
        aggregate()
    else:
        print("use --only <gid> / --all-cross to compute, --audit to inspect, "
              "or --aggregate to assemble")


if __name__ == "__main__":
    main()
