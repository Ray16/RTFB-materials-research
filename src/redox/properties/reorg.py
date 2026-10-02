"""Inner-sphere reorganization energy lambda_i (Nelsen 4-point) for each redox couple.

For a couple O + e- -> R (O = higher charge, R = one less), on the two adiabatic surfaces:

    lambda_i = [E_O(q_R) - E_O(q_O)] + [E_R(q_O) - E_R(q_R)]

where E_X(q_Y) is the energy of species X (its own charge+spin) evaluated at the OPTIMIZED
geometry of species Y. Two of the four points are already stored (each state's own energy);
the two CROSS points E_O(q_R) and E_R(q_O) are extra single points at the other geometry.

Why it computes reliably: all four points are the SAME molecule at two geometries, so basis-
set/functional error cancels strongly (same reasoning as the stability Delta-Gs). lambda_i is
by construction >= 0 (both brackets are distortion penalties); a negative value flags a broken
geometry pairing, wrong charge/spin, or atom-order mismatch.

Convention: computed GAS-PHASE (the inner/geometric part; the outer/solvent part is a separate
Marcus-continuum term). Uses the stored gas single-point level (wb97m-v/def2-tzvp) on the
SMD-optimized geometries, so it is consistent with the rest of the pipeline.

Low lambda_i => fast, reversible electron transfer (good). Large geometric reorganization also
tends to correlate with fragility, so it doubles as a stability signal.

  PYTHONPATH=src CUDA_VISIBLE_DEVICES=0 python -m redox.properties.reorg --only viologen --backend gpu
  PYTHONPATH=src python -m redox.properties.reorg --aggregate
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json

from redox.core.common import DFT, EV_KJ, EV_MEV, RESULTS, ROOT, state_names, group_ids
from redox.core.common import read_result as _res
from redox.core.protocol import ACTIVE_SP, protocol_hash

# --- cross-point cache validity -------------------------------------------------------------
# A cache entry used to be accepted on the sole basis of containing `e_gas_eV`. That silently
# pinned every pre-existing entry at whatever protocol wrote it: caches written before the
# gas-HOMO / unbound-anion diagnostics were added have no `anion_unbound` key, so the QC test
# `bool(d.get("anion_unbound"))` evaluated bool(None) -> False and the unbound-anion screen
# never fired on them. They would also never be regenerated, because the e_gas_eV test passed.
#
# A cache entry is now reusable only if it records the same SCHEMA, the same level of theory,
# the same charge/multiplicity, the same geometry (content hash of the xyz it was evaluated
# on), and carries every required diagnostic field. Anything else is recomputed.
CACHE_SCHEMA = 3          # v3: protocol-addressed file name + level-of-theory check
REQUIRED_DIAG = ("gas_homo_eV", "anion_unbound")


def _geom_hash(path):
    """Content hash of a geometry file — ties a cache entry to the exact coordinates used."""
    try:
        return hashlib.sha1(path.read_bytes()).hexdigest()[:16]
    except OSError:
        return None


def _cache_ok(d, q, m, geom_hash, proto=ACTIVE_SP):
    """True if a cached cross-point entry is reusable under the current protocol: same schema,
    level of theory (xc + basis + NLC/dispersion), charge/mult, geometry, diagnostics, and a
    converged SCF."""
    if not isinstance(d, dict) or d.get("e_gas_eV") is None:
        return False
    if (d.get("sp_xc") or "").lower() != proto["xc"] or (d.get("sp_basis") or "").lower() != proto["basis"]:
        return False
    if (d.get("sp_nlc") or None) != proto["nlc"] or (d.get("sp_disp") or None) != proto["disp"]:
        return False
    if d.get("converged_gas") is not True:      # missing flag = unverified = not usable
        return False
    if int(d.get("cache_schema", 1)) != CACHE_SCHEMA:
        return False
    if any(k not in d for k in REQUIRED_DIAG):
        return False
    if int(d.get("charge", -999)) != int(q) or int(d.get("mult", -999)) != int(m):
        return False
    if geom_hash is not None and d.get("geom_sha1") != geom_hash:
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


def _cache_path(gid, species_state, at_geom_state, q, m):
    """Protocol-addressed cross-point cache: the hash covers the geometry it is evaluated on,
    charge, mult and the active level, so a protocol change addresses a NEW file instead of
    overwriting (or silently reusing) the old one."""
    h = protocol_hash(_geom_hash(DFT / gid / at_geom_state / "opt.xyz") or "-", q, m)
    return DFT / gid / "reorg" / f"{species_state}_at_{at_geom_state}__{h}.json"


def _cross_energy_gas(gid, species_state, q, m, at_geom_state, backend):
    """Gas single-point of (q, m) at the optimized geometry of `at_geom_state`. Cached.

    Reuses a cached value only when `_cache_ok` accepts it (schema + level + charge/mult +
    geometry hash + required diagnostics); otherwise recomputes and rewrites the entry.
    """
    geom = DFT / gid / at_geom_state / "opt.xyz"
    ghash = _geom_hash(geom)
    cache = _cache_path(gid, species_state, at_geom_state, q, m)
    if cache.exists():
        try:
            d = json.loads(cache.read_text())
            if _cache_ok(d, q, m, ghash):
                return float(d["e_gas_eV"])
        except Exception:
            pass
    # ADOPT a legacy (un-hashed, schema-2) entry only if it was evaluated on this exact
    # geometry, charge/mult, at the active xc + basis, with the diagnostics present (e.g. the
    # anion cross points were already def2-TZVPD). Copied to the new path; legacy untouched.
    legacy = DFT / gid / "reorg" / f"{species_state}_at_{at_geom_state}.json"
    if legacy.exists():
        try:
            d = json.loads(legacy.read_text())
            if ((d.get("sp_xc") or "").lower() == ACTIVE_SP["xc"]
                    and (d.get("sp_basis") or "").lower() == ACTIVE_SP["basis"]
                    and d.get("geom_sha1") == ghash and d.get("e_gas_eV") is not None
                    and d.get("converged_gas") is True
                    and int(d.get("charge", -999)) == int(q) and int(d.get("mult", -999)) == int(m)
                    and all(k in d for k in REQUIRED_DIAG)):
                d.update(cache_schema=CACHE_SCHEMA, sp_nlc=ACTIVE_SP["nlc"],
                         sp_disp=ACTIVE_SP["disp"], protocol=ACTIVE_SP,
                         adopted_from=legacy.name)
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_text(json.dumps(d, indent=2))
                return float(d["e_gas_eV"])
        except Exception:
            pass
    import redox.qm.dft as D
    # gas single point only (inner-sphere lambda needs the gas energy; skip the SMD SCF)
    res = D.dft_smd(geom, q, m, do_opt=False, do_gas=True, do_smd=False, do_freq=False,
                    backend=backend)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"cache_schema": CACHE_SCHEMA,
                                 "e_gas_eV": res.get("e_gas_eV"),
                                 "gas_homo_eV": res.get("gas_homo_eV"),
                                 "gas_s_squared": res.get("gas_s_squared"),
                                 "anion_unbound": res.get("anion_unbound"),
                                 "species_state": species_state, "at": at_geom_state,
                                 "charge": q, "mult": m,
                                 "geom_sha1": ghash,
                                 "sp_xc": res.get("sp_xc"), "sp_basis": res.get("sp_basis"),
                                 "sp_nlc": res.get("sp_nlc"), "sp_disp": res.get("sp_disp"),
                                 "converged_gas": res.get("converged_gas"),
                                 "protocol": ACTIVE_SP},
                                indent=2))
    return res.get("e_gas_eV")


def stale_cross_points(gid):
    """Cross-point cache entries for `gid` that the current protocol will NOT reuse.
    Returns [(path, reason)] — used by --audit to size a regeneration before running it."""
    out = []
    for (sO, qO, mO), (sR, qR, mR) in _couples(gid):
        for (st, q, m, at) in ((sO, qO, mO, sR), (sR, qR, mR, sO)):
            c = _cache_path(gid, st, at, q, m)
            if not c.exists():
                out.append((c, "missing")); continue
            try:
                d = json.loads(c.read_text())
            except Exception:
                out.append((c, "unreadable")); continue
            if _cache_ok(d, q, m, _geom_hash(DFT / gid / at / "opt.xyz")):
                continue
            if int(d.get("cache_schema", 1)) != CACHE_SCHEMA:
                out.append((c, f"schema v{d.get('cache_schema', 1)}"))
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

    Why this and not the plain Kabsch value: a 180 deg rotation of a 2-fold symmetric aryl
    ring (para-phenylene, the two halves of a viologen, a symmetric diimide) maps the molecule
    onto ITSELF. It is the same physical conformer, but an index-wise RMSD sees every atom
    move and reports a large number. Measured here: ndi_ammonium ox->red1 falls 2.622 -> 0.123
    A, mophquinone neu->red1 1.010 -> 0.063 A, methyl_viologen ox1->neu 1.281 -> 0.059 A.
    Using the naive value made `conformer_jump` fire on couples whose lambda was perfectly
    fine. Returns None if bond perception fails, so callers can fall back."""
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
    for (sO, qO, mO), (sR, qR, mR) in _couples(gid):
        print(f"[reorg] {gid}: couple {sO}(q{qO:+d})/{sR}(q{qR:+d})", flush=True)
        # two cross points (the two own-geometry points are already stored)
        _cross_energy_gas(gid, sO, qO, mO, sR, backend)   # E_O at geom_R
        _cross_energy_gas(gid, sR, qR, mR, sO, backend)   # E_R at geom_O
    print(f"[reorg] {gid} cross points done", flush=True)


def lambda_for_couple(gid, O, R):
    (sO, qO, mO), (sR, qR, mR) = O, R
    # own-geometry points: ACTIVE-protocol records; cross points: protocol-addressed caches.
    # All four at one level -> each half is a single-basis surface difference.
    E_O_at_O = _res(gid, sO)["e_gas_eV"]
    E_R_at_R = _res(gid, sR)["e_gas_eV"]
    cO = _cache_path(gid, sO, sR, qO, mO)
    cR = _cache_path(gid, sR, sO, qR, mR)
    if E_O_at_O is None or E_R_at_R is None or not (cO.exists() and cR.exists()):
        return None
    for c, q_, m_, at in ((cO, qO, mO, sR), (cR, qR, mR, sO)):
        if not _cache_ok(json.loads(c.read_text()), q_, m_, _geom_hash(DFT / gid / at / "opt.xyz")):
            return None
    E_O_at_R = json.loads(cO.read_text())["e_gas_eV"]
    E_R_at_O = json.loads(cR.read_text())["e_gas_eV"]
    lam = (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)     # eV
    relax_ox = (E_O_at_R - E_O_at_O)                        # eV, each half-relaxation
    relax_red = (E_R_at_O - E_R_at_R)
    # QC guards (see the D3TaLES thiosuccinimide diagnosis): a valid inner-sphere lambda is >=0,
    # each half is a small distortion penalty (~0.1-0.7 eV), the two states differ by only local
    # inner-sphere modes (small heavy-atom RMSD), and a reduced state's gas anion must be BOUND.
    #   - large RMSD  -> the two states optimized into different CONFORMERS (lambda contaminated)
    #   - anion_unbound (gas HOMO>0) -> the gas cross-point energy is meaningless (score in SMD)
    gO, gR = DFT / gid / sO / "opt.xyz", DFT / gid / sR / "opt.xyz"
    rmsd_naive = _heavy_rmsd(gO, gR)
    # Prefer the automorphism-aware value; fall back to plain Kabsch if perception fails.
    rmsd_sym = _sym_heavy_rmsd(gO, gR, charge=qO)
    rmsd = rmsd_sym if rmsd_sym is not None else rmsd_naive
    anion_homo = json.loads(cR.read_text()).get("gas_homo_eV")   # E_R at neutral geom (the anion)
    anion_unbound = bool(json.loads(cR.read_text()).get("anion_unbound")) if qR < 0 else False
    # Record EVERY failed check ('+'-joined): a single precedence-ordered label let
    # conformer_jump hide anion_unbound, so a couple looked rescuable by a conformer-matched
    # recompute when its gas cross point is ill-defined regardless of conformer.
    reasons = []
    if lam < 0:
        reasons.append("negative_lambda")
    if min(relax_ox, relax_red) < -0.02:
        reasons.append("negative_half")
    if rmsd is None:
        reasons.append("rmsd_unavailable")      # conformer check could not run -> not QC-clean
    elif rmsd > 0.40:
        reasons.append("conformer_jump")
    if anion_unbound:
        reasons.append("anion_unbound")
    if max(relax_ox, relax_red) > 1.0:
        reasons.append("large_half>1eV")
    flag = "+".join(reasons)
    return dict(id=gid, couple=f"{sO}->{sR}", q_ox=qO, q_red=qR,
                lambda_i_eV=round(lam, 4), lambda_i_meV=round(lam * EV_MEV, 1),
                lambda_i_kJmol=round(lam * EV_KJ, 2),
                relax_ox_meV=round(relax_ox * EV_MEV, 3),
                relax_red_meV=round(relax_red * EV_MEV, 3),
                rmsd_A=round(rmsd, 3) if rmsd is not None else "",
                rmsd_naive_A=round(rmsd_naive, 3) if rmsd_naive is not None else "",
                rmsd_is_sym_aware=bool(rmsd_sym is not None),
                anion_homo_eV=round(anion_homo, 3) if anion_homo is not None else "",
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
        print("no lambda values yet (run --only <gid> to compute cross points)"); return

    # --- D3TaLES computed-lambda cross-check (identical 4-point definition, different level) ---
    smi = _our_smiles()
    d3 = _d3tales_reorg_map()
    for r in rows:
        r["d3tales_lambda_eV"] = ""
        r["d3tales_type"] = ""
        c = _canon(smi.get(r["id"], "")) if smi.get(r["id"]) else None
        if c and c in d3:
            hole, elec = d3[c]
            # hole reorg = neutral<->cation (q_ox=+1,q_red=0); electron = neutral<->anion (0,-1)
            if r["q_ox"] == 1 and r["q_red"] == 0 and pd_ok(hole):
                r["d3tales_lambda_eV"] = round(float(hole), 4); r["d3tales_type"] = "hole"
            elif r["q_ox"] == 0 and r["q_red"] == -1 and pd_ok(elec):
                r["d3tales_lambda_eV"] = round(float(elec), 4); r["d3tales_type"] = "electron"
    hdr = (f"{'id':22s} {'couple':12s} {'lam_i(eV)':>9s} {'D3TaLES(eV)':>11s} "
           f"{'type':>8s} {'|diff|':>7s} {'ok?':>4s}")
    print(hdr); print("-" * len(hdr))
    diffs = []
    for r in sorted(rows, key=lambda x: x["lambda_i_meV"]):
        ok = "yes" if r["lambda_i_eV"] >= 0 else "NEG!"
        d3 = r.get("d3tales_lambda_eV", "")
        if d3 != "":
            diff = abs(r["lambda_i_eV"] - d3); diffs.append(diff)
            d3s, dfs = f"{d3:11.3f}", f"{diff:7.3f}"
        else:
            d3s, dfs = f"{'-':>11s}", f"{'-':>7s}"
        print(f"{r['id']:22s} {r['couple']:12s} {r['lambda_i_eV']:9.3f} {d3s} "
              f"{r.get('d3tales_type',''):>8s} {dfs} {ok:>4s}")
    print("\nlambda_i = inner-sphere reorganization energy (gas, 4-point Nelsen; identical to")
    print("D3TaLES's ReorganizationCalc). Lower = faster/more reversible ET. Must be >= 0.")
    flagged = [r for r in rows if r.get("flag")]
    if flagged:
        print(f"\n[!] {len(flagged)} couple(s) flagged as likely broken/wrong-conformer geometry "
              f"— re-run with `python -m redox.qm.dft --only <id> --force --torsion-scan`:")
        for r in flagged:
            print(f"    {r['id']:22s} {r['couple']:12s} lam={r['lambda_i_eV']:.3f} "
                  f"relax_ox={r['relax_ox_meV']:.0f}meV relax_red={r['relax_red_meV']:.0f}meV "
                  f"[{r['flag']}]")
    if diffs:
        import statistics
        print(f"\nD3TaLES cross-check: n_matched={len(diffs)}  MAD={statistics.mean(diffs):.3f} eV"
              f"  (same 4-point FORMULA, but different LEVEL: ours=wB97M-V/def2-TZVPD gas vs "
              f"D3TaLES=IP-tuned LC-wHPBE/def2-svp gas — a functional+basis difference, see "
              f"FINDINGS #10 & validate_reorg_worker_d3tales.py for the level-matched comparison)")
    else:
        print("\nD3TaLES cross-check: no exact-SMILES matches among current molecules "
              "(most of ours are functionalized; run bare cores for a direct comparison).")
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "reorganization.csv"
    cols = ["id", "couple", "q_ox", "q_red", "lambda_i_eV", "lambda_i_meV",
            "lambda_i_kJmol", "relax_ox_meV", "relax_red_meV", "rmsd_A",
            "rmsd_naive_A", "rmsd_is_sym_aware",
            "anion_homo_eV", "flag", "d3tales_lambda_eV", "d3tales_type"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")


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
        for gid in gids:
            if not stale_cross_points(gid):
                print(f"[skip] {gid} (cache current)")
                continue
            try:
                compute_group(gid, a.backend)
            except Exception as exc:
                print(f"[fail] {gid}: {exc}", flush=True)
    elif a.aggregate:
        aggregate()
    else:
        print("use --only <gid> / --all-cross to compute, --audit to inspect, "
              "or --aggregate to assemble")


if __name__ == "__main__":
    main()
