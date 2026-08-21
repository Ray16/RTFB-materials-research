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

  PYTHONPATH=src CUDA_VISIBLE_DEVICES=0 python -m redox.reorg --only viologen --backend gpu
  PYTHONPATH=src python -m redox.reorg --aggregate
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
RESULTS = ROOT / "results"
EV_KJ = 96.485
EV_MEV = 1000.0


def _res(gid, state):
    p = DFT / gid / state / "result.json"
    return json.loads(p.read_text()) if p.exists() else None


def _couples(gid):
    """Adjacent-charge couples [(O_state, qO, mO), (R_state, qR, mR)] for a group."""
    states = []
    d = DFT / gid
    if not d.exists():
        return []
    for sd in d.iterdir():
        r = _res(gid, sd.name)
        if r and r.get("e_gas_eV") is not None and (sd / "opt.xyz").exists():
            states.append((sd.name, int(r["charge"]), int(r["mult"])))
    states.sort(key=lambda x: -x[1])
    out = []
    for (sO, qO, mO), (sR, qR, mR) in zip(states, states[1:]):
        if qO - qR == 1:
            out.append(((sO, qO, mO), (sR, qR, mR)))
    return out


def _cross_energy_gas(gid, species_state, q, m, at_geom_state, backend):
    """Gas single-point of (q, m) at the optimized geometry of `at_geom_state`. Cached."""
    cache = DFT / gid / "reorg" / f"{species_state}_at_{at_geom_state}.json"
    if cache.exists():
        try:
            v = json.loads(cache.read_text()).get("e_gas_eV")
            if v is not None:
                return float(v)
        except Exception:
            pass
    import redox.dft as D
    geom = DFT / gid / at_geom_state / "opt.xyz"
    # gas single point only (inner-sphere lambda needs the gas energy; skip the SMD SCF)
    res = D.dft_smd(geom, q, m, do_opt=False, do_gas=True, do_smd=False, do_freq=False,
                    backend=backend)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({"e_gas_eV": res.get("e_gas_eV"),
                                 "species_state": species_state, "at": at_geom_state,
                                 "charge": q, "mult": m}, indent=2))
    return res.get("e_gas_eV")


def compute_group(gid, backend="gpu"):
    for (sO, qO, mO), (sR, qR, mR) in _couples(gid):
        print(f"[reorg] {gid}: couple {sO}(q{qO:+d})/{sR}(q{qR:+d})", flush=True)
        # two cross points (the two own-geometry points are already stored)
        _cross_energy_gas(gid, sO, qO, mO, sR, backend)   # E_O at geom_R
        _cross_energy_gas(gid, sR, qR, mR, sO, backend)   # E_R at geom_O
    print(f"[reorg] {gid} cross points done", flush=True)


def lambda_for_couple(gid, O, R):
    (sO, qO, mO), (sR, qR, mR) = O, R
    E_O_at_O = _res(gid, sO)["e_gas_eV"]
    E_R_at_R = _res(gid, sR)["e_gas_eV"]
    cO = DFT / gid / "reorg" / f"{sO}_at_{sR}.json"
    cR = DFT / gid / "reorg" / f"{sR}_at_{sO}.json"
    if not (cO.exists() and cR.exists()):
        return None
    E_O_at_R = json.loads(cO.read_text())["e_gas_eV"]
    E_R_at_O = json.loads(cR.read_text())["e_gas_eV"]
    lam = (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)     # eV
    return dict(id=gid, couple=f"{sO}->{sR}", q_ox=qO, q_red=qR,
                lambda_i_eV=round(lam, 4), lambda_i_meV=round(lam * EV_MEV, 1),
                lambda_i_kJmol=round(lam * EV_KJ, 2),
                relax_ox_meV=round((E_O_at_R - E_O_at_O) * EV_MEV, 1),
                relax_red_meV=round((E_R_at_O - E_R_at_R) * EV_MEV, 1))


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
    gids = sorted({p.parent.parent.name for p in DFT.glob("*/*/result.json")})
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
    if diffs:
        import statistics
        print(f"\nD3TaLES cross-check: n_matched={len(diffs)}  MAD={statistics.mean(diffs):.3f} eV"
              f"  (same 4-point formula; our geoms are SMD-opt vs D3TaLES gas-opt, hence not exact)")
    else:
        print("\nD3TaLES cross-check: no exact-SMILES matches among current molecules "
              "(most of ours are functionalized; run bare cores for a direct comparison).")
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "reorganization.csv"
    cols = ["id", "couple", "q_ox", "q_red", "lambda_i_eV", "lambda_i_meV",
            "lambda_i_kJmol", "relax_ox_meV", "relax_red_meV",
            "d3tales_lambda_eV", "d3tales_type"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="compute cross points for one group id")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--backend", default="gpu")
    a = ap.parse_args()
    if a.only:
        compute_group(a.only, a.backend)
    elif a.aggregate:
        aggregate()
    else:
        print("use --only <gid> to compute, or --aggregate to assemble")


if __name__ == "__main__":
    main()
