"""Reversibility / structural-integrity screen for redox couples.

A usable flow-battery redox molecule must undergo a REVERSIBLE electron transfer: the charged
state must be (a) electronically BOUND (a real stable ion, not a resonance) and (b) structurally
INTACT (no bond dissociation or ring-opening on charging). Molecules that fail this don't have a
thermodynamic redox potential at all — they undergo dissociative electron attachment (e.g. CCl4 +
e- -> CCl3. + Cl-) or form unbound radical anions. Comparing a computed adiabatic dG to their
irreversible peak potentials is a category error, not a method failure.

This check is therefore three things at once:
  1. a HARD screening filter (reject non-reversible candidates),
  2. a STABILITY signal (a charged state that dissociates/unbinds is unstable),
  3. a gate applied BEFORE accuracy/ranking is assessed.

Per couple O + e- -> R (O = higher charge, R = one less):
  - bound?      : for anion products (q_red < 0), gas-phase EA = E_gas(O) - E_gas(R) must be > 0
                  (the added electron is bound; else the anion is a gas-phase artifact).
  - intact?     : Kabsch-aligned heavy-atom RMSD(O_geom, R_geom) < RMSD_MAX (no dissociation /
                  large rearrangement). Geometries are optimized in independent frames, so we
                  align first (same molecule, atom order preserved -> direct 1:1 Kabsch).
  Verdict: reversible  |  unbound  |  dissociative.

  PYTHONPATH=src python -m redox.reversibility
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
RESULTS = ROOT / "results"

RMSD_MAX = 0.8          # A, heavy-atom; above this = dissociation / big rearrangement
EA_MIN = 0.0            # eV; anion must be gas-phase bound


def _res(gid, state):
    p = DFT / gid / state / "result.json"
    return json.loads(p.read_text()) if p.exists() else None


def kabsch_rmsd_heavy(xyz_a: Path, xyz_b: Path):
    """Heavy-atom RMSD after optimal (Kabsch) superposition. Reported as a diagnostic ONLY —
    NOT used for the verdict, because conformational flexibility (floppy tails, ring pucker,
    methyl rotation) inflates whole-molecule RMSD without any bond breaking. Use bond-graph
    change (below) to detect dissociation instead."""
    from ase.io import read
    a, b = read(str(xyz_a)), read(str(xyz_b))
    sym = np.array(a.get_chemical_symbols())
    if len(a) != len(b) or list(sym) != list(b.get_chemical_symbols()):
        return None
    m = sym != "H"
    P = a.get_positions()[m]; Q = b.get_positions()[m]
    if len(P) < 2:
        return 0.0
    P = P - P.mean(0); Q = Q - Q.mean(0)
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rrot = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(((P @ Rrot.T - Q) ** 2).sum(1).mean()))


def _bond_set(atoms):
    """Covalent bond graph as a set of frozenset({i,j}) pairs (distance < 1.25*(r_i+r_j))."""
    from ase.data import covalent_radii
    pos = atoms.get_positions()
    Z = atoms.get_atomic_numbers()
    n = len(atoms)
    bonds = set()
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(pos[i] - pos[j])
            if d < 1.25 * (covalent_radii[Z[i]] + covalent_radii[Z[j]]):
                bonds.add(frozenset((i, j)))
    return bonds


def connectivity_change(xyz_a: Path, xyz_b: Path):
    """Number of covalent bonds that break or form between the two geometries (0 = intact).
    Robust to conformational change; catches dissociation/ring-opening. None if unusable."""
    from ase.io import read
    a, b = read(str(xyz_a)), read(str(xyz_b))
    if len(a) != len(b) or a.get_chemical_symbols() != b.get_chemical_symbols():
        return None
    ba, bb = _bond_set(a), _bond_set(b)
    return len(ba ^ bb)     # symmetric difference = bonds broken + formed


def _couples(gid):
    states = []
    d = DFT / gid
    if not d.exists():
        return []
    for sd in d.iterdir():
        r = _res(gid, sd.name)
        if r and r.get("e_gas_eV") is not None and (sd / "opt.xyz").exists():
            states.append((sd.name, int(r["charge"])))
    states.sort(key=lambda x: -x[1])
    return [((sO, qO), (sR, qR)) for (sO, qO), (sR, qR) in zip(states, states[1:])
            if qO - qR == 1]


def assess(gid):
    rows = []
    for (sO, qO), (sR, qR) in _couples(gid):
        rO, rR = _res(gid, sO), _res(gid, sR)
        xO, xR = DFT / gid / sO / "opt.xyz", DFT / gid / sR / "opt.xyz"
        # BINDING: use SOLVATED EA (e_smd) — gas-phase dianions are ~always unbound but are
        # perfectly bound and reversible in solution, so gas EA gives false "unbound".
        ea_solv = rO["e_smd_eV"] - rR["e_smd_eV"]           # solvated EA for O + e- -> R
        dbonds = connectivity_change(xO, xR)               # bonds broken/formed (0 = intact)
        rmsd = kabsch_rmsd_heavy(xO, xR)                    # diagnostic only
        unbound = (qR < 0) and (ea_solv <= EA_MIN)         # anion product unbound in solution
        dissoc = (dbonds is None) or (dbonds > 0)          # any bond change = dissociation/ring-open
        verdict = "unbound" if unbound else ("dissociative" if dissoc else "reversible")
        rows.append(dict(id=gid, couple=f"{sO}->{sR}", q_ox=qO, q_red=qR,
                         EA_solv_eV=round(ea_solv, 3), d_bonds=dbonds,
                         rmsd_heavy_A=(round(rmsd, 3) if rmsd is not None else None),
                         reversible=(verdict == "reversible"), verdict=verdict))
    return rows


def main():
    gids = sorted({p.parent.parent.name for p in DFT.glob("*/*/result.json")})
    rows = []
    for gid in gids:
        rows.extend(assess(gid))
    if not rows:
        print("no couples found"); return
    hdr = f"{'id':22s} {'couple':12s} {'EA_solv(eV)':>11s} {'dbonds':>6s} {'RMSD(A)':>8s} {'verdict':>13s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: (x["verdict"] != "reversible", x["id"])):
        rm = f"{r['rmsd_heavy_A']:.3f}" if r["rmsd_heavy_A"] is not None else "n/a"
        db = r["d_bonds"] if r["d_bonds"] is not None else "n/a"
        print(f"{r['id']:22s} {r['couple']:12s} {r['EA_solv_eV']:11.3f} {str(db):>6s} {rm:>8s} "
              f"{r['verdict']:>13s}")
    nrev = sum(r["reversible"] for r in rows)
    print(f"\n{nrev}/{len(rows)} couples reversible. Non-reversible = reject as candidate AND flag "
          f"unstable. Dissociation = bond-graph change (robust to flexibility); RMSD is diagnostic only.")
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "reversibility.csv"
    cols = ["id", "couple", "q_ox", "q_red", "EA_solv_eV", "d_bonds", "rmsd_heavy_A",
            "reversible", "verdict"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
