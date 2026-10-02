"""Charged-state INTEGRITY screen for redox couples (a necessary-condition gate).

What this checks, per couple O + e- -> R (O = higher charge):
  - bound : for anion products (q_red < 0) the SOLVATED electron attachment
            E_smd(O) - E_smd(R) > 0 (gas EA falsely flags dianions, bound only in solution);
  - intact: no covalent bond broken/formed between the O and R optimized geometries
            (distance-threshold bond graph; robust to torsional flexibility). Kabsch
            heavy-atom RMSD is a diagnostic only.
  Verdict: intact_bound | unbound | bond_change | incomplete (missing active-protocol energy).

What this does NOT establish: electrochemical REVERSIBILITY. Passing is necessary, not
sufficient — it misses follow-up chemistry with unchanged connectivity or with a second
molecule: proton transfer, radical addition, nucleophilic attack, dimerization, low-barrier
rearrangement. "Reversible" is reserved for experimental evidence or validated kinetic /
decomposition models. Use this as a hard pre-filter (a state that dissociates or cannot bind
the electron has no thermodynamic E° at all — e.g. CCl4 dissociative attachment, FINDINGS #7),
never as a stability verdict.

  PYTHONPATH=src python -m redox.properties.integrity      # -> results/state_integrity.csv
"""
from __future__ import annotations
import csv

import numpy as np

from redox.core.common import DFT, RESULTS, state_names, group_ids
from redox.core.common import read_result as _res

EA_MIN = 0.0            # eV; anion must be gas-phase bound


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
    states = [(s, int(_res(gid, s, raw=True)["charge"])) for s in state_names(gid)]
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
        dbonds = connectivity_change(xO, xR)               # bonds broken/formed (0 = intact)
        rmsd = kabsch_rmsd_heavy(xO, xR)                    # diagnostic only
        if rO.get("e_smd_eV") is None or rR.get("e_smd_eV") is None:
            ea_solv, verdict = None, "incomplete"          # no active-protocol energy yet
        else:
            ea_solv = rO["e_smd_eV"] - rR["e_smd_eV"]       # solvated EA for O + e- -> R
            unbound = (qR < 0) and (ea_solv <= EA_MIN)     # anion product unbound in solution
            changed = (dbonds is None) or (dbonds > 0)     # any bond change = dissociation/ring-open
            verdict = "unbound" if unbound else ("bond_change" if changed else "intact_bound")
        rows.append(dict(id=gid, couple=f"{sO}->{sR}", q_ox=qO, q_red=qR,
                         EA_solv_eV=(round(ea_solv, 3) if ea_solv is not None else None),
                         d_bonds=dbonds,
                         rmsd_heavy_A=(round(rmsd, 3) if rmsd is not None else None),
                         intact_bound=(verdict == "intact_bound"), verdict=verdict))
    return rows


def main():
    gids = group_ids()
    rows = []
    for gid in gids:
        rows.extend(assess(gid))
    if not rows:
        print("no couples found"); return
    hdr = f"{'id':22s} {'couple':12s} {'EA_solv(eV)':>11s} {'dbonds':>6s} {'RMSD(A)':>8s} {'verdict':>13s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: (x["verdict"] != "intact_bound", x["id"])):
        rm = f"{r['rmsd_heavy_A']:.3f}" if r["rmsd_heavy_A"] is not None else "n/a"
        db = r["d_bonds"] if r["d_bonds"] is not None else "n/a"
        ea = f"{r['EA_solv_eV']:11.3f}" if r["EA_solv_eV"] is not None else f"{'n/a':>11s}"
        print(f"{r['id']:22s} {r['couple']:12s} {ea} {str(db):>6s} {rm:>8s} "
              f"{r['verdict']:>13s}")
    nok = sum(r["intact_bound"] for r in rows)
    print(f"\n{nok}/{len(rows)} couples intact_bound (necessary condition only — NOT proof of "
          f"electrochemical reversibility). bond_change/unbound = reject; incomplete = no "
          f"active-protocol energy yet.")
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "state_integrity.csv"
    cols = ["id", "couple", "q_ox", "q_red", "EA_solv_eV", "d_bonds", "rmsd_heavy_A",
            "intact_bound", "verdict"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
