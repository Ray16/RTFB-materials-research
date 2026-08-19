"""Structural (and reorganization-energy) descriptors between redox states.

- RMSD between adjacent redox states (heavy-atom, Kabsch-aligned) and max atom displacement
  -> the "structural change on ox/red" descriptor. Uses UMA-relaxed geometries.
- Inner-sphere reorganization energy lambda_i via the 4-point scheme needs each state's
  energy at BOTH geometries; the two cross single-points are computed with UMA (see
  lambda_inner, which requires a FAIRChem calculator). RMSD needs no new calculations.

  python -m redox.descriptors            # writes results/structure_descriptors.csv (RMSD)
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
UMA = ROOT / "calcs" / "uma"
RESULTS = ROOT / "results"


def _read_xyz(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0])
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0]); xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return syms, np.array(xyz)


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD between P and Q (same ordering) after optimal rotation+translation."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    Pr = Pc @ R
    return float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean()))


def heavy_rmsd(gid: str, sO: str, sR: str):
    fO = UMA / gid / sO / "relaxed.xyz"
    fR = UMA / gid / sR / "relaxed.xyz"
    if not (fO.exists() and fR.exists()):
        return None
    symsO, PO = _read_xyz(fO)
    symsR, PR = _read_xyz(fR)
    heavy = [i for i, s in enumerate(symsO) if s != "H"]
    P, Q = PO[heavy], PR[heavy]
    rmsd = kabsch_rmsd(P, Q)
    # max single-atom displacement after alignment (heavy atoms)
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    maxd = float(np.sqrt(((Pc @ R - Qc) ** 2).sum(1)).max())
    return dict(rmsd_heavy=round(rmsd, 4), max_disp_heavy=round(maxd, 4), n_heavy=len(heavy))


def lambda_inner(calc, gid, sO, sR, qO, mO, qR, mR):
    """4-point inner-sphere reorganization energy (eV). Needs a FAIRChem calculator.
    lambda_i = [E_O(gR) - E_O(gO)] + [E_R(gO) - E_R(gR)]."""
    from ase.io import read
    def sp(geom_state, q, m):
        atoms = read(str(UMA / gid / geom_state / "relaxed.xyz"))
        atoms.info["charge"] = q; atoms.info["spin"] = m
        atoms.calc = calc
        return atoms.get_potential_energy()
    E_O_gO = json.loads((UMA / gid / sO / "result.json").read_text())["energy_eV"]
    E_R_gR = json.loads((UMA / gid / sR / "result.json").read_text())["energy_eV"]
    E_O_gR = sp(sR, qO, mO)   # oxidized charge/spin at reduced geometry
    E_R_gO = sp(sO, qR, mR)   # reduced charge/spin at oxidized geometry
    return (E_O_gR - E_O_gO) + (E_R_gO - E_R_gR)


def main():
    with (ROOT / "library" / "manifest.csv").open() as f:
        rows = list(csv.DictReader(f))
    groups = {}
    for r in rows:
        groups.setdefault(r["id"], []).append((r["state"], int(r["charge"])))

    out = []
    for gid, states in groups.items():
        states = sorted(states, key=lambda x: -x[1])
        for (sO, qO), (sR, qR) in zip(states, states[1:]):
            if qO - qR != 1:
                continue
            d = heavy_rmsd(gid, sO, sR)
            if d:
                out.append(dict(id=gid, event=f"{sO}->{sR}", **d))

    RESULTS.mkdir(exist_ok=True)
    cols = ["id", "event", "rmsd_heavy", "max_disp_heavy", "n_heavy"]
    with (RESULTS / "structure_descriptors.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(out)
    for r in out:
        print(f"{r['id']:16s} {r['event']:10s} RMSD={r['rmsd_heavy']:.3f} Å "
              f"maxdisp={r['max_disp_heavy']:.3f} Å")
    print(f"\nWrote results/structure_descriptors.csv ({len(out)} events)")


if __name__ == "__main__":
    main()
