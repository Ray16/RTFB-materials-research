"""The 'free' scorecard axes — computed with no new DFT, to complete the candidate picture.

  - CAPACITY: n accessible 1e- redox events, molecular weight, and specific capacity
      C_spec = n * F / MW  [mAh/g]   (F = 96485 C/mol; 1 mAh = 3.6 C)
    Exact bookkeeping (n from the resolved redox couples; MW from the formula).
  - SYNTHETIC ACCESSIBILITY: RDKit SA_Score from the SMILES (1 = easy ... 10 = hard).
    A cheap heuristic PROXY for synthesizability/cost, not a cost model.
  - SOLUBILITY PROXY: implicit solvation free energy dG_solv = E_SMD - E_gas of the NEUTRAL
    state (the usual solubility-limiting form). More negative = better solvated. RELATIVE
    within a family only — absolute log S needs the solid/lattice term (not computed).

Each axis is labelled by trust level so the Pareto scorecard treats it correctly
(capacity = trusted objective; SA + solubility = proxy filters).

  PYTHONPATH=src python -m redox.capacity_and_proxies
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
RESULTS = ROOT / "results"

F_C_PER_MOL = 96485.0
C_PER_MAH = 3.6

# helper / non-candidate species to skip (ion-pair, cluster, dimer, free-ion references)
SKIP = {"pf6", "mv_ip1", "mv_ip2", "mv_solv", "mv_dimer", "mv_dimer_cfg"}


def _sascorer():
    import sys
    from rdkit.Chem import RDConfig
    sys.path.append(f"{RDConfig.RDContribDir}/SA_Score")
    import sascorer
    return sascorer


def _manifest_smiles():
    out = {}
    mf = ROOT / "library" / "manifest.csv"
    if mf.exists():
        with mf.open() as f:
            for r in csv.DictReader(f):
                if r.get("smiles"):
                    out.setdefault(r["id"], r["smiles"])
    return out


def _states(gid):
    d = DFT / gid
    st = []
    if d.exists():
        for sd in d.iterdir():
            p = sd / "result.json"
            if p.exists():
                try:
                    r = json.loads(p.read_text())
                    st.append((sd.name, int(r["charge"]), r))
                except Exception:
                    pass
    return sorted(st, key=lambda x: -x[1])


def _n_events(states):
    """Number of adjacent-charge 1e- couples (accessible electrons)."""
    qs = sorted({q for _, q, _ in states})
    return sum(1 for a, b in zip(qs, qs[1:]) if b - a == 1)


def _dG_solv_neutral(states):
    """dG_solv of the neutral (charge 0) state; fall back to least-|charge| state."""
    cand = [r for s, q, r in states if q == 0] or \
           [r for s, q, r in sorted(states, key=lambda x: abs(x[1]))]
    for r in cand:
        v = r.get("dG_solv_eV")
        if isinstance(v, (int, float)):
            return float(v)
    return None


def main():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    sascorer = _sascorer()
    smi = _manifest_smiles()

    gids = sorted({p.parent.parent.name for p in DFT.glob("*/*/result.json")})
    rows = []
    for gid in gids:
        if gid in SKIP:
            continue
        states = _states(gid)
        if not states:
            continue
        s = smi.get(gid)
        mol = Chem.MolFromSmiles(s) if s else None
        if mol is None or "." in (s or ""):
            continue   # need a single well-defined neutral molecule
        mw = Descriptors.MolWt(mol)
        n = _n_events(states)
        cspec = n * F_C_PER_MOL / (mw * C_PER_MAH) if mw else None   # mAh/g
        sa = sascorer.calculateScore(mol)
        dgs = _dG_solv_neutral(states)
        rows.append(dict(
            id=gid, n_electrons=n, MW=round(mw, 1),
            specific_capacity_mAh_g=(round(cspec, 1) if cspec else None),
            SA_score=round(sa, 2),
            dGsolv_neutral_eV=(round(dgs, 3) if dgs is not None else None),
        ))

    if not rows:
        print("no candidate molecules found"); return
    hdr = (f"{'id':22s} {'n':>2s} {'MW':>7s} {'Cspec(mAh/g)':>12s} {'SA':>5s} "
           f"{'dGsolv_neu(eV)':>14s}")
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -(x["specific_capacity_mAh_g"] or 0)):
        cs = f"{r['specific_capacity_mAh_g']:.1f}" if r["specific_capacity_mAh_g"] else "n/a"
        dg = f"{r['dGsolv_neutral_eV']:.3f}" if r["dGsolv_neutral_eV"] is not None else "n/a"
        print(f"{r['id']:22s} {r['n_electrons']:2d} {r['MW']:7.1f} {cs:>12s} "
              f"{r['SA_score']:5.2f} {dg:>14s}")
    print("\nTrust: capacity = EXACT (objective). SA = heuristic proxy. dGsolv = RELATIVE "
          "solubility proxy (neutral; absolute logS needs the solid-state term).")
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "capacity_and_proxies.csv"
    cols = ["id", "n_electrons", "MW", "specific_capacity_mAh_g", "SA_score", "dGsolv_neutral_eV"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
