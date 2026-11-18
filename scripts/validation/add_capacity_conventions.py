#!/usr/bin/env python
"""Derive the two capacity conventions in the scored pool CSV (idempotent, CPU only).

    cap_active_mAh_g = n F / (3.6 M_rep)                    primary (manuscript, Figs 5/6)
    cap_salt_mAh_g   = n F / (3.6 (M_rep + n m_ion))        secondary (SI only)

M_rep is the neutral polystyrene repeat-unit mass: the 4-methylbenzyl proxy's CH3 becomes
one backbone carbon per graft site (CH3-Ar -> -[CH2-CH(Ar)]-), i.e. MW + 12.011 * n_graft
(n_graft = 1 for every pool member; checked against properties/capacity_and_proxies.py).
m_ion is PF6- for catholytes (oxidized states need anions) and TBA+ for anolytes.

The salt-loaded convention penalizes anolytes more than catholytes (TBA+ 242 vs PF6- 145
g/mol), so it must not be used to compare capacities across pools. Within a pool it
leaves Pareto-front membership unchanged.

The legacy column cap_mAh_g (= salt-loaded on the capped-monomer MW) is replaced.

    python scripts/validation/add_capacity_conventions.py [csv]
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    ROOT / "results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv"

FARADAY = 96485.33212   # C/mol
C_MASS = 12.011         # g/mol, one backbone carbon per graft site
N_GRAFT = 1             # every pool member carries exactly one 4-methylbenzyl graft
M_PF6 = 144.962         # g/mol, RDKit MolWt of F[P-](F)(F)(F)(F)F
M_TBA = 242.471         # g/mol, RDKit MolWt of tetra-n-butylammonium

COLS = ["id", "name", "family", "pool", "n", "E_mean", "E_couples",
        "cap_active_mAh_g", "cap_salt_mAh_g", "MW", "MW_repeat_unit"]


def main() -> None:
    rows = list(csv.DictReader(open(CSV)))
    for r in rows:
        n = int(r["n"])
        m_rep = float(r["MW"]) + C_MASS * N_GRAFT
        m_ion = M_PF6 if r["pool"] == "catholyte" else M_TBA
        r["MW_repeat_unit"] = f"{m_rep:.1f}"
        r["cap_active_mAh_g"] = f"{n * FARADAY / (3.6 * m_rep):.1f}"
        r["cap_salt_mAh_g"] = f"{n * FARADAY / (3.6 * (m_rep + n * m_ion)):.1f}"
    with open(CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {CSV}")


if __name__ == "__main__":
    main()
