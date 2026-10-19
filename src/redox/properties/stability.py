"""Disproportionation thermodynamics of redox intermediates.

Selecting a good flow-battery redox molecule needs BOTH a target potential AND a durable
charged/radical state. This module checks one balanced thermodynamic pathway;
other decomposition reactions, barriers and cycling conditions remain separate.

Reactions scored (per reactive intermediate R^q, the state that actually cycles):

  (1) DISPROPORTIONATION   2 R^q -> R^(q+1) + R^(q-1)
      dG_disp = [G(R^q+1) + G(R^q-1)] - 2 G(R^q)
      dG_disp > 0  => disproportionation is thermodynamically uphill for R^q.
      This is an exact arithmetic identity for the selected redox-state energies,
      with atom and electron conservation. Same-protocol errors can partly cancel;
      physical accuracy and cycling lifetime require separate validation. Equivalent wave
      spacing: dG_disp = F*(E_high - E_low); K_comp = exp(+dG_disp/RT).

  (2) DIMERIZATION         2 R^q -> (R2)^(2q)
      dG_dim = G(dimer) - 2 G(R^q). Not evaluated by this module. A qualified
      exact-molecule dimer, matched solution treatment and association standard state
      are required. An uphill association alone does not prove cycling durability.

Sign convention everywhere: dG_decomp > 0 == STABLE (decomposition uphill).
Free energy per species: G = e_smd_eV + the admitted thermal-policy correction.
Native solution DFT thermals and historical GFN2-xTB estimates remain separate.

  PYTHONPATH=src python -m redox.properties.stability            # disproportionation table
"""
from __future__ import annotations
import csv
import json
import math

from redox.core.common import (DFT, EV_KJ, KT_EV, RESULTS, free_energy, group_ids, read_result,
                               state_names)

F_EV_PER_V = 1.0        # 1 electron: dG(eV) per volt is numerically 1


def G(gid, state):
    """Free energy (eV) = e_smd (active protocol) + g_thermal; None (INCOMPLETE) if either
    is missing — a missing thermal term is never replaced by 0."""
    return free_energy(read_result(gid, state))[0]


def _states_by_charge(gid):
    """Return [(state, charge)] present in calcs/dft for gid, sorted high->low charge,
    using the charge recorded in each result.json (authoritative)."""
    out = []
    for st in state_names(gid):
        q = read_result(gid, st, raw=True).get("charge")
        if q is not None:
            out.append((st, int(q)))
    return sorted(out, key=lambda x: -x[1])


def disproportionation(gid):
    """dG_disp for every interior state R^q that has BOTH neighbors R^(q+1), R^(q-1).

    Returns list of dicts. The interior state is the reactive intermediate (radical /
    semiquinone / radical cation) whose persistence we care about.
    """
    states = _states_by_charge(gid)          # high -> low charge, adjacent differ by 1e-
    by_q = {q: s for s, q in states}
    rows = []
    for s, q in states:
        if (q + 1) in by_q and (q - 1) in by_q:
            g_mid = G(gid, s)
            g_hi = G(gid, by_q[q + 1])
            g_lo = G(gid, by_q[q - 1])
            if None in (g_mid, g_hi, g_lo):
                print(f"[incomplete] {gid}/{s}: missing active-protocol energy or thermal "
                      f"term -> dG_disp not reported", flush=True)
                continue
            dG = (g_hi + g_lo) - 2.0 * g_mid          # eV
            records = [read_result(gid,st) for st in (by_q[q+1],s,by_q[q-1])]
            qc=[free_energy(record)[1] for record in records]
            definitions = {record.get('thermal_definition', 'legacy_gfn2_gas_rrho') for record in records}
            valid=all(value=='ok' for value in qc) and len(definitions)==1
            thermal_qc = 'thermal_definition_mismatch' if len(definitions)!=1 else (
                'ok' if valid else next(value for value in qc if value!='ok'))
            row=dict(
                id=gid, intermediate=s, q=q,
                reaction=f"2 {s}(q={q:+d}) -> {by_q[q+1]}(q={q+1:+d}) + {by_q[q-1]}(q={q-1:+d})",
                status='ok' if valid else 'THERMAL_UNVALIDATED',thermal_qc=thermal_qc)
            values=dict(dG_disp_eV=round(dG,4),dG_disp_kJmol=round(dG*EV_KJ,1),
                        logK_disp=round(-dG/(KT_EV*math.log(10)),2),stable_vs_disprop=bool(dG>0))
            row.update(values if valid else {key+'_estimate':value for key,value in values.items()})
            rows.append(row)
    return rows


def main():
    # every group in calcs/dft that has >=3 charge states can have an interior intermediate
    gids = group_ids()
    all_rows = []
    for gid in gids:
        all_rows.extend(disproportionation(gid))

    if not all_rows:
        print("no interior redox states found (need 3 consecutive charge states)")
        return

    hdr = f"{'id':22s} {'interm':6s} {'q':>3s} {'dG_disp(kJ/mol)':>15s} {'log10K':>7s} {'stable?':>8s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(all_rows, key=lambda x: -(x.get("dG_disp_kJmol") or x.get('dG_disp_kJmol_estimate') or 0)):
        if r['status']!='ok':
            print(f"{r['id']:22s} {r['intermediate']:6s} THERMAL_UNVALIDATED (estimate retained)");continue
        print(f"{r['id']:22s} {r['intermediate']:6s} {r['q']:+3d} "
              f"{r['dG_disp_kJmol']:15.1f} {r['logK_disp']:7.2f} "
              f"{'YES' if r['stable_vs_disprop'] else 'no':>8s}")
    print("\nSign convention: dG_disp > 0  =>  intermediate STABLE against disproportionation.")
    print("(Exact arithmetic for admitted redox energies; decomposition pathways and cycling lifetime remain separate.)")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "stability_disproportionation.csv"
    cols = ["id", "intermediate", "q", "reaction", "dG_disp_eV", "dG_disp_kJmol",
            "logK_disp", "stable_vs_disprop","status","thermal_qc",
            "dG_disp_eV_estimate","dG_disp_kJmol_estimate","logK_disp_estimate","stable_vs_disprop_estimate"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
