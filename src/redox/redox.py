"""Compute redox potentials from DFT+SMD energies (and UMA gas-phase, for comparison).

For each 1-electron event O + e- -> R (O = higher charge, R = one less):
    dG        = G(R) - G(O)                       [eV]   (electron free energy in ref)
    E_abs     = -dG / n                            [V]    (absolute, vs free electron)
    E_vs_Fc   = E_abs - (SHE_abs + Fc_vs_SHE)      [V]    (referenced to Fc/Fc+)

G = E_smd + G_thermal: the SMD-solvated electronic energy PLUS a GFN2-xTB RRHO thermal
free-energy correction (dft._thermal_correction; g_thermal_eV per state, 0.0 if not computed
for a state). We use xtb RRHO rather than the DFT Hessian because gpu4pyscf's open-shell
(UKS) analytic Hessian is broken for radicals (FINDINGS.md #3-4). Absolute potentials remain
PROVISIONAL until the validation gate (§V) is passed; compared to measurement before any
ranking is trusted.

  python -m redox.redox            # writes results/redox_potentials.csv
"""
from __future__ import annotations
import csv

from redox.common import DFT, RESULTS, UMA, load_config, read_manifest, read_result


def _cfg(name, attr):
    return getattr(load_config(name), attr)


def _energy(root, gid, state, key):
    r = read_result(gid, state, root=root)
    return r.get(key) if r is not None else None


def _fc_reference():
    """Absolute potential of Fc/Fc+ to subtract, preferring the level-matched value.

    Order: (1) our own ferrocene computed at this level (FC_ABS_COMPUTED_V), which
    cancels systematic DFT/SMD error — physics, not fitting; (2) the thermodynamic
    SHE_abs + Fc_vs_SHE constant as a fallback. Returns (value, source_label).
    """
    ref = _cfg("electrolyte", "REFERENCE")
    fc_computed = _cfg("electrolyte", "FC_ABS_COMPUTED_V")
    if fc_computed is not None:
        return float(fc_computed), "Fc-level-matched(computed)"
    return ref["she_abs_V"] + ref["fc_vs_she_V"], "SHE_abs+Fc_vs_SHE(thermo)"


def main():
    e_ref_abs, ref_source = _fc_reference()
    print(f"[ref] Fc/Fc+ absolute reference = {e_ref_abs:.3f} V  ({ref_source})")

    rows = read_manifest()
    # group -> {state: charge}, preserving metadata
    groups = {}
    for r in rows:
        groups.setdefault(r["id"], {"name": r["name"], "family": r["family"], "states": []})
        groups[r["id"]]["states"].append((r["state"], int(r["charge"])))

    out = []
    for gid, g in groups.items():
        states = sorted(g["states"], key=lambda x: -x[1])  # high charge -> low
        for (sO, qO), (sR, qR) in zip(states, states[1:]):
            if qO - qR != 1:
                continue  # not a 1e step
            g_smd_O = _energy(DFT, gid, sO, "e_smd_eV")
            g_smd_R = _energy(DFT, gid, sR, "e_smd_eV")
            # thermal free-energy correction (G = E_smd + G_thermal); 0.0 if not computed
            # yet, which reproduces the prior electronic-energy-only behavior exactly.
            gth_O = _energy(DFT, gid, sO, "g_thermal_eV") or 0.0
            gth_R = _energy(DFT, gid, sR, "g_thermal_eV") or 0.0
            e_uma_O = _energy(UMA, gid, sO, "energy_eV")
            e_uma_R = _energy(UMA, gid, sR, "energy_eV")

            row = dict(id=gid, name=g["name"], family=g["family"],
                       event=f"{sO}->{sR}", q_ox=qO, q_red=qR, ref=ref_source)
            if g_smd_O is not None and g_smd_R is not None:
                dG = (g_smd_R + gth_R) - (g_smd_O + gth_O)
                E_abs = -dG
                row.update(dG_smd_eV=round(dG, 4),
                           g_thermal_ox_eV=round(gth_O, 4),
                           g_thermal_red_eV=round(gth_R, 4),
                           E_abs_V=round(E_abs, 3),
                           E_vs_Fc_V=round(E_abs - e_ref_abs, 3))
            if e_uma_O is not None and e_uma_R is not None:
                dG_gas = e_uma_R - e_uma_O
                row.update(dG_gas_uma_eV=round(dG_gas, 4),
                           E_vs_Fc_gas_V=round(-dG_gas - e_ref_abs, 3))
            out.append(row)

    RESULTS.mkdir(exist_ok=True)
    cols = ["id", "name", "family", "event", "q_ox", "q_red",
            "dG_smd_eV", "g_thermal_ox_eV", "g_thermal_red_eV",
            "E_abs_V", "E_vs_Fc_V", "dG_gas_uma_eV", "E_vs_Fc_gas_V", "ref"]
    outfile = RESULTS / "redox_potentials.csv"
    with outfile.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in out:
            w.writerow({c: r.get(c, "") for c in cols})

    # console summary
    print(f"{'group':16s} {'event':10s} {'E_vs_Fc(SMD)':>13s} {'E_vs_Fc(gas)':>13s}")
    for r in out:
        print(f"{r['id']:16s} {r['event']:10s} "
              f"{str(r.get('E_vs_Fc_V','--')):>13s} {str(r.get('E_vs_Fc_gas_V','--')):>13s}")
    print(f"\nWrote {outfile}  ({len(out)} redox events)")


if __name__ == "__main__":
    main()
