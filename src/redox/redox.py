"""Compute redox potentials from DFT+SMD energies (and UMA gas-phase, for comparison).

For each 1-electron event O + e- -> R (O = higher charge, R = one less):
    dG        = G(R) - G(O)                       [eV]   (electron free energy in ref)
    E_abs     = -dG / n                            [V]    (absolute, vs free electron)
    E_vs_Fc   = E_abs - (SHE_abs + Fc_vs_SHE)      [V]    (referenced to Fc/Fc+)

G is approximated by the SMD-solvated electronic energy (thermal corrections deferred —
see docs/PLAN.md). This is PROVISIONAL until the validation gate (§V) is passed; potentials
are compared to measurement before any ranking is trusted.

  python -m redox.redox            # writes results/redox_potentials.csv
"""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
UMA = ROOT / "calcs" / "uma"
RESULTS = ROOT / "results"


def _cfg(name, attr):
    spec = importlib.util.spec_from_file_location(name, ROOT / "config" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return getattr(mod, attr)


def read_manifest():
    with (ROOT / "library" / "manifest.csv").open() as f:
        return list(csv.DictReader(f))


def _energy(root: Path, gid, state, key):
    p = root / gid / state / "result.json"
    if not p.exists():
        return None
    return json.loads(p.read_text()).get(key)


def main():
    ref = _cfg("electrolyte", "REFERENCE")
    e_ref_abs = ref["she_abs_V"] + ref["fc_vs_she_V"]   # abs potential of Fc/Fc+

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
            e_uma_O = _energy(UMA, gid, sO, "energy_eV")
            e_uma_R = _energy(UMA, gid, sR, "energy_eV")

            row = dict(id=gid, name=g["name"], family=g["family"],
                       event=f"{sO}->{sR}", q_ox=qO, q_red=qR)
            if g_smd_O is not None and g_smd_R is not None:
                dG = g_smd_R - g_smd_O
                E_abs = -dG
                row.update(dG_smd_eV=round(dG, 4),
                           E_abs_V=round(E_abs, 3),
                           E_vs_Fc_V=round(E_abs - e_ref_abs, 3))
            if e_uma_O is not None and e_uma_R is not None:
                dG_gas = e_uma_R - e_uma_O
                row.update(dG_gas_uma_eV=round(dG_gas, 4),
                           E_vs_Fc_gas_V=round(-dG_gas - e_ref_abs, 3))
            out.append(row)

    RESULTS.mkdir(exist_ok=True)
    cols = ["id", "name", "family", "event", "q_ox", "q_red",
            "dG_smd_eV", "E_abs_V", "E_vs_Fc_V", "dG_gas_uma_eV", "E_vs_Fc_gas_V"]
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
