"""Compute redox potentials from DFT+SMD energies (and UMA gas-phase, for comparison).

For each 1-electron event O + e- -> R (O = higher charge, R = one less):
    dG        = G(R) - G(O)                       [eV]   (electron free energy in ref)
    E_abs     = -dG / n                            [V]    (absolute, vs free electron)
    E_vs_Fc   = E_abs - (SHE_abs + Fc_vs_SHE)      [V]    (referenced to Fc/Fc+)

G = E_smd + G_thermal: the SMD-solvated electronic energy PLUS a GFN2-xTB RRHO thermal
free-energy correction (dft._thermal_correction; g_thermal_eV per state). A state missing
either term makes the couple INCOMPLETE (status column) — no electronic-only fallback. E_smd is
the ACTIVE-protocol energy (redox.core.protocol: uniform def2-TZVPD) via common.read_result. We use xtb RRHO rather than the DFT Hessian because gpu4pyscf's open-shell
(UKS) analytic Hessian is broken for radicals (FINDINGS.md #3-4). Absolute potentials remain
PROVISIONAL until the validation gate (§V) is passed; compared to measurement before any
ranking is trusted.

  python -m redox.properties.potentials            # writes results/redox_potentials.csv
                                   #      + results/multielectron_stability.csv
"""
from __future__ import annotations
import csv
import math

from redox.core.common import DFT, RESULTS, UMA, free_energy, load_config, read_manifest, read_result

# RT/F at 298.15 K (V): thermal voltage for the comproportionation equilibrium. No
# temperature is stored in project.json, so we fix the standard 298.15 K used everywhere
# else in the pipeline (xtb RRHO thermal corrections are also at 298.15 K).
RT_OVER_F_V = 8.314462618 * 298.15 / 96485.33212  # = 0.025693 V
LOG10_FACTOR_V = RT_OVER_F_V * math.log(10.0)      # = 0.059160 V (Nernst slope)


def _cfg(name, attr):
    return getattr(load_config(name), attr)


def _energy(root, gid, state, key):
    r = read_result(gid, state, root=root)
    return r.get(key) if r is not None else None


def fc_live():
    """Fc+/Fc absolute potential (V) from the ferrocene states at the ACTIVE protocol with the
    SAME free-energy definition as every molecule (G = E_smd + g_thermal). None if incomplete."""
    Gox, _ = free_energy(read_result("ferrocene", "ox"))
    Gneu, _ = free_energy(read_result("ferrocene", "neu"))
    return (Gox - Gneu) if (Gox is not None and Gneu is not None) else None


def _fc_reference():
    """Absolute potential of Fc/Fc+ to subtract.

    Computed LIVE from the ferrocene calculation at the active protocol, so the reference can
    never drift from the level/free-energy definition of the molecules (a stored constant did:
    it was electronic-only while molecules used G, a 79 mV offset). The stored
    FC_ABS_COMPUTED_V is only a cross-check; a mismatch > 5 mV is reported. Refuses to fall
    back to the thermodynamic SHE constant silently: that raises unless explicitly allowed.
    """
    fc = fc_live()
    stored = _cfg("electrolyte", "FC_ABS_COMPUTED_V")
    if fc is not None:
        if stored is not None and abs(stored - fc) > 0.005:
            print(f"[ref] NOTE stored FC_ABS_COMPUTED_V={stored:.4f} differs from live "
                  f"{fc:.4f} V — run scripts/pipeline/set_fc_reference.py to sync the record",
                  flush=True)
        return fc, "Fc-level-matched(live,G)"
    raise SystemExit("ferrocene has no active-protocol free energy (run: python -m redox.qm.sp "
                     "--only ferrocene) — refusing to reference E° to a mismatched constant")


def _multielectron_summary(events, window):
    """Per-molecule multi-electron stability gate from the 1e reduction potentials.

    Pure post-processing of the E_vs_Fc already computed for each 1e step -- NO extra QM.
    For a molecule with >=2 consecutive 1e reductions (O -> I -> R), reports:

      E1, E2   : first / second reduction potential (V vs Fc); E1 is the less-negative step
      dE12     : E1 - E2 (>0 = normal ordering, second electron harder)
      logK_comp: log10 K for comproportionation  O + R <=> 2 I,  = dE12 / (RT/F ln10).
                 K>>1 (dE12>0) -> the 1e intermediate is thermodynamically stable (two
                 resolved waves); K<1 (dE12<0, potential inversion) -> the intermediate
                 disproportionates and both electrons pass in a single 2e wave.
      usable_e : 2 if BOTH potentials fall inside the electrolyte window, else 1 -- this is
                 the capacity gate: a second electron below the cathodic limit is not
                 chargeable, so "2e on paper" collapses to 1e usable capacity.

    Only reductions are chained here (multi-e anolytes); single-event couples (e.g. a lone
    1e oxidation) have no dE12 and are skipped.
    """
    win_lo, win_hi = window
    # group 1e events by molecule, in reduction order (descending oxidized-state charge)
    by_id = {}
    for r in events:
        if r.get("E_vs_Fc_V") is None:
            continue
        by_id.setdefault(r["id"], {"name": r["name"], "family": r["family"], "ev": []})
        by_id[r["id"]]["ev"].append(r)
    rows = []
    for gid, g in by_id.items():
        ev = sorted(g["ev"], key=lambda r: -r["q_ox"])   # most-oxidized step first
        if len(ev) < 2:
            continue
        E1 = ev[0]["E_vs_Fc_V"]          # first reduction (less negative)
        E2 = ev[1]["E_vs_Fc_V"]          # second reduction (more negative, normally)
        dE12 = round(E1 - E2, 3)
        logK = round(dE12 / LOG10_FACTOR_V, 1)
        in1 = win_lo <= E1 <= win_hi
        in2 = win_lo <= E2 <= win_hi
        usable = 2 if (in1 and in2) else (1 if in1 else 0)
        if not in2:
            note = (f"2nd e- at {E2:+.2f} V is outside window [{win_lo},{win_hi}] "
                    f"-> effectively 1e usable capacity")
        elif dE12 >= 0.20:
            note = "two resolved 1e waves; radical-anion intermediate stable (K_comp>>1)"
        elif dE12 <= -0.05:
            note = "potential-inverted single 2e wave; intermediate disproportionates"
        else:
            note = "near-degenerate 2e plateau"
        rows.append(dict(id=gid, name=g["name"], family=g["family"],
                         first_event=ev[0]["event"], second_event=ev[1]["event"],
                         E1_vs_Fc_V=E1, E2_vs_Fc_V=E2, dE12_V=dE12,
                         logK_comp=logK, E1_in_window=in1, E2_in_window=in2,
                         usable_electrons=usable, note=note))
    return rows


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
            # G = E_smd (active protocol) + G_thermal. Missing energy OR thermal term ->
            # the couple is INCOMPLETE (no E° reported); never substitute 0 for thermal.
            rO, rR = read_result(gid, sO), read_result(gid, sR)
            GO, tqO = free_energy(rO)
            GR, tqR = free_energy(rR)
            g_smd_O = rO.get("e_smd_eV") if rO else None
            g_smd_R = rR.get("e_smd_eV") if rR else None
            gth_O = rO.get("g_thermal_eV") if rO else None
            gth_R = rR.get("g_thermal_eV") if rR else None
            e_uma_O = _energy(UMA, gid, sO, "energy_eV")
            e_uma_R = _energy(UMA, gid, sR, "energy_eV")

            row = dict(id=gid, name=g["name"], family=g["family"],
                       event=f"{sO}->{sR}", q_ox=qO, q_red=qR, ref=ref_source)
            miss = [s for s, r_, G_ in ((sO, rO, GO), (sR, rR, GR)) if G_ is None]
            row["status"] = ("ok" if not miss else
                             "INCOMPLETE:" + ",".join(
                                 f"{s}({'no_energy' if (read_result(gid, s) or {}).get('e_smd_eV') is None else 'no_thermal'})"
                                 for s in miss))
            row["thermal_qc"] = ("imag_modes" if "imag_modes" in (tqO, tqR) else
                                 ("ok" if not miss else ""))
            if GO is not None and GR is not None:
                # thermal contribution to E (sensitivity of E° to the RRHO model)
                row["dE_thermal_V"] = round(-(gth_R - gth_O), 4)
                dG = GR - GO
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
            "E_abs_V", "E_vs_Fc_V", "status", "thermal_qc", "dE_thermal_V",
            "dG_gas_uma_eV", "E_vs_Fc_gas_V", "ref"]
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

    # --- multi-electron stability gate (free post-processing of the same potentials) -------
    window = _cfg("electrolyte", "WINDOW_V_VS_FC")
    me = _multielectron_summary(out, window)
    me_cols = ["id", "name", "family", "first_event", "second_event",
               "E1_vs_Fc_V", "E2_vs_Fc_V", "dE12_V", "logK_comp",
               "E1_in_window", "E2_in_window", "usable_electrons", "note"]
    me_file = RESULTS / "multielectron_stability.csv"
    with me_file.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=me_cols); w.writeheader()
        for r in me:
            w.writerow({c: r.get(c, "") for c in me_cols})
    print(f"\nmulti-electron gate (window {tuple(window)} V vs Fc):")
    print(f"{'group':16s} {'E1':>7s} {'E2':>7s} {'dE12':>7s} {'logK':>6s} {'use_e':>6s}  note")
    for r in me:
        print(f"{r['id']:16s} {r['E1_vs_Fc_V']:>7.2f} {r['E2_vs_Fc_V']:>7.2f} "
              f"{r['dE12_V']:>7.2f} {r['logK_comp']:>6.1f} {r['usable_electrons']:>6d}  {r['note']}")
    print(f"\nWrote {me_file}  ({len(me)} multi-electron molecules)")


if __name__ == "__main__":
    main()
