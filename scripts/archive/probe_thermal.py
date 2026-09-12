#!/usr/bin/env python
"""DIAGNOSTIC (not pipeline): how big is the thermal free-energy correction for the
viologen couples?

The redox table currently uses G ~ E_elec(SMD) only (thermal corrections deferred).
For the two viologen reductions MV2+ -> MV+. -> MV0 the computed E deg is ~0.35-0.5 V
too negative. Thermal corrections are the roadmap's next term; this probe QUANTIFIES
them before we commit to implementing them, so we don't fix the smallest lever.

Method (cheap, decisive on magnitude): finite-difference Hessian from the UMA MLIP on
each state's relaxed gas-phase geometry (correct charge+spin), -> harmonic vib energies
-> ASE IdealGasThermo Gibbs correction at 298.15 K. Then recompute the two couples with
  G = E_smd(DFT)  +  G_thermal(UMA)
and compare to the electronic-only numbers and to experiment.

This writes nothing into results/ or calcs/ -- it only prints a table.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from redox.uma import make_calculator  # noqa: E402

T = 298.15
GID = "methyl_viologen"
STATES = ["ox2", "ox1", "neu"]          # q = 2, 1, 0
# experimental anchors (V vs Fc/Fc+, MeCN) -- flagged "approx, verify" in config
EXP = {"ox2->ox1": -0.45, "ox1->neu": -0.88}


def _load(state):
    uj = json.loads((ROOT / "calcs" / "uma" / GID / state / "result.json").read_text())
    dj = json.loads((ROOT / "calcs" / "dft" / GID / state / "result.json").read_text())
    return uj, dj


def _thermal_correction(calc, atoms, charge, mult):
    """Return Gibbs thermal correction (eV) = G(T) - E_elec, via UMA harmonic Hessian."""
    from ase.vibrations import Vibrations
    from ase.thermochemistry import IdealGasThermo
    import numpy as np, tempfile, shutil

    a = atoms.copy()
    a.info["charge"] = int(charge)
    a.info["spin"] = int(mult)
    a.calc = calc

    tmp = Path(tempfile.mkdtemp(prefix="vib_"))
    try:
        vib = Vibrations(a, name=str(tmp / "vib"))
        vib.run()
        energies = vib.get_energies()          # complex ndarray (eV)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # IdealGasThermo(nonlinear) wants EXACTLY 3N-6 vib modes. Vibrations returns 3N
    # (complex) energies incl. 6 near-zero trans/rot. Drop the 6 smallest-magnitude
    # modes, then floor tiny/imaginary internal modes to 12 cm^-1 (~0.0015 eV) so the
    # count is preserved (standard low-frequency raising; those modes barely move G).
    FLOOR_eV = 0.0015
    mags = np.array([abs(e) for e in energies])
    keep = np.argsort(mags)[6:]                 # drop 6 trans/rot (lowest |E|)
    n_imag = int(sum(1 for i in keep if abs(energies[i].imag) > 1e-6))
    vib_energies = []
    for i in keep:
        e = energies[i]
        val = e.real if abs(e.imag) < 1e-9 else 0.0
        vib_energies.append(max(val, FLOOR_eV))
    vib_energies = np.sort(np.array(vib_energies))
    assert len(vib_energies) == 3 * len(a) - 6, (len(vib_energies), 3*len(a)-6)

    spin = (int(mult) - 1) / 2.0
    thermo = IdealGasThermo(vib_energies=list(vib_energies), geometry="nonlinear",
                            atoms=a, symmetrynumber=1, spin=spin)
    # potentialenergy defaults to atoms' energy; we want the correction only, so pass 0
    thermo.potentialenergy = 0.0
    g_corr = thermo.get_gibbs_energy(temperature=T, pressure=101325.0, verbose=False)
    return g_corr, len(vib_energies), n_imag


def main():
    from ase.io import read
    print(f"[load] UMA calculator ...", flush=True)
    calc = make_calculator("uma-s-1p2p1", "cuda")

    rows = {}
    for st in STATES:
        uj, dj = _load(st)
        q, mult = int(dj["charge"]), int(dj["mult"])
        geom = ROOT / "calcs" / "uma" / GID / st / "relaxed.xyz"
        atoms = read(str(geom))
        print(f"[freq] {GID}/{st} q={q} m={mult} ({len(atoms)} atoms) ...", flush=True)
        g_corr, nmodes, nimag = _thermal_correction(calc, atoms, q, mult)
        rows[st] = dict(q=q, mult=mult, e_smd_eV=dj["e_smd_eV"],
                        g_corr_eV=g_corr, nmodes=nmodes, nimag=nimag)
        print(f"       G_thermal_corr = {g_corr:+.4f} eV  "
              f"(modes={nmodes}, imag_dropped={nimag})", flush=True)

    print("\n=== per-state ===")
    for st in STATES:
        r = rows[st]
        print(f"  {st:4s} q={r['q']:+d} m={r['mult']}  "
              f"E_smd={r['e_smd_eV']:.3f} eV  Gcorr={r['g_corr_eV']:+.4f} eV")

    print("\n=== couples: E_vs_abs shift from adding thermal (n=1 e-) ===")
    print("  (E_abs = E_O - E_R ; only the DIFFERENCE of Gcorr matters)")
    pairs = [("ox2", "ox1", "ox2->ox1"), ("ox1", "neu", "ox1->neu")]
    for O, R, tag in pairs:
        dE_elec = rows[O]["e_smd_eV"] - rows[R]["e_smd_eV"]          # electronic-only E_abs
        dG_corr = rows[O]["g_corr_eV"] - rows[R]["g_corr_eV"]        # thermal shift to E_abs
        dG_tot = dE_elec + dG_corr
        print(f"  {tag:10s}  E_abs(elec)={dE_elec:+.4f}  "
              f"deltaGcorr={dG_corr:+.4f} eV  E_abs(G)={dG_tot:+.4f}")
        print(f"             -> thermal moves this couple by {dG_corr*1000:+.0f} mV")
    print("\nNOTE: sign/reference to Fc applied in redox.py; here we isolate the "
          "thermal magnitude. If |deltaGcorr| << 0.35 V, thermal is NOT the fix.")


if __name__ == "__main__":
    main()
