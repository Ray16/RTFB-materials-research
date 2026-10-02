"""Compute the level-matched Fc/Fc+ absolute reference from the ferrocene DFT+SMD run and
write it into config/project.json (reference.fc_abs_computed_V, exposed as
electrolyte.FC_ABS_COMPUTED_V). This is both the reference the
redox scale needs AND a validation check: the value should land near the OROP MeCN
published number for the same functional (B3LYP ~4.63 V) — see electrolyte.FC_ABS_REF_MeCN_V.

  python scripts/pipeline/set_fc_reference.py            # read calcs/dft/ferrocene, patch project.json
  python scripts/pipeline/set_fc_reference.py --dry-run  # print only

The reference uses the SAME free energy as every molecule (G = E_smd + g_thermal, redox.py);
an electronic-only reference shifts every E_vs_Fc by the Fc/Fc+ thermal difference (-0.08 V).
"""
import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft" / "ferrocene"
PROJECT = ROOT / "config" / "project.json"


def g_solv(state):
    """G = E_smd (ACTIVE protocol) + g_thermal (eV) — identical to redox.py."""
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from redox.core.common import free_energy, read_result
    G, _ = free_energy(read_result("ferrocene", state))
    if G is None:
        raise SystemExit(f"ferrocene/{state} lacks an active-protocol energy or thermal term "
                         f"(run: python -m redox.qm.sp --only ferrocene)")
    return G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    g_ox = g_solv("ox")     # Fc+  (higher charge)
    g_neu = g_solv("neu")   # Fc   (reduced)
    # Fc+ + e- -> Fc :  dG = G(Fc) - G(Fc+);  E_abs = -dG = G(ox) - G(neu)
    fc_abs = g_ox - g_neu

    from redox.core.protocol import ACTIVE_SP
    ref_json = dict(xc=ACTIVE_SP["xc"], basis=ACTIVE_SP["basis"], solvent=ACTIVE_SP["solvent"])
    xc = ref_json.get("xc", "?")
    published = {"b3lyp": 4.63210987710688, "b3lyp-d3": 4.66224854035885,
                "wb97x-d3": 4.618265132412174}.get(xc.lower())
    print(f"Fc/Fc+ absolute reference (computed, {xc}/{ref_json.get('basis','?')}, "
          f"SMD-{ref_json.get('solvent','?')}): {fc_abs:.4f} V")
    if published is not None:
        print(f"OROP published ({xc}, MeCN):        {published:.4f} V  "
              f"(delta {fc_abs - published:+.3f} V)")

    if args.dry_run:
        return

    # in-place edit of the one value (keeps the hand-formatted JSON intact)
    txt = PROJECT.read_text()
    old = json.loads(txt)["reference"].get("fc_abs_computed_V")
    new, n = re.subn(r'("fc_abs_computed_V"\s*:\s*)(null|[0-9.eE+-]+)',
                     lambda m: m.group(1) + repr(fc_abs), txt, count=1)
    if n != 1:
        raise SystemExit(f"could not find reference.fc_abs_computed_V in {PROJECT}")
    PROJECT.write_text(new)
    print(f"[patched] {PROJECT} reference.fc_abs_computed_V: {old} -> {fc_abs:.4f}")

if __name__ == "__main__":
    main()
