"""Compute the level-matched Fc/Fc+ absolute reference from the ferrocene DFT+SMD run and
write it into config/electrolyte.py (FC_ABS_COMPUTED_V). This is both the reference the
redox scale needs AND a validation check: the value should land near the OROP MeCN
published number for the same functional (B3LYP ~4.63 V) — see electrolyte.FC_ABS_REF_MeCN_V.

  python scripts/set_fc_reference.py            # read calcs/dft/ferrocene, patch electrolyte.py
  python scripts/set_fc_reference.py --dry-run  # print only
"""
import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DFT = ROOT / "calcs" / "dft" / "ferrocene"
ELEC = ROOT / "config" / "electrolyte.py"


def e_smd(state):
    p = DFT / state / "result.json"
    if not p.exists():
        raise SystemExit(f"missing {p} — run: python -m redox.dft --only ferrocene")
    return json.loads(p.read_text())["e_smd_eV"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    e_ox = e_smd("ox")     # Fc+  (higher charge)
    e_neu = e_smd("neu")   # Fc   (reduced)
    # Fc+ + e- -> Fc :  dG = G(Fc) - G(Fc+);  E_abs = -dG = e_smd(ox) - e_smd(neu)
    fc_abs = e_ox - e_neu

    ref_json = json.loads((DFT / "ox" / "result.json").read_text())
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

    txt = ELEC.read_text()
    new = re.sub(r"FC_ABS_COMPUTED_V\s*=\s*None",
                 f"FC_ABS_COMPUTED_V = {fc_abs!r}", txt, count=1)
    if new == txt:
        # already set to a number; replace the numeric value
        new = re.sub(r"FC_ABS_COMPUTED_V\s*=\s*[0-9.eE+-]+",
                     f"FC_ABS_COMPUTED_V = {fc_abs!r}", txt, count=1)
    ELEC.write_text(new)
    print(f"[patched] {ELEC} FC_ABS_COMPUTED_V = {fc_abs:.4f}")


if __name__ == "__main__":
    main()
