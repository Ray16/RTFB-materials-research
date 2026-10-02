#!/usr/bin/env python
"""Re-run ONE phase of a candidate SCF check in a fresh process and merge it into the record
(calcs/validation/scf_level/candidates/<id>/<name>.json). For large open-shell systems the
9-guess sequence can fragment GPU memory so the second phase runs out of memory although
the same SCF fits when run first (as production does).

  python scripts/validation/reference/scf_check_phase.py ndi_ammonium red1 solv --backend gpu
"""
import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scf_level_check as S  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("id"); ap.add_argument("name"); ap.add_argument("phase", choices=["gas", "solv"])
ap.add_argument("--backend", default="gpu")
a = ap.parse_args()
p = S.CAND_OUT / a.id / f"{a.name}.json"
rec = json.loads(p.read_text())
t0 = time.time()
ph = S._multi_guess(Path(rec["xyz"]), rec["q"], rec["mult"], S.LEVELS["prod"], (a.phase,),
                    a.backend, tag=f"cand {a.id}/{a.name}")[a.phase]
pe = rec["prod"].get(a.phase)
ph["production_minus_best_eV"] = (pe - ph["best"]["e_Ha"] * S.HARTREE_EV
                                  if (ph["best"] and pe is not None) else None)
ph["rerun_wall_s"] = round(time.time() - t0, 1)
rec["phases"][a.phase] = ph
p.write_text(json.dumps(rec, indent=1))
print("merged", p, "best", ph["best"])
