#!/usr/bin/env python
"""Positive control for the open-shell SCF guess policy (redox.backends.dft._guess_policy).

The OROP-170 radical anion in SMD(MeCN) has two SCF solutions: minao and atom converge 0.112 eV
above the ground state that only the huckel guess finds (stored 9-guess search:
calcs/validation/scf_level/prod/orop_170/red.json). Re-score that exact geometry through the
PRODUCTION dft_smd and check that the new policy (minao + huckel, escalating to atom when they
disagree)
  * detects the split in SMD (smd_scf_escalated True, spread ~0.112 eV),
  * returns the ground state found by the 9-guess search, in SMD and in gas phase,
  * does not escalate in gas phase, where every guess agrees.

Single points only (no optimization, no Hessian). Writes one JSON record.

  gpu_reserve run <idx> -- env OMP_NUM_THREADS=4 PYTHONPATH=src \\
      python scripts/validation/scf_policy_control.py
"""
from __future__ import annotations

import hashlib
import json
import socket
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.backends import dft  # noqa: E402
from redox.core.common import CALCS, HARTREE_EV  # noqa: E402
from redox.core.io import atomic_write_json  # noqa: E402
from redox.core.protocol import software_versions  # noqa: E402

GEOM = CALCS / "orop" / "170" / "red" / "opt.xyz"
REFERENCE = CALCS / "validation" / "scf_level" / "prod" / "orop_170" / "red.json"
OUT = ROOT / "manuscript" / "evidence" / "2026-11-18_gpu_checks" / "orop170_scf_policy.json"
TOL_HA = 1e-6          # agreement with the stored search (same level, same geometry)


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    ref = json.loads(REFERENCE.read_text())
    assert ref["charge"] == -1 and ref["mult"] == 2
    best = {ph: min(x["e_Ha"] for x in v["guesses"] if x.get("converged"))
            for ph, v in ref["phases"].items()}
    t0 = time.perf_counter()
    res = dft.dft_smd(GEOM, -1, 2, do_opt=False, do_freq=False, backend="gpu")
    wall = time.perf_counter() - t0
    checks = dict(
        smd_escalated=res.get("smd_scf_escalated") is True,
        smd_spread_eV_about_0p112=abs((res.get("smd_scf_guess_spread_eV") or 0) - 0.112) < 0.01,
        smd_ground_state=abs(res["e_smd_Ha"] - best["solv"]) < TOL_HA,
        gas_not_escalated=res.get("gas_scf_escalated") is False,
        gas_ground_state=abs(res["e_gas_Ha"] - best["gas"]) < TOL_HA,
        converged=bool(res.get("converged_smd")) and bool(res.get("converged_gas")),
    )
    record = dict(
        purpose="positive control: open-shell SCF guess policy (minao+huckel, escalate on split)",
        geometry=str(GEOM), geometry_sha256=_sha256(GEOM),
        reference_search=str(REFERENCE), reference_sha256=_sha256(REFERENCE),
        reference_best_Ha=best, tolerance_Ha=TOL_HA,
        computed=dict(e_smd_Ha=res["e_smd_Ha"], e_gas_Ha=res["e_gas_Ha"],
                      scf_guesses=res.get("scf_guesses"),
                      scf_guess_energies_smd=res.get("scf_guess_energies_smd"),
                      scf_guess_energies_gas=res.get("scf_guess_energies_gas"),
                      smd_scf_guess_spread_eV=res.get("smd_scf_guess_spread_eV"),
                      smd_scf_escalated=res.get("smd_scf_escalated"),
                      gas_scf_guess_spread_eV=res.get("gas_scf_guess_spread_eV"),
                      gas_scf_escalated=res.get("gas_scf_escalated"),
                      smd_s_squared=res.get("smd_s_squared"), gas_s_squared=res.get("gas_s_squared"),
                      energy_protocol=res.get("energy_protocol"), timing_s=res.get("timing_s")),
        delta_vs_reference_eV=dict(smd=(res["e_smd_Ha"] - best["solv"]) * HARTREE_EV,
                                   gas=(res["e_gas_Ha"] - best["gas"]) * HARTREE_EV),
        checks=checks, passed=all(checks.values()),
        host=socket.gethostname(), wall_s=round(wall, 1), software=software_versions(),
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(OUT, record)
    print(json.dumps(dict(checks=checks, passed=record["passed"], wall_s=record["wall_s"],
                          delta_vs_reference_eV=record["delta_vs_reference_eV"]), indent=1))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
