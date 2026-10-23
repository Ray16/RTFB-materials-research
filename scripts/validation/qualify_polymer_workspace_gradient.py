#!/usr/bin/env python
"""Apply the separately qualified gradient contractions to a full resource pilot.

Reuse the existing tensor inventory, fixed-geometry energy/gradient evaluation
and continuity gates. Only the explicitly controlled execution adapter changes.
"""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    original = ROOT / "scripts/validation/qualify_polymer_workspace.py"
    module = api.load_module("preserved_full_workspace_qualifier", original)
    module.ADAPTER = ROOT / "scripts/pipeline/smd_gradient_workspace.py"
    previous = module.checked_controls

    def checked(paths, budget):
        sources = previous(paths, budget)
        for value in paths:
            p = Path(value).parent / "compact/execution/status.json"
            report = json.loads(p.read_text())
            if (report.get("gradient_solver_calls", 0) < 1
                    or report.get("direct_LU_solver_retained") is not True
                    or report.get("gradient_contraction_branch_only") is not True):
                raise ValueError("Direct-LU bounded-gradient control was not exercised")
        return sources

    module.checked_controls = checked
    module.main()
    if sys.argv[1] == "prepare":
        campaign = Path(sys.argv[sys.argv.index("--campaign") + 1]).resolve()
        plan = api.read_plan(campaign)
        for p in (Path(__file__).resolve(), original):
            plan["producer_sha256"][str(p)] = api.digest(p)
        plan["execution_profile_qualifier"] = str(Path(__file__).resolve())
        plan["bounded_gradient_controls_required"] = True
        atomic_write_json(campaign / "plan.json", plan)


if __name__ == "__main__":
    main()
