#!/usr/bin/env python
"""Reuse the fixed-density response diagnostic with the controlled R profile."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


ADAPTER = ROOT / "scripts/pipeline/smd_first_response_workspace_v2.py"


def checked_controls(folder, budget=28):
    if budget != 28:
        raise ValueError("This profile requires its 28 GiB controls")
    original = ROOT / "scripts/validation/qualify_pcm_first_response.py"
    base = api.load_module("first_response_control_checks", original)
    base.ADAPTER = ADAPTER
    proof = base.checked_controls(folder)
    for mult in (1,2):
        path = folder / "controls" / str(mult) / "compact/execution/contracted_R.json"
        data = json.loads(path.read_text())
        if (data["R_applications"] < 1 or data["equations_changed"] is not False
                or any(v < 1 for v in data["transformed_calls"].values())):
            raise ValueError("Every contracted Hessian response path must be exercised")
        api.verify_sources(data["source_sha256"])
        proof[str(path)] = api.digest(path)
    return proof


def main():
    original = ROOT / "scripts/validation/qualify_pcm_first_response.py"
    module = api.load_module("preserved_fixed_density_response_diagnostic", original)
    module.ADAPTER = ADAPTER
    module.checked_controls = checked_controls
    module.main()
    if sys.argv[1] == "prepare":
        out = Path(sys.argv[sys.argv.index("--out")+1]).resolve()
        plan = api.read_plan(out)
        for p in (Path(__file__).resolve(), original):
            plan["producer_sha256"][str(p)] = api.digest(p)
        plan["execution_profile"] = "Blocked first derivatives and contracted Hessian R; fixed-density response pilot only"
        atomic_write_json(out / "plan.json", plan)


if __name__ == "__main__":
    main()
