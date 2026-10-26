#!/usr/bin/env python
"""Place one prospective hybrid-energy pilot with CPU-only dependency probes."""
import inspect
from pathlib import Path
from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


if __name__=='__main__':
    module=api.load_module('prospective_smooth_energy_reserved_placement',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    module.SCRIPT=Path(__file__).resolve()
    module.WORKER=ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'
    body=inspect.getsource(module.preflight)
    old="helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
    if body.count(old)!=1:raise ValueError('CPU dependency helper template changed')
    body=body.replace(old,"helper = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'",1)
    local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(module)),local)
    module.preflight=local['preflight'];module.main()
