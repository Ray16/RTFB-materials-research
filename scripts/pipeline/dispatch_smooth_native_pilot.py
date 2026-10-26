#!/usr/bin/env python
"""Place one alternate-native pilot with existing reserved runtime guards."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    source=ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py'
    module=api.load_module('smooth_native_pilot_resource_contract',source)
    module.SCRIPT=Path(__file__).resolve()
    module.WORKER=ROOT/'scripts/pipeline/smooth_polymer_native.py'
    body=inspect.getsource(module.preflight)
    old="helper = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'"
    if body.count(old)!=1:raise ValueError('CPU-safe node preflight template changed')
    body=body.replace(old,"helper = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'",1)
    local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(module)),local)
    module.preflight=local['preflight']
    body=inspect.getsource(module.wait)
    changes={"manifest['stage'] != 'electronic'":"manifest['stage'] != 'native'",
        'experimental_energy/status.json':'experimental_native/status.json',
        'experimental_electronic_scientific_stage_started':'experimental_native_scientific_stage_started'}
    for old,new in changes.items():
        if body.count(old)!=1:raise ValueError('Native pilot placement template changed')
        body=body.replace(old,new,1)
    local={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(module)),local)
    module.wait=local['wait']
    module.main()


if __name__=='__main__':main()
