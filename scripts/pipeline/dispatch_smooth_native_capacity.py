#!/usr/bin/env python
"""Place one native pilot with bound-aware CPU preflight and reserved execution."""
import inspect
from pathlib import Path

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def main():
    helper=api.load_module('single_native_capacity_contract',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    native=api.load_module('single_native_capacity_placement',ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py')
    helper.SCRIPT=Path(__file__).resolve();helper.WORKER=ROOT/'scripts/pipeline/smooth_polymer_native.py'
    helper.placement_module=native.native_placement
    # Keep the unmodified electronic preflight available to the native adapter.
    original=api.load_module('original_single_native_preflight',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py')
    helper.preflight=lambda batch,host,profile:native.cpu_preflight(original,batch,host,profile)
    body=inspect.getsource(helper.wait)
    for old,new in {"manifest['stage'] != 'electronic'":"manifest['stage'] != 'native'",
                    'experimental_energy/status.json':'experimental_native/status.json',
                    'experimental_electronic_scientific_stage_started':'experimental_native_scientific_stage_started'}.items():
        if body.count(old)!=1:raise ValueError('Native single-pilot placement template changed')
        body=body.replace(old,new,1)
    scope={};exec(compile(body,str(Path(__file__).resolve()),'exec'),dict(vars(helper)),scope)
    helper.wait=scope['wait'];helper.main()


if __name__=='__main__':main()
