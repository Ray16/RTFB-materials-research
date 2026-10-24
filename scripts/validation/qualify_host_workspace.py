#!/usr/bin/env python
"""Profile the checked host-backed adapter on a smaller physical device class.

Reuse the exact resource-pilot implementation and require a matched large-device
reference. This qualifies a starting-geometry footprint, never a minimum or G.
"""
import json
import math
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


ADAPTER=ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py'
CHECKER=ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py'


def comparison(plan,directory):
    previous=Path(plan['comparison_reference_profile'])/'energy/pilot/energy_and_gradient.json'
    current=directory/'energy_and_gradient.json'
    a,b=(json.loads(p.read_text()) for p in (previous,current))
    differences=dict(maximum_energy_difference_eV=max(abs(a['result'][k]-b['result'][k])
        for k in ('e_smd_eV','e_gas_eV')),
        gradient_max_difference_Ha_Bohr=abs(a['gradient_max_Ha_Bohr']-b['gradient_max_Ha_Bohr']))
    if (any(not math.isfinite(v) for v in differences.values())
            or differences['maximum_energy_difference_eV']>1e-5
            or differences['gradient_max_difference_Ha_Bohr']>1e-6):
        raise ValueError('Smaller-device point differs from its exact large-device reference')
    return dict(compared=True,differences=differences,
        source_sha256={str(p):api.digest(p) for p in (previous,current)},
        scope='Identical starting geometry/model; not optimized-state or G qualification')


def main():
    original=ROOT/'scripts/validation/qualify_polymer_workspace.py'
    module=api.load_module('host_workspace_original_resource_pilot',original)
    module.ADAPTER=ADAPTER
    def checked(paths,budget):
        roots={Path(p).parents[2] for p in paths}
        if len(paths)!=2 or len(roots)!=1:
            raise ValueError('One matched full-control spin pair required')
        return api.load_module('host_workspace_full_controls',CHECKER).checked_controls(roots.pop(),budget)
    module.checked_controls=checked
    reference=None;maximum=None
    if sys.argv[1]=='prepare':
        for flag in ('--reference-profile','--maximum-peak-GPU-MiB'):
            if flag not in sys.argv:raise ValueError('Declare reference and physical peak ceiling before preparation')
        index=sys.argv.index('--reference-profile');reference=Path(sys.argv[index+1]).resolve();del sys.argv[index:index+2]
        index=sys.argv.index('--maximum-peak-GPU-MiB');maximum=int(sys.argv[index+1]);del sys.argv[index:index+2]
        minimum=int(sys.argv[sys.argv.index('--minimum-GPU-MiB')+1])
        if maximum<=0 or minimum-maximum<1024:
            raise ValueError('Keep at least one GiB physical device headroom')
    loader=api.load_module
    def load(name,path):
        loaded=loader(name,path)
        if name=='workspace_continuity_check':loaded.compare_previous_pilot=comparison
        return loaded
    api.load_module=load
    try:module.main()
    finally:api.load_module=loader
    if sys.argv[1]=='prepare':
        campaign=Path(sys.argv[sys.argv.index('--campaign')+1]).resolve();plan=api.read_plan(campaign)
        old=api.read_plan(reference);api.verify_sources(old['producer_sha256']);api.verify_sources(old['input_sha256'])
        if (plan['tasks']!=old['tasks'] or any(plan[k]!=old[k] for k in ('quantum_settings','thermochemistry','solvent'))):
            raise ValueError('Reference belongs to a different exact task or model')
        key=plan['tasks'][0]['key'];a=campaign/'live_inputs'/key/'input.xyz';b=reference/'live_inputs'/key/'input.xyz'
        if api.digest(a)!=api.digest(b):raise ValueError('Reference geometry differs')
        status=reference/'energy/pilot/status.json';record=json.loads(status.read_text())
        if record['status']!='memory_pilot_passed':raise ValueError('Complete large-device point required')
        api.verify_sources(record['output_sha256'])
        plan.update(comparison_reference_profile=str(reference),maximum_memory_pilot_MiB=maximum,
            resource_policy=dict(workspace_estimate_GiB=plan['execution_workspace_GiB'],
                measured_peak_ceiling_MiB=maximum,minimum_device_MiB=minimum,
                minimum_device_headroom_MiB=minimum-maximum,scientific_criteria_changed=False))
        for p in (reference/'plan.json',status,status.parent/'energy_and_gradient.json',b):
            plan['input_sha256'][str(p)]=api.digest(p)
        for p in (Path(__file__).resolve(),original,ADAPTER,CHECKER):plan['producer_sha256'][str(p)]=api.digest(p)
        atomic_write_json(campaign/'plan.json',plan)


if __name__=='__main__':main()
