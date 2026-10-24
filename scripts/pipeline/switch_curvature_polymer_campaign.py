#!/usr/bin/env python
"""Native stages with tested host switch curvature and explicit device headroom."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER=ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py'


def checked_controls(folder,budget=28):
    base=api.load_module('switch_underlay_control_checks',ROOT/'scripts/validation/qualify_pcm_first_response_v2.py')
    proof=base.checked_controls(folder,budget)
    for mult in (1,2):
        control=folder/'controls'/str(mult)
        stats_path=control/'compact/execution/switch_curvature.json'
        stats=json.loads(stats_path.read_text())
        if (stats['factories']<1 or stats['primitive_calls']<1 or stats['host_tensor_bytes']<=0
                or stats['solute_checkpoints']<1 or stats['live_binding_refreshes']<1
                or stats['equations_changed'] or stats['precision_changed'] or stats['surface_grid_changed']):
            raise ValueError('Host switch curvature, checkpoints and live bindings must all be exercised')
        api.verify_sources(stats['source_sha256'])
        guard=folder/'residency_guards'/str(mult)
        owner_path,state_path=guard/'owner.json',guard/'state.json'
        owner,state=(json.loads(p.read_text()) for p in (owner_path,state_path))
        reference=json.loads((control/'status.json').read_text())
        if (owner['pid']!=reference['pid'] or owner['process_start']!=reference['process_start']
                or owner.get('guard_stopped') is not True or state['foreign_resident_PIDs']
                or not state['own_residency_observed'] or (guard/'interruption.json').exists()
                or (guard/'observation_error.json').exists()):
            raise ValueError('Full controls lack an uninterrupted residency audit')
        for p in (stats_path,owner_path,state_path):proof[str(p)]=api.digest(p)
        proof.update(stats['source_sha256'])
    return proof


def device_limit(maximum,minimum):
    if not isinstance(maximum,int) or not isinstance(minimum,int) or maximum<=0 or minimum-maximum<1024:
        raise ValueError('Declared minimum device must retain at least one GiB above the measured-peak ceiling')
    return dict(maximum_native_peak_MiB=maximum,minimum_device_MiB=minimum,
        minimum_device_headroom_MiB=minimum-maximum,workspace_estimate_GiB=28,
        resource_and_scientific_criteria_separate=True)


def main():
    own=Path(__file__).resolve()
    first_path=ROOT/'scripts/pipeline/first_response_polymer_campaign.py'
    first=api.load_module('preserved_response_campaign_stages',first_path)
    first.ADAPTER=ADAPTER
    first.CHECKER=own
    first.__file__=str(own)
    ceiling=None
    if sys.argv[1]=='prepare-native':
        flag='--maximum-native-GPU-MiB'
        if flag not in sys.argv or '--minimum-GPU-MiB' not in sys.argv:
            raise ValueError('Declare both the physical device class and native measured-peak ceiling')
        index=sys.argv.index(flag);maximum=int(sys.argv[index+1]);del sys.argv[index:index+2]
        minimum=int(sys.argv[sys.argv.index('--minimum-GPU-MiB')+1])
        ceiling=device_limit(maximum,minimum)
    first.main()
    if sys.argv[1]=='prepare-native':
        out=Path(sys.argv[sys.argv.index('--out')+1]).resolve()
        plan=api.read_plan(out)
        plan['maximum_memory_pilot_MiB']=ceiling['maximum_native_peak_MiB']
        plan['native_device_limit_policy']=ceiling
        plan['purpose']='Fresh full native resource/minimum/partition qualification; exact host switch curvature, original physical gates'
        for p in (own,first_path,ROOT/'scripts/fleet/exclusive_python_worker.py'):
            plan['producer_sha256'][str(p)]=api.digest(p)
        atomic_write_json(out/'plan.json',plan)


if __name__=='__main__':
    main()
