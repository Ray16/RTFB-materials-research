#!/usr/bin/env python
"""Reconstruct the retained first-case repeatability failure from saved calls."""
import argparse
import itertools
import json
import math
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def difference(left,right):
    if len(left['forces_eV_A'])!=len(right['forces_eV_A']):raise ValueError('Force shapes differ')
    energy=abs(left['energy_eV']-right['energy_eV'])
    force=max(abs(a-b) for x,y in zip(left['forces_eV_A'],right['forces_eV_A']) for a,b in zip(x,y))
    if not math.isfinite(energy) or not math.isfinite(force):raise ValueError('Nonfinite inference difference')
    return energy,force


def review(campaign,raw,out):
    if out.exists():raise FileExistsError('Preserve prior reviews')
    plan=api.read_plan(campaign);api.verify_sources(plan['input_sha256']);api.verify_sources(plan['producer_sha256'])
    case=plan['cases'][0];sources={str(p):api.digest(p) for p in (campaign/'plan.json',raw/'status.json',Path(__file__).resolve())}
    modes={}
    for mode in ('reference','nitro'):
        folder=raw/mode;points=json.loads((folder/'evaluation_diagnostics.json').read_text());stored=json.loads((folder/'repeat_diagnostic.json').read_text());state=json.loads((folder/'status.json').read_text())
        api.verify_sources(points['source_sha256']);api.verify_sources(stored['source_sha256'])
        if state['mode']!=mode or stored['case_id']!=case['id'] or state['plan_sha256']!=api.digest(campaign/'plan.json'):
            raise ValueError('Wrong inference identity')
        calls=points['calls']
        for i,row in enumerate(calls):
            if (row['call']!=i or row['case_id']!=case['id'] or not 0<=row['frame_index']<4
                    or len(row['forces_eV_A'])!=case['atoms']
                    or any(len(f)!=3 for f in row['forces_eV_A'])
                    or not all(math.isfinite(v) for f in row['forces_eV_A'] for v in f)
                    or not math.isfinite(row['energy_eV'])):
                raise ValueError('Malformed point-call inventory')
        same=[r for r in calls if r['frame_index']==0]
        if len(same)<3:raise ValueError('Multiple exact-frame calls required')
        energy,force=difference(same[-1],same[0]);etol=plan['repeat_energy_tolerance_eV_per_atom']*case['atoms'];ftol=plan['repeat_force_tolerance_eV_A']
        if (abs(energy-stored['repeat_energy_error_eV'])>1e-12 or abs(force-stored['repeat_force_error_eV_A'])>1e-12
                or stored['energy_tolerance_eV']!=etol or stored['force_tolerance_eV_A']!=ftol):
            raise ValueError('Repeat diagnostic does not reconstruct')
        later=[difference(a,b) for a,b in itertools.combinations(same[1:],2)]
        modes[mode]=dict(calls=len(calls),exact_frame_repeats=len(same),initial_to_final_energy_error_eV=energy,
                         initial_to_final_force_error_eV_A=force,energy_tolerance_eV=etol,force_tolerance_eV_A=ftol,
                         repeat_gate_passed=energy<=etol and force<=ftol,
                         maximum_later_pair_force_difference_eV_A=max(x[1] for x in later),
                         checkpoint_sha256=state['model_identity']['checkpoint_sha256'],GPU=state['GPU'],host=state['host'],
                         device_name=state['device_name'],source_paths=[str(folder/n) for n in ('evaluation_diagnostics.json','repeat_diagnostic.json','status.json')])
        sources.update({str(folder/name):api.digest(folder/name) for name in ('evaluation_diagnostics.json','repeat_diagnostic.json','status.json')})
    if (not modes['reference']['repeat_gate_passed'] or modes['nitro']['repeat_gate_passed']
            or any(modes['reference'][k]!=modes['nitro'][k] for k in ('checkpoint_sha256','GPU','host','device_name'))):
        raise ValueError('Expected same-device reference pass and Nitro failure were not demonstrated')
    result=dict(status='first_case_Nitro_repeatability_failure_independently_reconstructed',checked=time.time(),
                case_id=case['id'],atoms=case['atoms'],modes=modes,source_sha256=sources,
                production_calculator_qualified=False,full_panel_qualified=False,
                scope='Retained first-case stability failure under declared repeat thresholds. Does not establish the implementation cause or qualify altered tolerances.')
    api.verify_sources(sources);out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('campaign','raw','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();r=review(a.campaign.resolve(),a.raw.resolve(),a.out.resolve());print(json.dumps(dict(status=r['status'],modes=r['modes'])))
