#!/usr/bin/env python
"""Apply the shared all-axis policy to the exact, source-reviewed TEMPO graft."""
import argparse
import json
from pathlib import Path

from redox.screening.feasibility import assess_candidate
from redox.workflows import polymer_campaign as api


def assess(route_path, screen_path, potential_path, policy_path):
    route=json.loads(route_path.read_text());screen=json.loads(screen_path.read_text())
    potential=json.loads(potential_path.read_text());policy=json.loads(policy_path.read_text())
    entries=screen['candidates']
    if (len(entries)!=1 or entries[0]['id']!='tempo' or entries[0]['route_status']!='proceed'
            or not entries[0]['proceed_to_full_redox'] or entries[0]['evidence']!=route['evidence']):
        raise ValueError('Exact route review and prospective gate differ')
    if (potential['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed'
            or potential['product_smiles']!=route['candidate']['product_smiles']):
        raise ValueError('Native potential belongs to a different local graft')
    api.verify_sources(potential['source_sha256'])
    key=potential['product_inchikey']
    rules={a['id']:a for a in policy['required_axes']}
    route_sources=[dict(path=str(p),sha256=api.digest(p)) for p in (route_path,screen_path)]
    axes={
        'attachment_route':dict(product_key=key,scope=rules['attachment_route']['scope'],status='qualified',
            value='proceed',sources=route_sources,
            context=dict(solvent='DMF or THF',base='NaH',support='Merrifield resin',source='US6441243B1 Example A4')),
        'redox_window':dict(product_key=key,scope=rules['redox_window']['scope'],status='partial',
            reason='Exact local model potential reconstructed; application window and predictive uncertainty are unresolved',
            model_potential_V_vs_Fc=potential['potential_V_vs_Fc_by_rotor_cutoff']['50'],
            sources=[dict(path=str(potential_path),sha256=api.digest(potential_path))])}
    candidate=dict(id='tempo',product_key=key,role='candidate',axes=axes,electrode_pools=['catholyte'])
    decision=assess_candidate(candidate,policy)
    if decision['status']!='pending' or decision['all_required_axes_pass']:
        raise ValueError('Missing material evidence must remain pending')
    return dict(status='exact_TEMPO_route_and_all_axis_decision_checked',candidate=candidate,decision=decision,
        source_sha256={str(p):api.digest(p) for p in (route_path,screen_path,potential_path,policy_path,Path(__file__).resolve())},
        primary_decision_ledger_modified=False,quantum_jobs_submitted=0,
        interpretation='Route proceed permits further evaluation only; it is not a synthesis or material shortlist.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('route','screen','potential','policy','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=assess(a.route.resolve(),a.screen.resolve(),a.potential.resolve(),a.policy.resolve())
    if a.out.exists():raise FileExistsError('Preserve an existing decision')
    a.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['decision']['status'],'pending_axes':len(result['decision']['pending_axes'])}))
