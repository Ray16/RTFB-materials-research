#!/usr/bin/env python3
"""Collect all declared electronic controls and reconstruct only their shifts."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

from redox.core.common import ROOT,HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def collect(manifest_path,output):
    manifest=json.loads(manifest_path.read_text());output.mkdir(parents=True,exist_ok=True)
    definition=importlib.util.spec_from_file_location('independent_screening_energy_review',ROOT/'scripts/validation/review_screening_energy.py')
    checker=importlib.util.module_from_spec(definition);definition.loader.exec_module(checker)
    checked=[];pending=[];errors=[];used={str(manifest_path):sha256_file(manifest_path)}
    for row in manifest['states']:
        namespace=Path(row['namespace']);root=namespace.parent;case=root/'case';path=case/'status.json'
        review_path=output/(row['key']+'.json')
        if review_path.exists():
            result=json.loads(review_path.read_text())
            if any(sha256_file(p)!=h for p,h in result['source_sha256'].items()):raise ValueError('Collected backing changed')
        else:
            if not path.exists():pending.append(dict(key=row['key'],status='not_started'));continue
            stage=json.loads(path.read_text())
            if 'finished' not in stage:pending.append(dict(key=row['key'],status='execution_pending'));continue
            if stage['status']!='screening_energy_control_finished_review_pending':
                errors.append(dict(key=row['key'],status=stage['status'],error=stage.get('error')));continue
            guard=case/'residency_guard/owner.json'
            if not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
                pending.append(dict(key=row['key'],status='guard_closure_pending'));continue
            placements=[p for p in root.glob('placement_*') if (p/'completion.json').exists()
                        and json.loads((p/'completion.json').read_text()).get('returncode')==0]
            if len(placements)!=1:errors.append(dict(key=row['key'],status='ambiguous_normal_placement'));continue
            try:result=checker.review(namespace,placements[0])
            except (ValueError,OSError,KeyError) as error:
                errors.append(dict(key=row['key'],status='independent_review_failed',error=repr(error)));continue
            atomic_write_json(review_path,result)
        if result['key']!=row['key'] or result['status']!='benchmark_screening_energy_control_independently_checked':
            raise ValueError('Wrong collected state or scope')
        used[str(review_path)]=sha256_file(review_path);checked.append(result)
    complete=len(checked)==manifest['target_states'] and not errors and not pending
    report=dict(status='all_screening_energy_states_independently_checked' if complete else 'screening_energy_panel_collection_incomplete',
        target_states=manifest['target_states'],checked_states=len(checked),pending=pending,errors=errors,
        source_sha256=used,native_G_qualified=False,new_redox_potentials_qualified=False)
    if complete:
        aliases={a:r for r in checked for a in r['source_state']['aliases']}
        benchmark=json.loads((manifest_path.parent/'inputs/benchmark.json').read_text());endpoints=[]
        for target in benchmark['outliers']:
            left,right=target['event'].split('->')
            a=aliases['benchmark:'+target['id']+':'+left];b=aliases['benchmark:'+target['id']+':'+right]
            if a['source_state']['charge']!=b['source_state']['charge']+1:
                raise ValueError('One-electron reduction orientation differs')
            fo=aliases['benchmark:ferrocene:ox'];fr=aliases['benchmark:ferrocene:neu']
            shifts={k.replace('_Ha','_V'):HARTREE_EV*(a['differences'][k]-b['differences'][k]-fo['differences'][k]+fr['differences'][k])
                    for k in ('installed_minus_stored_Ha','corrected_minus_stored_Ha','corrected_minus_installed_Ha')}
            endpoints.append(dict(id=target['id'],event=target['event'],electronic_contribution_to_potential_shift_V=shifts))
        report.update(target_endpoints=benchmark['admitted_events'],endpoints=endpoints,
            maximum_state_energy_differences_Ha={k:max(abs(r['differences'][k]) for r in checked) for k in checked[0]['differences']},
            maximum_electronic_potential_shifts_V={k:max(abs(r['electronic_contribution_to_potential_shift_V'][k]) for r in endpoints) for k in endpoints[0]['electronic_contribution_to_potential_shift_V']},
            all_declared_agreement_gates_passed=all(r['declared_agreement_passed'] for r in checked),
            interpretation='Electronic implementation sensitivity only. No new potential combines these energies with retained thermal terms.')
    report['producer_sha256']=sha256_file(__file__)
    atomic_write_json(output/'summary.json',report)
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=collect(args.manifest.resolve(),args.out.resolve())
    print(json.dumps({k:v for k,v in result.items() if k in ('status','checked_states','target_states','errors','maximum_state_energy_differences_Ha','maximum_electronic_potential_shifts_V')}))
