#!/usr/bin/env python3
"""Review exact component controls and release only independently checked pilots."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.component_sensitivity import reconstruct


def collect(manifest_path,output):
    manifest=json.loads(manifest_path.read_text());output.mkdir(parents=True,exist_ok=True)
    definition=importlib.util.spec_from_file_location('component_energy_reviewer',ROOT/'scripts/validation/review_screening_energy.py')
    reviewer=importlib.util.module_from_spec(definition);definition.loader.exec_module(reviewer)
    rows=manifest['states'];checked=[];pending=[];failures=[]
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:
            raise ValueError('Manifest-bound namespace changed')
        root=Path(row['namespace']).parent;review_path=output/(row['key']+'.json')
        if review_path.exists():
            result=json.loads(review_path.read_text())
            if any(sha256_file(p)!=h for p,h in result['source_sha256'].items()):raise ValueError('Reviewed source changed')
        else:
            status=root/'case/status.json'
            if not status.exists():pending.append(dict(key=row['key'],status='not_started'));continue
            stage=json.loads(status.read_text())
            if 'finished' not in stage:pending.append(dict(key=row['key'],status='running_not_reviewed'));continue
            if stage['status']!='screening_energy_control_finished_review_pending':
                failures.append(dict(key=row['key'],status=stage['status'],error=stage.get('error')));continue
            placements=[p for p in root.glob('placement_*') if (p/'completion.json').exists()
                        and json.loads((p/'completion.json').read_text()).get('returncode')==0]
            if len(placements)!=1:pending.append(dict(key=row['key'],status='normal_placement_closure_pending'));continue
            try:result=reviewer.review(Path(row['namespace']),placements[0])
            except (OSError,ValueError,KeyError) as error:
                failures.append(dict(key=row['key'],status='independent_review_failed',error=repr(error)));continue
            atomic_write_json(review_path,result)
        if (result['key']!=row['key'] or result['status']!='benchmark_screening_energy_control_independently_checked'
                or Path(result['namespace']).resolve()!=Path(row['namespace']).resolve()
                or result['source_state']!=row['source_state']):
            raise ValueError('Wrong state or independent-review contract')
        checked.append(result)
    pilot_keys={max([r for r in rows if (r['source_state']['mult']>1)==opened],key=lambda r:(r['nao'],-abs(r['source_state']['charge'])))['key'] for opened in (False,True)}
    if pilot_keys<={r['key'] for r in checked}:
        receipt=manifest_path.parent/'resource_pilots.json'
        if not receipt.exists():atomic_write_json(receipt,reviewer.pilots([output/(key+'.json') for key in sorted(pilot_keys)],manifest_path))
    complete=len(checked)==manifest['target_states'] and not pending and not failures
    report=dict(status='exact_component_controls_independently_checked' if complete else 'component_controls_incomplete',
        workload=manifest.get('workload','exact_viologen_package'),target_states=manifest['target_states'],checked_states=len(checked),pending=pending,failures=failures,
        checked=time.time(),source_sha256={str(manifest_path):sha256_file(manifest_path),**{str(output/(r['key']+'.json')):sha256_file(output/(r['key']+'.json')) for r in checked}},
        new_free_energy_qualified=False, observer_owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
    if any(r['key'] in pilot_keys for r in failures):
        report['status']='resource_pilot_failed_remaining_states_held'
    elif failures and not pending:
        report['status']='component_controls_terminal_with_failures'
    if complete:
        package=json.loads((manifest_path.parent/'inputs/package_manifest.json').read_text());energies={}
        for row in checked:
            case=Path(row['namespace']).parent/'case';plan=json.loads((case/'plan.json').read_text())
            energies[row['key']]=dict(stored_Ha=plan['reference_energy_Ha'],
                installed_Ha=json.loads((case/'installed/result.json').read_text())['e_smd_Ha'],
                corrected_Ha=json.loads((case/'corrected/result.json').read_text())['e_smd_Ha'],
                charge=row['source_state']['charge'],mult=row['source_state']['mult'])
        report['sensitivity']=reconstruct(package,energies)
        report['all_declared_agreement_gates_passed']=all(r['declared_agreement_passed'] for r in checked)
        atomic_write_json(output/'scalar_inputs.json',dict(package=package,energies=energies))
    atomic_write_json(output/'summary.json',report);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--watch-hours',type=float,default=0);a=p.parse_args()
    if not 0<=a.watch_hours<=24:raise ValueError('Bounded 0--24 hour observation required')
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    deadline=time.time()+a.watch_hours*3600
    while True:
        result=collect(a.manifest.resolve(),a.out.resolve())
        print(json.dumps({k:result[k] for k in ('status','checked_states','target_states','failures')}),flush=True)
        if result['status'] in ('exact_component_controls_independently_checked','resource_pilot_failed_remaining_states_held','component_controls_terminal_with_failures') or time.time()>=deadline:break
        time.sleep(10)
