#!/usr/bin/env python3
"""Review a pre-science GPU-residency interruption without changing its files."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def classify(coordinator,launch,guard,state,interruption,runtime,log,raw_exists,stage_exists):
    if (coordinator.get('status')!='launched_scientific_ownership_requires_review'
            or coordinator.get('scientific_starts')!=0 or not coordinator.get('finished')
            or launch.get('status')!='gate_started_scientific_start_pending'
            or guard.get('status')!='program_interrupted_or_failed_preserved'
            or runtime.get('status')!='runtime_or_task_failed_preserved' or not guard.get('finished')
            or not runtime.get('finished') or raw_exists or stage_exists):
        raise ValueError('Attempt is not a terminal pre-science interruption')
    if (not interruption.get('foreign_resident_PIDs') or interruption.get('foreign_processes_signaled') is not False
            or interruption.get('signal_target')!='own reserved process only'
            or interruption.get('pid')!=guard.get('pid') or runtime.get('pid')!=guard.get('pid')
            or state.get('own_residency_observed') is not False
            or state.get('foreign_resident_PIDs')!=interruption.get('foreign_resident_PIDs')
            or 'import gpu4pyscf.lib.cutensor as engine' not in log):
        raise ValueError('Foreign-residency interruption is not demonstrated before numerical work')
    if (not isinstance(guard.get('peak_RSS_bytes'),int) or not isinstance(guard.get('free_RAM_bytes'),int)
            or guard['peak_RSS_bytes']>=guard['maximum_RSS_GiB']*2**30
            or guard['free_RAM_bytes']<=guard['minimum_free_GiB']*2**30):
        raise ValueError('Sampled host-memory evidence is ambiguous')
    return dict(status='foreign_residency_pre_science_stop_recovery_reviewed',
        foreign_processes_signaled=False,own_residency_observed=False,scientific_starts=0,
        sampled_host_RSS_bytes=guard['peak_RSS_bytes'],sampled_available_RAM_bytes=guard['free_RAM_bytes'],
        explanation='Residency guard recorded a foreign process; the inner SIGTERM handler emitted a misleading host-RAM message. No sampled RAM bound breach was recorded.',
        scope='Only a proven pre-science interruption is recoverable; unknown or partially started attempts remain held.')


def local_owner_dead(owner):
    if owner['host']!=os.uname().nodename.split('.')[0]:raise ValueError('Coordinator is not on this host')
    path=Path('/proc')/str(int(owner['pid']))/'stat'
    try:fields=path.read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:return True
    return fields[19]!=str(owner['process_start']) or fields[0] in ('Z','X')


def remote_owners_dead(host,owners):
    if not re.fullmatch(r'lambda(?:0|1|2|4|5|6|11|12|13)',host):raise ValueError('Unrecognized host')
    values=[]
    for owner in owners:
        pid=int(owner['pid']);start=str(owner['process_start'])
        if pid<=0 or not start.isdigit():raise ValueError('Invalid process identity')
        values.extend([str(pid),start])
    code="""import pathlib,sys\nfor i in range(1,len(sys.argv),2):\n p=pathlib.Path('/proc')/sys.argv[i]/'stat'\n try: f=p.read_text().rsplit(')',1)[1].split()\n except FileNotFoundError: continue\n if f[19]==sys.argv[i+1] and f[0] not in ('Z','X'): sys.exit(4)\n"""
    command=['python3','-c',code,*values]
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(command)],
                          capture_output=True,text=True,timeout=20)
    if result.returncode!=0:raise ValueError('Remote process ownership is live or unverifiable')
    return True


def review(case,out):
    case=case.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve previous reviews')
    files={name:case/name for name in ('coordinator.json','launch.json','host_guard/status.json',
        'residency_guard/owner.json','residency_guard/state.json','residency_guard/interruption.json','tensor_runtime/status.json','worker.log','plan.json')}
    data={name:json.loads(p.read_text()) for name,p in files.items() if p.suffix=='.json'}
    plan=data['plan.json'];api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    result=classify(data['coordinator.json'],data['launch.json'],data['host_guard/status.json'],
        data['residency_guard/state.json'],data['residency_guard/interruption.json'],data['tensor_runtime/status.json'],
        files['worker.log'].read_text(),Path(plan['raw_root']).exists(),(case/'counterpoise').exists())
    coordinator=data['coordinator.json']['owner'];launch=data['launch.json']['owner'];worker=data['residency_guard/owner.json']
    if (worker['host']!=launch['host'] or worker['pid']!=data['host_guard/status.json']['pid']
            or worker['GPU']!=data['residency_guard/state.json']['GPU'] or not worker.get('finished')):
        raise ValueError('Guard and gate ownership differ')
    if not local_owner_dead(coordinator):raise ValueError('Coordinator still live')
    remote_owners_dead(launch['host'],[launch,worker])
    sources={str(p):api.digest(p) for p in [*files.values(),Path(__file__).resolve()]}
    result.update(checked=time.time(),case=str(case),original_gate_owner=launch,original_worker_owner=worker,
                  original_coordinator_owner=coordinator,all_original_owners_stopped=True,source_sha256=sources,
                  primary_candidate_admission=False,solvent_counterpoise_computed=False)
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.case,a.out);print(json.dumps(dict(status=result['status'],scientific_starts=result['scientific_starts'])))
