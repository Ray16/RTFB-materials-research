#!/usr/bin/env python3
"""Prepare and dispatch the ten resting charged occupancy/configuration graphs.

One bounded dispatch pass uses every suitable free slot needed by the roster.
Reruns preserve uncertain/failed attempts; they never reclaim an owner by age.
Each scientific worker has its own gate, CPU/RAM checks and output budget.
"""
import argparse
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def prepare(preparation,agreement,existing,archive,output,storage_peers):
    from redox.workflows.charged_polymer_pilot import prepare_optimization
    if output.exists():raise FileExistsError('Preserve the existing panel')
    old=verify(existing);old_plan=json.loads((Path(old['case'])/'plan.json').read_text())
    if not old_plan.get('charged_polymer_optimization'):raise ValueError('Existing largest geometry pilot required')
    data=json.loads(preparation.read_text());selected=[r for r in data['cases'] if r['case']['occupied_sites']]
    if len(selected)!=10:raise ValueError('Two configurations and five occupied graphs each required')
    output.mkdir(parents=True);rows=[]
    for row in selected:
        key=row['case']['key']
        if key==old_plan['prepared_case']['key']:namespace=existing
        else:
            folder=output/'states'/key
            prepare_optimization(archive,preparation,agreement,folder,key)
            namespace=folder/'namespace.json'
        spec=verify(namespace);plan=json.loads((Path(spec['case'])/'plan.json').read_text())
        if plan['prepared_case']!=row['case']:raise ValueError('Exact graph/state differs from the roster')
        rows.append(dict(key=key,namespace=str(namespace),namespace_sha256=sha256_file(namespace),
                         charge=plan['charge'],mult=plan['mult'],atoms=plan['atoms']))
    peers=[]
    for path in storage_peers:
        verify(path);peers.append(dict(namespace=str(path),namespace_sha256=sha256_file(path)))
    manifest=dict(status='ten_exact_resting_ionic_geometry_models_prepared',target_states=10,states=rows,
        storage_peers=peers,preparation=str(preparation),preparation_sha256=sha256_file(preparation),
        agreement=str(agreement),agreement_sha256=sha256_file(agreement),
        conformer_sampling_qualified=False,native_G_qualified=False,material_admission=False,
        producer_sha256=sha256_file(__file__))
    atomic_write_json(output/'manifest.json',manifest);return manifest


def inspect(row,alive):
    root=Path(row['namespace']).parent;case=root/'case';stage=case/'status.json'
    if stage.exists():
        data=json.loads(stage.read_text())
        if 'finished' in data:
            return ('finished_review_pending' if data['status']=='charged_polymer_solution_geometry_finished_independent_review_pending'
                    else 'failed_preserved'),data
        live=alive(data['owner'])
        return ('running' if live is True else 'stopped_preserved' if live is False else 'ownership_unverified'),data
    guard=case/'residency_guard/owner.json'
    if guard.exists():
        data=json.loads(guard.read_text());live=alive(data)
        return ('guard_started' if live is True else 'stopped_preserved' if live is False else 'ownership_unverified'),data
    claim=root/'dispatch_claim'
    if not claim.exists():return 'pending',{}
    # A submitted/uncertain attempt is retained even if its controller vanished;
    # a separate evidence-based prelaunch recovery is required to reopen it.
    receipts=sorted(claim.glob('try_*.json'))
    if not receipts:return 'uncertain_submission_preserved',{}
    receipt=json.loads(receipts[-1].read_text());live=alive(receipt['owner'])
    return ('submitted' if live is True else 'ownership_unverified' if live is None else 'closed_submission_preserved'),receipt


def storage_headroom(rows,peers,additional,alive,free_bytes):
    """Reserve each running/uncertain worker's remaining bounded output once."""
    remaining=0;seen=set()
    for row in [*rows,*peers]:
        namespace=Path(row['namespace'])
        if namespace in seen:continue
        seen.add(namespace);spec=verify(namespace);case=Path(spec['case']);root=case.parent
        status=case/'status.json';active=False
        if status.exists():
            data=json.loads(status.read_text());active=('finished' not in data and alive(data['owner']) is not False)
        elif (case/'residency_guard/owner.json').exists():
            active=alive(json.loads((case/'residency_guard/owner.json').read_text())) is not False
        elif (root/'dispatch_claim').exists():active=True
        if active:
            budget=json.loads((case/'plan.json').read_text())['maximum_durable_output_bytes']
            written=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
            remaining+=max(0,budget-written)
    return dict(free_bytes=free_bytes,remaining_existing_budget_bytes=remaining,
        additional_budget_bytes=additional,minimum_margin_bytes=32*2**30,
        sufficient=free_bytes>=32*2**30+remaining+additional)


def dispatch(manifest_path,interpreter):
    utilities=load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py','ionic_panel_capacity')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','ionic_panel_owners')
    remote=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','ionic_panel_remote_launcher')
    data=json.loads(manifest_path.read_text());rows=data['states'];peers=data['storage_peers']
    if len(rows)!=10 or len({r['key'] for r in rows})!=10:raise ValueError('Complete distinct occupancy roster required')
    for row in [*rows,*peers]:
        if sha256_file(row['namespace'])!=row['namespace_sha256']:raise ValueError('Prepared namespace changed')
        verify(Path(row['namespace']))
    policy=ROOT.parent/'.config/gpu_reserve/config.env';match=re.search(r'^ALLOWED_HOSTS="([^"]+)"',policy.read_text(),re.M)
    if match is None:raise ValueError('System gate policy required')
    output=manifest_path.parent;record=dict(started=time.time(),manifest_sha256=sha256_file(manifest_path),
        source_sha256={str(p):sha256_file(p) for p in [Path(__file__),policy,Path(remote.__file__),Path(utilities.__file__)]},
        launches=[],refusals=[],native_G_qualified=False)
    with (output/'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        slots,errors=utilities.discover_slots(match.group(1).split());record['availability_errors']=errors
        pending=[r for r in rows if inspect(r,owners.alive)[0]=='pending'];blocked=set()
        for host,gpu in slots:
            if not pending:break
            if host in blocked:continue
            row=pending[0];spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
            storage=storage_headroom(rows,peers,plan['maximum_durable_output_bytes'],owners.alive,shutil.disk_usage(output).free)
            if not storage['sufficient']:
                record['refusals'].append(dict(reason='aggregate_durable_storage_budget',**storage));break
            try:capacity=utilities.capacity(host,gpu,interpreter)
            except (OSError,subprocess.SubprocessError) as error:
                record['refusals'].append(dict(host=host,GPU=gpu,error=repr(error)));continue
            if not utilities.fits(plan,capacity):continue
            claim=Path(row['namespace']).parent/'dispatch_claim';claim.mkdir()
            atomic_write_json(claim/'owner.json',dict(queue=str(output),host=os.uname().nodename,pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
            pending.pop(0)
            try:
                receipt=remote.launch(row,host,gpu,interpreter,claim,1)
                record['launches'].append(dict(key=row['key'],host=host,GPU=gpu,receipt=receipt,storage=storage))
            except (OSError,subprocess.SubprocessError) as error:
                atomic_write_json(claim/'uncertain_001.json',dict(host=host,GPU=gpu,error=repr(error)))
                record['refusals'].append(dict(key=row['key'],host=host,error=repr(error)));blocked.add(host)
            atomic_write_json(output/'dispatch_current.json',record)
        record.update(finished=time.time(),status='bounded_dispatch_pass_finished_workers_preserved',
            states=[dict(key=r['key'],status=inspect(r,owners.alive)[0]) for r in rows])
        atomic_write_json(output/('dispatch_'+str(time.time_ns())+'.json'),record)
        atomic_write_json(output/'dispatch_current.json',record)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','dispatch'))
    for name in ('preparation','agreement','existing','archive','out','manifest'):p.add_argument('--'+name,type=Path)
    p.add_argument('--storage-peer',type=Path,action='append',default=[]);p.add_argument('--python')
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=(prepare(a.preparation.resolve(),a.agreement.resolve(),a.existing.resolve(),a.archive.resolve(),a.out.resolve(),[x.resolve() for x in a.storage_peer])
        if a.action=='prepare' else dispatch(a.manifest.resolve(),a.python))
    print(json.dumps({k:v for k,v in result.items() if k in ('status','target_states','launches','refusals','states')}))
