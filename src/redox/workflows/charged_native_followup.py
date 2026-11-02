"""Release prepared small charged-native components after their exact pilot passes.

Uses the existing gate, node probe/launcher, capacity utility and aggregate
storage accounting. One bounded dispatch pass considers every suitable free GPU;
uncertain, occupied or failed attempts remain preserved for explicit review.
"""
import argparse
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import shutil
import time
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def prerequisite(manifest,review_path):
    from redox.workflows.native_component_resources import LARGE_CLASS,peak_limit,validate_plan
    large=manifest.get('native_resource_class')==LARGE_CLASS
    if manifest.get('native_resource_class') and not large:raise ValueError('Unknown native followup resource class')
    review=json.loads(review_path.read_text())
    if (review['status']!='charged_polymer_native_component_independently_reconstructed' or not review['native_G_qualified']
            or Path(review['namespace']).resolve()!=Path(manifest['prerequisite_namespace']).resolve()):
        raise ValueError('Exact completed first charged-native pilot required')
    if any(sha256_file(p)!=h for p,h in review['source_sha256'].items()):raise ValueError('Pilot backing changed')
    original=verify(Path(review['namespace']));pilot=json.loads((Path(original['case'])/'plan.json').read_text())
    raw=json.loads((Path(original['case'])/'thermodynamic_component.json').read_text())
    if large and review.get('native_resource_class')!=LARGE_CLASS:raise ValueError('Exact completed larger dimension pilot required')
    if large:validate_plan(pilot)
    if (review['charge']!=(4 if large else 2) or review['mult']!=1 or set(review['G_by_cutoff_eV'])!={'25','50','100'}
            or any(type(v) not in (int,float) or not math.isfinite(v) for v in review['G_by_cutoff_eV'].values())
            or not 0<review['peak_GPU_MiB']<=peak_limit(pilot) or not 0<=review['maximum_column_error']<=.01
            or not math.isfinite(review['minimum_frequency_cm']) or review['minimum_frequency_cm'] < -5
            or any(abs(v-raw['G_by_cutoff_eV'][k])>1e-8 for k,v in review['G_by_cutoff_eV'].items())):
        raise ValueError('Complete matching physical/resource pilot quantities required')
    rows=[]
    for row in manifest['prepared']:
        spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
        if (not plan.get('charged_polymer_native_component') or plan['response_nao']>pilot['response_nao']
                or plan['response_nocc']>pilot['response_nocc'] or plan['mult']!=1 or plan['charge'] not in ((4,) if large else (0,2))
                or plan['solvent_model']!=pilot['solvent_model'] or plan['numerical_settings']!=pilot['numerical_settings']):
            raise ValueError('Followup exceeds the checked method/dimension class')
        if large:validate_plan(plan)
        names=('worker.py','helper.py','kernel.cu','checkpoints.py','response_helper.py','restart_helper.py','response_resume_helper.py')
        if large:names+=('native_resource_helper.py',)
        for name in names:
            if sha256_file(Path(spec['sources'])/name)!=sha256_file(Path(original['sources'])/name):
                raise ValueError('Followup scientific implementation differs from the checked pilot')
        rows.append(dict(key=row['key'],namespace=row['namespace'],namespace_sha256=sha256_file(row['namespace'])))
    if len({r['key'] for r in rows})!=len(rows):raise ValueError('Distinct prepared components required')
    return rows


def dispatch(manifest_path,review_path,interpreter,avoid_hosts=()):
    manifest=json.loads(manifest_path.read_text());rows=prerequisite(manifest,review_path)
    capacity=load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py','charged_native_capacity')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','charged_native_owners')
    remote=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','charged_native_shared_launcher')
    storage=load(ROOT/'scripts/fleet/run_charged_polymer_geometry_panel.py','charged_native_shared_storage')
    policy=Path(os.environ.get('GPU_RESERVE_CONFIG',str(Path.home()/'.config/gpu_reserve/config.env'))).resolve()
    match=re.search(r'^ALLOWED_HOSTS="([^"]+)"',policy.read_text(),re.M)
    if match is None:raise ValueError('Authoritative system host policy required')
    output=manifest_path.parent;record=dict(status='checked_native_followup_dispatching',started=time.time(),
        source_sha256={str(p):sha256_file(p) for p in [manifest_path,review_path,Path(__file__),policy]},launches=[],refusals=[],native_G_qualified=False,avoided_hosts=sorted(set(avoid_hosts)))
    with (output/'followup_dispatch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        pending=[]
        for row in rows:
            root=Path(row['namespace']).parent
            if not (root/'case/status.json').exists() and not (root/'case/residency_guard/owner.json').exists() and not (root/'dispatch_claim').exists():pending.append(row)
        hosts=[host for host in match.group(1).split() if host not in avoid_hosts]
        slots,errors=capacity.discover_slots(hosts);record['availability_errors']=errors;blocked_hosts=set()
        for host,gpu in slots:
            if not pending:break
            if host in blocked_hosts:continue
            row=pending[0];spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
            space=storage.storage_headroom(rows,[],plan['maximum_durable_output_bytes'],owners.alive,shutil.disk_usage(output).free)
            if not space['sufficient']:record['refusals'].append(dict(reason='aggregate_storage',**space));break
            claim=None
            try:
                observed=capacity.capacity(host,gpu,interpreter)
                if observed['GPU_MiB']<plan.get('minimum_GPU_memory_MiB',31*1024) or not capacity.fits(plan,observed):continue
                candidate_claim=Path(row['namespace']).parent/'dispatch_claim';candidate_claim.mkdir();claim=candidate_claim
                atomic_write_json(claim/'owner.json',dict(queue=str(output),host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
                pending.pop(0)
                receipt=remote.launch(row,host,gpu,interpreter,claim,1)
                record['launches'].append(dict(key=row['key'],host=host,GPU=gpu,receipt=receipt))
            except Exception as error:
                record['refusals'].append(dict(key=row['key'],host=host,GPU=gpu,error=repr(error),uncertain_attempt_preserved=claim is not None))
                blocked_hosts.add(host)
                if claim is not None:
                    atomic_write_json(claim/'uncertain_001.json',dict(host=host,GPU=gpu,error=repr(error),checked=time.time()))
            atomic_write_json(output/'followup_dispatch_current.json',record)
        record.update(status='checked_native_followup_bounded_pass_finished',finished=time.time(),unsubmitted=[r['key'] for r in pending])
        atomic_write_json(output/('followup_dispatch_'+str(time.time_ns())+'.json'),record);atomic_write_json(output/'followup_dispatch_current.json',record)
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--pilot-review',type=Path,required=True);p.add_argument('--python',required=True);p.add_argument('--avoid-host',action='append',default=[],help='Observed unhealthy host for this pass; system policy is unchanged');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=dispatch(a.manifest.resolve(),a.pilot_review.resolve(),a.python,a.avoid_host);print(json.dumps(dict(status=r['status'],launches=len(r['launches']),refusals=r['refusals'],unsubmitted=r['unsubmitted'])))
