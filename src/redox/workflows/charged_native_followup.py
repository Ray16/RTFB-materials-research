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
import hashlib
import shlex
import subprocess
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify

SOURCE_BYTES=Path(__file__).read_bytes()


def priority_requests(path):
    """Read declared, immutable job identities; this grants no scientific admission."""
    path=Path(path)
    if not path.exists():return []
    data=json.loads(path.read_text())
    if data.get('schema')!=1 or not isinstance(data.get('requests'),list):
        raise ValueError('Declared dispatch-priority schema required')
    rows=data['requests']
    if len({r['namespace'] for r in rows})!=len(rows):raise ValueError('Duplicate priority job')
    for row in rows:
        if sha256_file(row['namespace'])!=row['namespace_sha256'] or row.get('deferred_slots')!=1:
            raise ValueError('Priority binds one exact namespace and one free slot')
    return rows


def priority_for_slot(requests,candidate,observed,already_deferred,alive,fits):
    """Leave one suitable free slot for a live pending lead; running jobs are untouched."""
    for row in requests:
        namespace=Path(row['namespace'])
        if str(namespace)==candidate['namespace'] or str(namespace) in already_deferred:continue
        spec=json.loads(namespace.read_text());case=Path(spec['case'])
        # A started, failed or finished scientific attempt cannot request an
        # extra GPU or be restarted by this advisory scheduling rule.
        if (case/'status.json').exists() or (case/'residency_guard/owner.json').exists():continue
        status=Path(row['controller_status'])
        if not status.exists():continue
        controller=json.loads(status.read_text())
        if ('finished' in controller or controller.get('namespace')!=str(namespace)
                or alive(controller['owner']) is not True):continue
        plan_path=case/'plan.json'
        if sha256_file(plan_path)!=spec['plan_sha256']:raise ValueError('Priority plan changed')
        plan=json.loads(plan_path.read_text())
        if (observed['GPU_MiB']<plan.get('minimum_GPU_memory_MiB',31*1024)
                or observed.get('GPU_compute_capability',0)<plan.get('minimum_GPU_compute_capability',0)
                or not fits(plan,observed)):continue
        return row,plan
    return None


def device_capability(host,gpu):
    argv=['nvidia-smi','--id',str(gpu),'--query-gpu=compute_cap','--format=csv,noheader,nounits']
    result=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(argv)],
                          capture_output=True,text=True,timeout=25,check=True)
    return float(result.stdout.strip())


def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def prerequisite(manifest,review_path):
    from redox.workflows.native_component_resources import LARGE_CLASS,EXPANDED_RHF_CLASS,peak_limit,validate_plan
    large=manifest.get('native_resource_class') in (LARGE_CLASS,EXPANDED_RHF_CLASS)
    if manifest.get('native_resource_class') and not large:raise ValueError('Unknown native followup resource class')
    review=json.loads(review_path.read_text())
    if (review['status']!='charged_polymer_native_component_independently_reconstructed' or not review['native_G_qualified']
            or Path(review['namespace']).resolve()!=Path(manifest['prerequisite_namespace']).resolve()):
        raise ValueError('Exact completed first charged-native pilot required')
    if any(sha256_file(p)!=h for p,h in review['source_sha256'].items()):raise ValueError('Pilot backing changed')
    original=verify(Path(review['namespace']));pilot=json.loads((Path(original['case'])/'plan.json').read_text())
    raw=json.loads((Path(original['case'])/'thermodynamic_component.json').read_text())
    if large and review.get('native_resource_class')!=manifest['native_resource_class']:raise ValueError('Exact completed larger dimension pilot required')
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
    return dispatch_rows(rows,manifest_path.parent,interpreter,[manifest_path,review_path],avoid_hosts)


def dispatch_rows(rows,output,interpreter,prerequisite_sources,avoid_hosts=(),storage_peers=(),eligible_keys=None):
    """Bounded shared launcher; callers first qualify their scientific class."""
    from redox.workflows.resources import GATE,free_slots
    capacity=load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py','charged_native_capacity')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','charged_native_owners')
    remote=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','charged_native_shared_launcher')
    storage=load(ROOT/'scripts/fleet/run_charged_polymer_geometry_panel.py','charged_native_shared_storage')
    policy=Path(os.environ.get('GPU_RESERVE_CONFIG',str(GATE.parent.parent/'.config/gpu_reserve/config.env'))).resolve()
    match=re.search(r'^ALLOWED_HOSTS="([^"]+)"',policy.read_text(),re.M)
    if match is None:raise ValueError('Authoritative system host policy required')
    priority_path=ROOT/'config/campaign_dispatch_priorities.json'
    priorities=priority_requests(priority_path)
    snapshots={}
    for source in [*prerequisite_sources,Path(__file__),policy,
                   *([priority_path] if priority_path.exists() else []),
                   *[Path(m.__file__) for m in (capacity,owners,remote,storage)]]:
        source=Path(source);content=SOURCE_BYTES if source==Path(__file__) else source.read_bytes()
        digest=hashlib.sha256(content).hexdigest();frozen=output/'source_snapshots'/digest/source.name
        if not frozen.exists():frozen.parent.mkdir(parents=True,exist_ok=True);frozen.write_bytes(content)
        if sha256_file(frozen)!=digest:raise ValueError('As-used dispatch source snapshot differs')
        snapshots[str(source)]=dict(snapshot=str(frozen),sha256=digest)
    record=dict(status='checked_native_followup_dispatching',started=time.time(),source_snapshots=snapshots,
        source_sha256={row['snapshot']:row['sha256'] for row in snapshots.values()},launches=[],refusals=[],native_G_qualified=False,avoided_hosts=sorted(set(avoid_hosts)))
    record['eligible_keys']=None if eligible_keys is None else sorted(eligible_keys)
    record['priority_deferrals']=[]
    with (output/'followup_dispatch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        pending=[]
        for row in rows:
            if eligible_keys is not None and row['key'] not in eligible_keys:continue
            root=Path(row['namespace']).parent
            condition,_=remote.inspect(row,root/'dispatch_claim',owners.alive)
            if condition in ('pending','prelaunch_refusal_no_scientific_attempt'):pending.append(row)
        hosts=[host for host in match.group(1).split() if host not in avoid_hosts]
        slots,errors=capacity.discover_slots(hosts);record['availability_errors']=errors;blocked_hosts=set()
        # The fast inventory can time out on a healthy many-GPU host. One
        # authoritative 60-second gate query distinguishes that observation
        # limit from an unavailable node; capacity and atomic reservation still
        # run for every returned device. No worker or claim is restarted here.
        record['availability_retries']=[]
        for error in errors:
            host=error['host']
            if host not in hosts:continue
            try:
                recovered=free_slots([host])
                if any(name!=host for name,_ in recovered):raise ValueError('Gate returned an unrequested host')
                slots.extend(slot for slot in recovered if slot not in slots)
                record['availability_retries'].append(dict(host=host,status='authoritative_gate_query_completed',slots=recovered))
            except Exception as retry_error:
                record['availability_retries'].append(dict(host=host,status='unverified_no_submission',error=repr(retry_error)))
        deferred=set()
        for host,gpu in slots:
            if not pending:break
            if host in blocked_hosts:continue
            row=pending[0]
            claim=None
            try:
                observed=capacity.capacity(host,gpu,interpreter)
                # A returned SSH receipt can precede the Python guard's start.
                # Account for its RAM until the live process scan can see that
                # worker, rather than assigning the same headroom twice.
                starting=sum(previous['maximum_host_RSS_GiB']*2**30 for previous in record['launches']
                    if previous['host']==host and not (Path(previous['receipt']['placement']).parent/'case/residency_guard/owner.json').exists())
                observed=dict(observed,promised=observed['promised']+starting)
                selected=None
                # A large input at the head of the queue must not hide a
                # smaller ready input that fits this node's remaining RAM.
                for index,candidate in enumerate(pending):
                    hint=json.loads(Path(candidate['namespace']).read_text())
                    candidate_plan=json.loads((Path(hint['case'])/'plan.json').read_text())
                    if observed['GPU_MiB']>=candidate_plan.get('minimum_GPU_memory_MiB',31*1024) and capacity.fits(candidate_plan,observed):
                        selected=(index,candidate);break
                if selected is None:continue
                index,row=selected
                if priorities and any(r['namespace']!=row['namespace'] for r in priorities):
                    observed['GPU_compute_capability']=device_capability(host,gpu)
                    priority=priority_for_slot(priorities,row,observed,deferred,owners.alive,capacity.fits)
                    if priority is not None:
                        request,priority_plan=priority
                        room=storage.storage_headroom(rows,storage_peers,priority_plan['maximum_durable_output_bytes'],
                            owners.alive,shutil.disk_usage(output).free)
                        if room['sufficient']:
                            deferred.add(request['namespace'])
                            record['priority_deferrals'].append(dict(host=host,GPU=gpu,
                                namespace=request['namespace'],controller_status=request['controller_status'],
                                reason=request['reason'],resources=observed,storage=room,
                                existing_workers_interrupted=False))
                            atomic_write_json(output/'followup_dispatch_current.json',record)
                            continue
                spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text())
                if observed['GPU_MiB']<plan.get('minimum_GPU_memory_MiB',31*1024) or not capacity.fits(plan,observed):
                    raise ValueError('Selected resource hints changed before verification')
                space=storage.storage_headroom(rows,storage_peers,plan['maximum_durable_output_bytes'],owners.alive,shutil.disk_usage(output).free)
                if not space['sufficient']:record['refusals'].append(dict(reason='aggregate_storage',**space));break
                candidate_claim=Path(row['namespace']).parent/'dispatch_claim'
                if candidate_claim.exists():
                    if json.loads((candidate_claim/'owner.json').read_text()).get('queue')!=str(output):
                        raise ValueError('A foreign dispatch claim remains preserved')
                    condition,_=remote.inspect(row,candidate_claim,owners.alive)
                    if condition!='prelaunch_refusal_no_scientific_attempt':raise ValueError('Only a proven closed prelaunch refusal may retry')
                else:
                    candidate_claim.mkdir()
                    atomic_write_json(candidate_claim/'owner.json',dict(queue=str(output),host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
                claim=candidate_claim
                pending.pop(index)
                number=len(list(claim.glob('try_*.json')))+1
                receipt=remote.launch(row,host,gpu,interpreter,claim,number)
                record['launches'].append(dict(key=row['key'],host=host,GPU=gpu,receipt=receipt,maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB']))
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
