"""Prepare and run a bounded, read-only orbital-response performance diagnostic.

The production 64-GiB disk gate remains unchanged. Short probes reuse large
input arrays and cap output at 64 MiB. The explicitly selected checkpointed
profile reserves a 2-GiB output budget above a 32-GiB shared-volume margin.
Its recreated grid is not an exact native-Hessian restart.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.qm.response_diagnostic import read_checkpoint
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows import molecular_native_namespace as execution
from redox.workflows.archived_native_namespace import verify_spec,remaining_host_promises


def prepare(namespace,output,debugger=None,checkpoint_response=False):
    source=verify_spec(namespace);reader=NamespaceSources(source);source_plan=Path(source['attempt'])/'plan.json'
    parent=json.loads(source_plan.read_text());reader.verify_all(parent['producer_sha256']);reader.verify_all(parent['input_sha256'])
    if len(parent['tasks'])!=1 or (parent['tasks'][0]['charge'],parent['tasks'][0]['mult'])!=(0,1):raise ValueError('One declared neutral singlet method control required')
    if output.exists() or any(output.is_relative_to(Path(source[k]).resolve()) for k in ('archive','prior_root','attempt','raw_parent')):raise FileExistsError('Use a new durable diagnostic namespace outside every scientific source')
    raw=Path(source['raw_parent'])/source['key'];checkpoint,_,_=read_checkpoint(raw)
    if checkpoint['source_sha256'][str(raw/'geometry.xyz')]!=source['geometry_sha256']:raise ValueError('Checkpoint geometry differs from the declared scored input')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    paths=[(ROOT/'scripts/validation/run_response_diagnostic.py','worker.py'),(ROOT/'src/redox/qm/response_diagnostic.py','helper.py'),(ROOT/'src/redox/workflows/owned_debugger.py','debugger.py'),(ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(Path(__file__),'preparer.py')]
    for path,name in paths:shutil.copyfile(path,sources/name)
    bindings={**{str(reader.path(p)):h for p,h in {**parent['producer_sha256'],**parent['input_sha256']}.items()},**checkpoint['source_sha256']}
    for path in [namespace,source_plan,*sources.glob('*'),*[Path(source['archive'])/name for name in ('scripts/fleet/exclusive_python_worker.py','scripts/fleet/tensor_runtime_worker.py','environments/cutensor22/runtime.json')]]:bindings[str(path)]=sha256_file(path)
    from redox.workflows.owned_debugger import configuration
    debug=None if debugger is None else configuration(debugger)
    limits=(7200,2*2**30,34) if checkpoint_response else (900,64*2**20,2)
    plan=dict(role='orbital_response_diagnostic',source_id='orbital_response_diagnostic',state='rest',source_namespace=str(namespace),source_key=source['key'],raw_source=str(raw),checkpoint=checkpoint,debugger=debug,checkpoint_response=checkpoint_response,
        symbols=[line.split()[0] for line in (raw/'geometry.xyz').read_text().splitlines()[2:]],source_sha256=bindings,
        maximum_host_RSS_GiB=80,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,maximum_wall_seconds=limits[0],maximum_durable_output_bytes=limits[1],
        minimum_diagnostic_free_GiB=limits[2],native_G_qualified=False,material_admission=False,original_grid_arrays_available=False,
        available_response_memory_cap_GiB=3,original_native_memory_policy_context_required=True,
        scope='Saved orbitals and first-derivative frame; recreated initial grid recipe, fixed-density energy check and traced response only. No new SCF or native-G admission.')
    atomic_write_json(case/'plan.json',plan)
    spec=dict(prior_root=source['prior_root'],archive=source['archive'],case=str(case),sources=str(sources),source_id=plan['source_id'],state='rest',plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=80,minimum_free_GiB=32,scientific_sources_read_only=True)
    atomic_write_json(output/'namespace.json',spec);reader.recheck();return spec


def probe(path):
    spec=execution.verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    program=('import ast,json,importlib.metadata as m,importlib.util,pathlib;import numpy,h5py,psutil,ase;'
        'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
        'sources=pathlib.Path('+repr(spec['sources'])+');'
        'modules=sorted({n.module for path in [sources/name for name in ["worker.py","helper.py","debugger.py","host_guard.py"]] for n in ast.walk(ast.parse(path.read_text())) if isinstance(n,ast.ImportFrom) and n.module and n.module.startswith("redox.")});'
        'assert all(importlib.util.find_spec(name) is not None for name in modules),modules;'
        'definition=importlib.util.spec_from_file_location("response_helper",'+repr(spec['sources']+'/helper.py')+');'
        'helper=importlib.util.module_from_spec(definition);definition.loader.exec_module(helper);'
        'plan=json.loads(pathlib.Path('+repr(spec['case']+'/plan.json')+').read_text());'
        'bundle,_,_=helper.read_checkpoint(pathlib.Path(plan["raw_source"]));assert bundle==plan["checkpoint"];'
        'print(json.dumps(dict(status="saved_response_frame_and_CPU_dependencies_checked",versions=wanted)))')
    run=subprocess.run(execution.command(spec)+[runtime['python'],'-c',program],capture_output=True,text=True,check=True,timeout=180)
    result=dict(host=os.uname().nodename.split('.')[0],namespace_sha256=sha256_file(path),checks=json.loads(run.stdout),GPU_started=False)
    atomic_write_json(path.with_name('node_probe.json'),result);return result


def launch(path,gpu,output):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation,gpu_suitability,free_slots,GATE
    spec=execution.verify(path);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());host=os.uname().nodename.split('.')[0]
    limits=(plan['maximum_wall_seconds'],plan['maximum_durable_output_bytes'],plan['minimum_diagnostic_free_GiB'])
    expected=(7200,2*2**30,34) if plan.get('checkpoint_response') else (900,64*2**20,2)
    if (plan['role']!='orbital_response_diagnostic' or limits!=expected
            or plan['native_G_qualified'] or not spec['scientific_sources_read_only']):raise ValueError('Only the explicitly bounded diagnostic has this storage admission')
    proof=json.loads(path.with_name('node_probe.json').read_text())
    if proof['host']!=host or proof['namespace_sha256']!=sha256_file(path):raise ValueError('Exact node-local preflight required')
    if output.exists() or any((case/name).exists() for name in ('status.json','residency_guard','host_guard')):raise FileExistsError('Preserve all previous diagnostics')
    runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text());prior=Path(spec['prior_root'])
    output.mkdir(parents=True);shutil.copyfile(Path(__file__),output/'launcher_source.py')
    argv=execution.command(spec);locks=GATE.parent.parent/'.gpu_locks';argv[-1:-1]=['--bind',str(locks),str(locks)]
    # This is a prospective stricter guard for the response-only diagnostic.
    # A successful run must remain below 40 GiB and also satisfies its original
    # 80-GiB plan bound. No running reservation is changed or retrospectively
    # shrunk; the launch record retains both values and withholds qualification.
    effective_bound=40
    argv += [str(GATE),'run',str(gpu),'--label','redox-response-diagnostic','--',runtime['python'],str(prior/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(case/'residency_guard'),str(Path(spec['sources'])/'host_guard.py'),'--out',str(case/'host_guard'),'--maximum-rss-gib',str(effective_bound),'--minimum-free-gib','32',str(prior/'scripts/fleet/tensor_runtime_worker.py'),'--profile',str(prior/'environments/cutensor22/runtime.json'),'--out',str(case/'tensor_runtime'),str(Path(spec['sources'])/'worker.py'),'--plan',str(case/'plan.json')]
    env=dict(os.environ,PYTHONPATH=str(prior/'src'),LD_LIBRARY_PATH=runtime['library_path'],PYTHONDONTWRITEBYTECODE='1',CUPY_CACHE_DIR='/tmp/cupy_kernel_cache',REDOX_HESSIAN_MEMORY_GIB='6',TMPDIR='/tmp',NUMPY_MADVISE_HUGEPAGE='0')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
    with allocation_guard():
        capacity=host_capacity(2,minimum_ram_gib=effective_bound+32);promises=remaining_host_promises();free=shutil.disk_usage(case).free
        if free<plan['minimum_diagnostic_free_GiB']*2**30 or capacity['available_RAM_bytes']-promises<(effective_bound+32)*2**30:raise RuntimeError('Diagnostic storage or host-RAM admission unavailable')
        cpus=cpu_allocation(2);device=gpu_suitability(gpu)
        if (host,gpu) not in free_slots([host]):raise RuntimeError('GPU occupied; no diagnostic started')
        with (output/'command.log').open('x') as log:
            child=subprocess.Popen(['taskset','--cpu-list',','.join(map(str,cpus)),*argv],env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,cpus))
        result=dict(status='bounded_response_diagnostic_submitted_start_unverified',host=host,GPU=gpu,pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=cpus,device=device,capacity=capacity,remaining_RAM_promises_bytes=promises,free_storage_bytes=free,namespace=str(path),namespace_sha256=sha256_file(path),argv=argv,started=time.time(),native_G_qualified=False,production_64_GiB_gate_changed=False,maximum_wall_seconds=900,maximum_durable_output_bytes=64*2**20,numpy_madvise_hugepage=0,shared_kernel_policy_changed=False,launcher_source_sha256=sha256_file(output/'launcher_source.py'))
        result.update(original_plan_maximum_host_RSS_GiB=80,effective_maximum_host_RSS_GiB=effective_bound,prospective_resource_profile_qualified=False,
            maximum_wall_seconds=plan['maximum_wall_seconds'],maximum_durable_output_bytes=plan['maximum_durable_output_bytes'],checkpoint_response=plan.get('checkpoint_response',False),
            budget_observer_owner=dict(host=host,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
        atomic_write_json(output/'launch.json',result)
    return child,result


def observe_budget(child,spec,result,output):
    """Signal only the exact owned scientific process if this probe exceeds budget."""
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    atomic_write_json(output/'observer.json',dict(owner=result['budget_observer_owner'],CPU_affinity=sorted(os.sched_getaffinity(0)),status='bounded_response_budget_observer_running'))
    case=Path(spec['case']);reason=None;requested=None;term_sent=False;kill_sent=False;debugger=None
    plan=json.loads((case/'plan.json').read_text());debug_module=None
    if plan.get('debugger') is not None:
        debug_module=execution.load(Path(spec['sources'])/'debugger.py','frozen_owned_response_debugger')
    while child.poll() is None:
        if debug_module is not None and debugger is None and (case/'debugger_wait.json').exists():
            debugger=debug_module.start(case,plan['debugger'],str(case/'plan.json'))
        allocated=0
        for path in case.parent.rglob('*'):
            try:
                if path.is_file():allocated+=path.stat().st_blocks*512
            except FileNotFoundError:
                # An owned atomic checkpoint rename can retire a temporary
                # filename between discovery and stat; the next poll sees its
                # published name. Do not lose the budget observer in that race.
                continue
        if reason is None and (time.time()-result['started']>plan['maximum_wall_seconds'] or allocated>plan['maximum_durable_output_bytes']):
            reason='wall_time_budget' if time.time()-result['started']>plan['maximum_wall_seconds'] else 'durable_output_budget';requested=time.monotonic()
            atomic_write_json(output/'budget_stop.json',dict(reason=reason,native_G_qualified=False,checked=time.time()))
        if reason:
            guard=case/'residency_guard/owner.json'
            if guard.exists():
                owner=json.loads(guard.read_text());pid=owner['pid'];proc=Path('/proc')/str(pid)
                if (proc.exists() and proc.stat().st_uid==os.getuid()
                        and proc.joinpath('stat').read_text().rsplit(')',1)[1].split()[19]==owner['process_start']
                        and str(case/'plan.json').encode() in proc.joinpath('cmdline').read_bytes()):
                    if not term_sent:os.kill(pid,signal.SIGTERM);term_sent=True
                    elif not kill_sent and time.monotonic()-requested>=15:os.kill(pid,signal.SIGKILL);kill_sent=True
            elif not any((case/name).exists() for name in ('status.json','tensor_runtime/status.json')):
                proc=Path('/proc')/str(child.pid)
                if (proc.exists() and proc.joinpath('stat').read_text().rsplit(')',1)[1].split()[19]==result['process_start']
                        and os.getpgid(child.pid)==child.pid):
                    if not term_sent:os.killpg(child.pid,signal.SIGTERM);term_sent=True
                    elif not kill_sent and time.monotonic()-requested>=15:os.killpg(child.pid,signal.SIGKILL);kill_sent=True
        time.sleep(2)
    if debugger is not None:
        try:debugger.wait(timeout=5)
        except subprocess.TimeoutExpired:
            debugger.terminate()
            try:debugger.wait(timeout=3)
            except subprocess.TimeoutExpired:debugger.kill();debugger.wait(timeout=3)
        result['debugger_returncode']=debugger.returncode
        stop=case/'debugger_stop.json'
        if stop.exists():result['debugger_stop']=json.loads(stop.read_text())
    result.update(status='bounded_response_diagnostic_returned_analysis_pending',returncode=child.returncode,finished=time.time(),budget_stop=reason)
    atomic_write_json(output/'completion.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','probe','launch'));p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path);p.add_argument('--gpu',type=int);p.add_argument('--debugger',type=Path);p.add_argument('--checkpoint-response',action='store_true');a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.spec.resolve(),a.out.resolve(),None if a.debugger is None else a.debugger.resolve(),a.checkpoint_response)
    elif a.action=='probe':result=probe(a.spec.resolve())
    else:
        spec=execution.verify(a.spec.resolve());child,result=launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True);result=observe_budget(child,spec,result,a.out.resolve())
    print(json.dumps(result))
