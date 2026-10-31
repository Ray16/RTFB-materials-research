"""Source-frozen small GPU controls for the upstream screening-barrier fix."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows import molecular_native_namespace as execution
from redox.workflows.archived_native_namespace import remaining_host_promises
from redox.workflows.response_diagnostic import observe_budget


def prepare(archive,output,water_reference=False,molecular_control=None):
    if molecular_control not in (None,'closed','open') or (water_reference and molecular_control):raise ValueError('Choose one declared numerical control')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):raise FileExistsError('Use a fresh durable control outside scientific checkouts')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    worker='run_screening_molecular_control.py' if molecular_control else 'run_screening_control.py'
    for path,name in [(ROOT/'scripts/validation'/worker,'worker.py'),(ROOT/'src/redox/qm/screening_barrier.py','helper.py'),(ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),(ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(Path(__file__),'preparer.py')]:shutil.copyfile(path,sources/name)
    if water_reference:shutil.copyfile(ROOT/'tests/fixtures/screening_water39.json',sources/'water.json')
    if molecular_control:shutil.copyfile(ROOT/'tests/fixtures/screening_molecular_controls.json',sources/'molecule.json')
    package=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
    files=[*sources.glob('*'),*[p for p in (archive/'src/redox/core').glob('*.py')],archive/'src/redox/workflows/resources.py',archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',archive/'environments/cutensor22/runtime.json',package/'dft/gen_grid.py',package/'dft/numint.py',package/'lib/libgdft.so']
    bound,seconds,storage=(4,1800,64*2**20) if water_reference else (2,600,32*2**20)
    if molecular_control:
        bound,seconds,storage=8,7200,512*2**20
        files.extend(archive/'src/redox/qm'/name for name in ('dft.py','cds_finite.py'))
    plan=dict(role='screening_molecular_control' if molecular_control else 'screening_barrier_control',source_id='screening_barrier_control',state='rest',source_sha256={str(p):sha256_file(p) for p in files},maximum_host_RSS_GiB=bound,minimum_free_GiB=32,maximum_wall_seconds=seconds,maximum_durable_output_bytes=storage,water_reference=water_reference,molecular_control=molecular_control,installed_library_files_modified=False,implementation_qualified=False,native_G_qualified=False)
    atomic_write_json(case/'plan.json',plan);spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=plan['source_id'],state='rest',plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=bound,minimum_free_GiB=32,scientific_sources_read_only=True)
    atomic_write_json(output/'namespace.json',spec);return spec


def probe(path):
    spec=execution.verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    program=('import ast,importlib.util,importlib.metadata as m,json,pathlib;import numpy,psutil,ase,h5py;from pyscf.gto import mole;'
        'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
        'root=pathlib.Path('+repr(spec['sources'])+');'
        'modules=sorted({n.module for name in ["worker.py","helper.py","host_guard.py","scientific_worker.py","checkpoints.py","response_helper.py","restart_helper.py"] if (root/name).exists() for n in ast.walk(ast.parse((root/name).read_text())) if isinstance(n,ast.ImportFrom) and n.module and n.module.startswith("redox.")});'
        'assert all(importlib.util.find_spec(name) is not None for name in modules),modules;'
        'print(json.dumps(dict(versions=wanted,worker_imports=modules,CUDA_context_started=False)))')
    run=subprocess.run(execution.command(spec)+[runtime['python'],'-c',program],capture_output=True,text=True,check=True,timeout=120)
    result=dict(host=os.uname().nodename.split('.')[0],namespace_sha256=sha256_file(path),checks=json.loads(run.stdout),GPU_started=False);atomic_write_json(path.with_name('node_probe.json'),result);return result


def prepare_polymer(namespace,output,variant):
    """Keep the existing three-point scientific diagnostic and all source ties."""
    from redox.workflows.polymer_gradient_diagnostic import prepare as prepare_gradient
    if variant not in ('installed','corrected'):raise ValueError('Explicit screening variant required')
    spec=prepare_gradient(namespace,output)
    case=Path(spec['case']);sources=Path(spec['sources']);plan_path=case/'plan.json'
    plan=json.loads(plan_path.read_text())
    shutil.copyfile(plan_path,case/'gradient_preparation.json')
    shutil.copyfile(sources/'worker.py',sources/'scientific_worker.py')
    for original,name in [(ROOT/'scripts/validation/run_screening_polymer_control.py','worker.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),
            (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),(Path(__file__),'screening_preparer.py')]:
        shutil.copyfile(original,sources/name)
    package=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
    files=[*sources.glob('*'),case/'gradient_preparation.json',package/'dft/gen_grid.py',package/'lib/libgdft.so']
    plan['source_sha256'].update({str(p):sha256_file(p) for p in files})
    plan.update(screening_variant=variant,maximum_wall_seconds=3600,maximum_durable_output_bytes=512*2**20,
                implementation_qualified=False,production_default_changed=False)
    atomic_write_json(plan_path,plan);spec['plan_sha256']=sha256_file(plan_path)
    atomic_write_json(output/'namespace.json',spec)
    return spec


def prepare_polymer_hessian(namespace,reference,output,resume_from=None):
    """A new corrected derivative control; never a mixed-implementation G."""
    from redox.validation.namespace_sources import read_checked_review,NamespaceSources
    from redox.workflows.archived_native_namespace import verify_spec
    checked,binding=read_checked_review(reference,'native')
    ref=verify_spec(reference);original=verify_spec(namespace)
    if (ref['key']!=original['key'] or ref['geometry_sha256']!=original['geometry_sha256']
            or ref['method_id']!=original['method_id']):raise ValueError('Exact matched native reference required')
    reader=NamespaceSources(ref);raw=reader.path(ref['original_raw_root'])
    spec=prepare_polymer(namespace,output,'corrected')
    case=Path(spec['case']);sources=Path(spec['sources']);plan_path=case/'plan.json'
    plan=json.loads(plan_path.read_text());shutil.copyfile(plan_path,case/'gradient_control_preparation.json')
    for path,name in [(ROOT/'scripts/validation/run_screening_polymer_hessian.py','worker.py'),
            (ROOT/'src/redox/qm/response_diagnostic.py','response_helper.py'),
            (ROOT/'src/redox/qm/hessian_restart.py','restart_helper.py')]:shutil.copyfile(path,sources/name)
    inputs=[reference,Path(binding['review']),raw/'hessian.npy',raw/'geometry.xyz',raw/'experimental_partition.json']
    files=[*sources.glob('*'),case/'gradient_control_preparation.json',*inputs]
    archive=Path(spec['archive'])
    for name in ('qrrho.py','vibrations.py','cds_finite.py','hessian_memory.py'):
        files.append(archive/'src/redox/qm'/name)
    plan['source_sha256'].update({str(p):sha256_file(p) for p in files})
    native=execution.load(archive/'scripts/pipeline/smooth_polymer_native.py','full_Hessian_control_original_settings')
    from pyscf import gto
    from redox.qm.hessian_restart import response_ceiling_gib
    center=json.loads(Path(plan['points'][0]['evaluation']).read_text())
    mol=gto.M(atom=list(zip(plan['symbols'],center['coordinates_Bohr'])),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0)
    response_cap,minimum_response=response_ceiling_gib(mol.nao,mol.nelectron//2)
    plan.update(role='screening_polymer_Hessian_control',full_polymer_Hessian_control=True,
        numerical_settings=native.NATIVE_SETTINGS,response_available_memory_GiB=response_cap,
        response_nao=mol.nao,response_nocc=mol.nelectron//2,minimum_aligned_response_bytes=minimum_response,maximum_wall_seconds=14400,
        maximum_durable_output_bytes=4*2**30,minimum_initial_storage_GiB=36,minimum_live_storage_GiB=32,
        checked_native_reference=dict(namespace=str(reference),review=binding,
            Hessian=str(raw/'hessian.npy'),partition=str(raw/'experimental_partition.json'),geometry=str(raw/'geometry.xyz')),
        native_G_qualified=False,chemical_reaction_G_assembled=False)
    if resume_from is not None:
        from redox.qm.hessian_restart import read_frame
        resume=execution.verify(resume_from);resume_case=Path(resume['case'])
        old=json.loads((resume_case/'plan.json').read_text());stage=json.loads((resume_case/'status.json').read_text())
        if (stage['status']!='screening_polymer_Hessian_control_failed_preserved'
                or stage.get('error')!="RuntimeError('GPU memory insufficient for solving CPHF equations')"
                or old['original_namespace']!=plan['original_namespace'] or old['screening_variant']!='corrected'
                or not json.loads((resume_case/'residency_guard/owner.json').read_text()).get('guard_stopped')):
            raise ValueError('Only a closed matched response-workspace refusal may supply these components')
        owners=execution.load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','restart_closed_owner')
        if owners.alive(stage['owner']) is not False:raise ValueError('Original Hessian owner must be proven gone')
        frame,_=read_frame(resume_case)
        plan['resume_components']=dict(namespace=str(resume_from),source_case=str(resume_case),frame=frame)
        plan['source_sha256'].update(frame['source_sha256'])
        for p in [resume_from,resume_case/'plan.json',resume_case/'status.json',resume_case/'screening_implementation/implementation.json']:
            plan['source_sha256'][str(p)]=sha256_file(p)
    atomic_write_json(plan_path,plan);spec['plan_sha256']=sha256_file(plan_path)
    atomic_write_json(output/'namespace.json',spec);reader.recheck()
    return spec


def launch(path,gpu,output):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation,gpu_suitability,free_slots,GATE
    spec=execution.verify(path);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text());host=os.uname().nodename.split('.')[0];proof=json.loads(path.with_name('node_probe.json').read_text())
    expected=(4,1800,64*2**20) if plan.get('water_reference') else (2,600,32*2**20)
    if plan.get('molecular_control'):expected=(8,7200,512*2**20)
    role='screening_molecular_control' if plan.get('molecular_control') else 'screening_barrier_control'
    if plan.get('screening_variant'):
        expected=(80,3600,512*2**20);role='polymer_gradient_diagnostic'
    if plan.get('full_polymer_Hessian_control'):
        expected=(80,14400,4*2**30);role='screening_polymer_Hessian_control'
        if (plan.get('minimum_initial_storage_GiB'),plan.get('minimum_live_storage_GiB'))!=(36,32):raise ValueError('Full-Hessian control storage margins changed')
    if plan.get('energy_control'):
        from redox.workflows.screening_energy_panel import output_budget
        budget=output_budget(plan['nao'],plan['source_state']['mult'],plan['has_iron'])
        expected=(40,7200,budget);role='screening_energy_control'
        if (plan.get('minimum_initial_storage_GiB'),plan.get('minimum_live_storage_GiB'))!=(32+budget/2**30,32):raise ValueError('Energy-control storage margins changed')
    if plan.get('charged_polymer_pilot') and not plan.get('charged_polymer_optimization'):
        expected=(40,7200,512*2**20);role='charged_polymer_SCF_gradient_resource_pilot'
        if (plan.get('minimum_initial_storage_GiB'),plan.get('minimum_live_storage_GiB'))!=(32.5,32):raise ValueError('Charged-pilot storage margins changed')
    if plan.get('charged_polymer_optimization'):
        expected=(40,28800,256*2**20);role='charged_polymer_solution_geometry_pilot'
        if (plan.get('minimum_initial_storage_GiB'),plan.get('minimum_live_storage_GiB'))!=(32.25,32):raise ValueError('Charged-geometry storage margins changed')
    if plan['role']!=role or (plan['maximum_host_RSS_GiB'],plan['maximum_wall_seconds'],plan['maximum_durable_output_bytes'])!=expected:raise ValueError('Exact bounded control resource class required')
    if proof['host']!=host or proof['namespace_sha256']!=sha256_file(path):raise ValueError('Exact node-local preflight required')
    if output.exists() or any((case/name).exists() for name in ('status.json','residency_guard','host_guard')):raise FileExistsError('Preserve earlier controls')
    output.mkdir(parents=True);shutil.copyfile(__file__,output/'launcher_source.py');shutil.copyfile(observe_budget.__code__.co_filename,output/'budget_observer_source.py');runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text());prior=Path(spec['prior_root'])
    argv=execution.command(spec);locks=GATE.parent.parent/'.gpu_locks';argv[-1:-1]=['--bind',str(locks),str(locks)]
    argv += [str(GATE),'run',str(gpu),'--label','redox-screening-control','--',runtime['python'],str(prior/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(case/'residency_guard'),str(Path(spec['sources'])/'host_guard.py'),'--out',str(case/'host_guard'),'--maximum-rss-gib',str(plan['maximum_host_RSS_GiB']),'--minimum-free-gib','32',str(prior/'scripts/fleet/tensor_runtime_worker.py'),'--profile',str(prior/'environments/cutensor22/runtime.json'),'--out',str(case/'tensor_runtime'),str(Path(spec['sources'])/'worker.py'),'--plan',str(case/'plan.json')]
    env=dict(os.environ,PYTHONPATH=str(prior/'src'),LD_LIBRARY_PATH=runtime['library_path'],PYTHONDONTWRITEBYTECODE='1',CUPY_CACHE_DIR='/tmp/cupy_kernel_cache',NUMPY_MADVISE_HUGEPAGE='0',TMPDIR='/tmp')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
    with allocation_guard():
        capacity=host_capacity(2,minimum_ram_gib=plan['maximum_host_RSS_GiB']+32);promises=remaining_host_promises()
        storage_floor=plan.get('minimum_initial_storage_GiB',2)
        if capacity['available_RAM_bytes']-promises<(plan['maximum_host_RSS_GiB']+32)*2**30 or shutil.disk_usage(case).free<storage_floor*2**30:raise RuntimeError('Control RAM or output headroom unavailable')
        cpus=cpu_allocation(2);device=gpu_suitability(gpu)
        if plan.get('charged_polymer_pilot') and device['memory_MiB']<plan.get('minimum_GPU_memory_MiB',31*1024):
            raise RuntimeError('The charged-polymer reference needs its declared larger-memory device')
        if (host,gpu) not in free_slots([host]):raise RuntimeError('GPU occupied; no screening control started')
        with (output/'command.log').open('x') as log:child=subprocess.Popen(['taskset','--cpu-list',','.join(map(str,cpus)),*argv],env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,cpus))
        result=dict(status='screening_control_submitted_start_unverified',host=host,GPU=gpu,pid=child.pid,process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=cpus,capacity=capacity,remaining_RAM_promises_bytes=promises,device=device,argv=argv,namespace=str(path),namespace_sha256=sha256_file(path),started=time.time(),maximum_wall_seconds=plan['maximum_wall_seconds'],maximum_durable_output_bytes=plan['maximum_durable_output_bytes'],budget_observer_owner=dict(host=host,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),budget_observer_source_sha256=sha256_file(output/'budget_observer_source.py'),native_G_qualified=False,production_storage_gate_changed=False)
        atomic_write_json(output/'launch.json',result)
    return child,result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','prepare-polymer','prepare-Hessian','probe','launch'));p.add_argument('--spec',type=Path);p.add_argument('--reference',type=Path);p.add_argument('--resume-components',type=Path);p.add_argument('--archive',type=Path);p.add_argument('--out',type=Path);p.add_argument('--gpu',type=int);p.add_argument('--water-reference',action='store_true');p.add_argument('--molecular-control',choices=('closed','open'));p.add_argument('--screening-variant',choices=('installed','corrected'));a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.archive.resolve(),a.out.resolve(),a.water_reference,a.molecular_control)
    elif a.action=='prepare-polymer':result=prepare_polymer(a.spec.resolve(),a.out.resolve(),a.screening_variant)
    elif a.action=='prepare-Hessian':result=prepare_polymer_hessian(a.spec.resolve(),a.reference.resolve(),a.out.resolve(),None if a.resume_components is None else a.resume_components.resolve())
    elif a.action=='probe':result=probe(a.spec.resolve())
    else:
        spec=execution.verify(a.spec.resolve());child,result=launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True);result=observe_budget(child,spec,result,a.out.resolve())
    print(json.dumps(result))
