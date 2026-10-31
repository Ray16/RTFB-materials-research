"""Complete exact route-admitted molecular states in isolated archive namespaces.

Scientific programs and their controls remain source-frozen. Only private DFT,
lock and bounded-repair directories are writable at the original recorded paths.
"""
import argparse
import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP,load_record,record_path
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.solution_pair_namespace import validate_layout
from redox.workflows.archived_native_namespace import remaining_host_promises

LAUNCHER_SOURCE_BYTES=Path(__file__).read_bytes()

def load(path,name):
    definition=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(definition);definition.loader.exec_module(module);return module


def safe_component(value):
    if not value or Path(value).name!=value or value in ('.','..'):
        raise ValueError('One safe nonempty molecular ID/state component required')


def source_state(archive,gid,state,product):
    """Validate the immutable source before copying any candidate output."""
    directory=archive/'calcs/dft'/gid/state
    raw=json.loads((directory/'result.json').read_text())
    with (archive/'library/manifest.csv').open() as stream:
        rows=[r for r in csv.DictReader(stream) if (r['id'],r['state'])==(gid,state)]
    if len(rows)!=1:raise ValueError('One exact source-manifest state required')
    row=rows[0];canonical=lambda s:Chem.MolToSmiles(Chem.MolFromSmiles(s))
    if (type(raw['charge']) is not int or type(raw['mult']) is not int
            or canonical(row['smiles'])!=canonical(product) or int(row['charge'])!=raw['charge']
            or int(row['mult_hint'])!=raw['mult'] or raw['id']!=gid or raw['state']!=state
            or raw.get('optimized') is not True
            or any(raw.get(k)!=v for k,v in dict(opt_xc='r2scan',opt_basis='def2-svp',opt_disp='d4',solvent='acetonitrile').items())
            or hydrogen_topology_status(product,directory/'opt.xyz')!='ok'):
        raise ValueError('Source graph/state/full connectivity or geometry protocol differs')
    energy=load_record(directory,raw['charge'],raw['mult'])
    spin=None if energy is None else energy.get('smd_spin_contam')
    if (energy is None or energy.get('protocol')!=ACTIVE_SP
            or energy.get('software')!={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'}
            or type(spin) not in (int,float) or not math.isfinite(spin) or abs(spin)>.1):
        raise ValueError('Checked matched diffuse-basis electronic parent required')
    if any((directory/name).exists() for name in ('native_thermal','projected_thermal')):
        raise FileExistsError('Existing native attempts require explicit reuse/review, never replacement')
    return row,raw,record_path(directory,raw['charge'],raw['mult'])


def prepare(archive,registry,routes,assessment,candidate_alias,source_id,states,output):
    archive=Path(archive).resolve();output=Path(output).resolve()
    safe_component(source_id)
    for state in states:safe_component(state)
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Use a fresh durable destination outside both source checkouts')
    records=json.loads(Path(registry).read_text())['candidates']
    matches=[r for r in records if candidate_alias in r['aliases'] and source_id in r['aliases']]
    if len(matches)!=1:raise ValueError('One exact graph must bind both the enumerated alias and source ID')
    candidate=matches[0]
    declarations=[d for d in candidate['declarations'] if d['id']==candidate_alias]
    if len(declarations)!=1:raise ValueError('One exact declared candidate required')
    gate=load(ROOT/'scripts/pipeline/screen_grafting_routes.py','molecular_native_prospective_route')
    route=gate.screen(dict(id=candidate_alias,frag=declarations[0]['fragment'],product_smiles=candidate['product_smiles']),
        json.loads(Path(routes).read_text())['candidates'],json.loads(Path(assessment).read_text()))
    if not route['proceed_to_full_redox'] or not route.get('exact_route_binding_verified'):
        raise ValueError('The exact prospective attachment route must pass before native calculations')
    declared={s[0]:(s[1],s[2]) for s in declarations[0]['states']}
    if len(states)!=len(set(states)) or not set(states)<=set(declared):raise ValueError('Distinct declared redox states required')
    prepared=[]
    for state in states:
        row,raw,energy=source_state(archive,source_id,state,candidate['product_smiles'])
        if (raw['charge'],raw['mult'])!=declared[state]:raise ValueError('Enumerated and source charge/spin differ')
        prepared.append((state,row,raw,energy))
    output.mkdir(parents=True)
    for state,row,raw,energy in prepared:
        folder=output/state;case=folder/'case';sources=folder/'sources';case.mkdir(parents=True);sources.mkdir()
        for name in ('dft','locks','repairs'): (case/name).mkdir()
        target=case/'dft'/source_id/state;target.mkdir(parents=True)
        origin=archive/'calcs/dft'/source_id/state
        for name in ('opt.xyz','result.json'):shutil.copyfile(origin/name,target/name)
        (target/'sp').mkdir();shutil.copyfile(energy,target/'sp'/energy.name)
        copied={}
        for path,name in [(ROOT/'scripts/pipeline/run_molecular_native.py','worker.py'),
                          (ROOT/'scripts/validation/review_molecular_native.py','reviewer.py'),
                          (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
                          (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
                          (Path(__file__),'launcher_source.py')]:
            dest=sources/name;shutil.copyfile(path,dest);copied[str(dest)]=sha256_file(dest)
        atomic_write_json(sources/'candidate.json',candidate)
        atomic_write_json(sources/'route_gate.json',route)
        for name in ('candidate.json','route_gate.json'):copied[str(sources/name)]=sha256_file(sources/name)
        scientific=[*sorted((archive/'src/redox/qm').glob('*.py')),
                    *sorted((archive/'src/redox/core').glob('*.py')),
                    *sorted((archive/'config').glob('*.py')),*sorted((archive/'config').glob('*.json')),
                    archive/'library/manifest.csv',archive/'scripts/validation/audit/native_partition.py',
                    archive/'scripts/pipeline/run_pi_stacking.py',
                    archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',
                    archive/'environments/cutensor22/runtime.json',origin/'opt.xyz',origin/'result.json',energy]
        plan=dict(status='prepared_exact_molecular_native_state',candidate_alias=candidate_alias,source_id=source_id,
            state=state,charge=raw['charge'],mult=raw['mult'],product_smiles=candidate['product_smiles'],
            product_inchikey=candidate['inchikey'],manifest_row=row,route_gate=route,source_sha256={**copied,
                **{str(p):sha256_file(p) for p in scientific}},source_geometry_sha256=sha256_file(origin/'opt.xyz'),
            source_registry=dict(path=str(Path(registry).resolve()),sha256=sha256_file(registry)),
            target_directory=str(ROOT/'calcs/dft'/source_id/state),
            maximum_host_RSS_GiB=16,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
            prospective_resource_pilot=True,allow_bounded_native_minimum_repair=True,
            global_conformer_convergence_proven=False,material_admission=False,
            source_geometry_scope='Previously optimized exact graft; no new conformer-search convergence inferred')
        atomic_write_json(case/'plan.json',plan)
        spec=dict(status=plan['status'],prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),
            plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=16,minimum_free_GiB=32,
            source_id=source_id,state=state)
        atomic_write_json(folder/'namespace.json',spec)
    result=dict(status='exact_molecular_native_inputs_prepared',candidate_alias=candidate_alias,
        product_inchikey=candidate['inchikey'],namespaces=[str(output/s/'namespace.json') for s in states],
        quantum_jobs_submitted=0,material_admission=False)
    atomic_write_json(output/'prepared.json',result);return result


def verify(path):
    spec=json.loads(Path(path).read_text());validate_layout(spec)
    for key in ('source_id','state'):
        safe_component(spec[key])
    plan=Path(spec['case'])/'plan.json'
    if sha256_file(plan)!=spec['plan_sha256']:raise ValueError('Prepared native-state plan changed')
    data=json.loads(plan.read_text())
    if data['source_id']!=spec['source_id'] or data['state']!=spec['state']:
        raise ValueError('Native namespace and exact state identity differ')
    if any(sha256_file(p)!=h for p,h in data['source_sha256'].items()):raise ValueError('Native source snapshot changed')
    return spec


def command(spec):
    validate_layout(spec);prior=Path(spec['prior_root']);case=Path(spec['case'])
    binary=ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    aliases=[]
    if not spec.get('scientific_sources_read_only',False):
        aliases=['--bind',str(case/'dft'),str(prior/'calcs/dft'),
            '--bind',str(case/'locks'),str(prior/'logs/pipeline_locks'),
            '--bind',str(case/'repairs'),str(prior/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs')]
    argv=[str(binary),'--die-with-parent','--ro-bind','/','/','--ro-bind',spec['archive'],spec['prior_root'],
        '--bind',str(case),str(case),*aliases,
        '--proc','/proc','--dev-bind','/dev','/dev','--tmpfs','/tmp','--chdir',str(prior),'--']
    if spec.get('library_overlay') is not None:
        library=Path(spec['library_overlay'])
        if library.resolve()!=(case/'library').resolve() or not library.is_dir():
            raise ValueError('A library overlay must be the private case/library directory')
        argv[-1:-1]=['--ro-bind',str(library),str(prior/'library')]
    return argv


def probe(path):
    spec=verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code=('import json,pathlib,hashlib,importlib.metadata as m;import ase,rdkit,numpy,psutil;'
        'p=json.loads(pathlib.Path('+repr(spec['case']+'/plan.json')+').read_text());'
        'assert all(hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()==v for k,v in p["source_sha256"].items());'
        'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
        'from redox.qm.native_thermal import implementation_controls;'
        'controls=implementation_controls();assert all(pathlib.Path(c["path"]).is_file() for c in controls);'
        'from redox.core.protocol import load_record;'
        'assert load_record(pathlib.Path(p["target_directory"]),p["charge"],p["mult"]) is not None;'
        'print(json.dumps(dict(versions=wanted,controls=controls,active_energy_readable=True)))')
    run=subprocess.run(command(spec)+[runtime['python'],'-c',code],capture_output=True,text=True,check=True,timeout=180)
    result=dict(status='molecular_native_node_sources_checked',host=os.uname().nodename.split('.')[0],
        namespace_sha256=sha256_file(path),checks=json.loads(run.stdout),GPU_started=False)
    atomic_write_json(Path(path).with_name('node_probe.json'),result);return result


def launch(path,gpu,output):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,gpu_suitability,free_slots,GATE
    from redox.workflows.native_resource_profiles import output_room
    spec=verify(path);case=Path(spec['case']);prior=Path(spec['prior_root']);bound=spec['maximum_host_RSS_GiB']
    proof=json.loads(Path(path).with_name('node_probe.json').read_text());host=os.uname().nodename.split('.')[0]
    if proof['host']!=host or proof['namespace_sha256']!=sha256_file(path):raise ValueError('Exact node-local preflight required')
    if output.exists() or any((case/n).exists() for n in ('status.json','residency_guard','host_guard')):
        raise FileExistsError('Preserve previous scientific and placement attempts')
    output.mkdir(parents=True);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    (output/'launcher_source.py').write_bytes(LAUNCHER_SOURCE_BYTES)
    argv=command(spec);locks=GATE.parent.parent/'.gpu_locks';argv[-1:-1]=['--bind',str(locks),str(locks)]
    argv += [str(GATE),'run',str(gpu),'--label','redox-molecular-native','--',runtime['python'],
        str(prior/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(case/'residency_guard'),
        str(Path(spec['sources'])/'host_guard.py'),'--out',str(case/'host_guard'),'--maximum-rss-gib',str(bound),
        '--minimum-free-gib','32',str(prior/'scripts/fleet/tensor_runtime_worker.py'),'--profile',
        str(prior/'environments/cutensor22/runtime.json'),'--out',str(case/'tensor_runtime'),
        str(Path(spec['sources'])/'worker.py'),'--plan',str(case/'plan.json')]
    env=dict(os.environ,PYTHONPATH=str(prior/'src'),LD_LIBRARY_PATH=runtime['library_path'],
        PYTHONDONTWRITEBYTECODE='1',CUPY_CACHE_DIR='/tmp/cupy_kernel_cache',REDOX_HESSIAN_MEMORY_GIB='6',TMPDIR='/tmp')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
    with allocation_guard():
        capacity=host_capacity(2,minimum_ram_gib=bound+32);promises=remaining_host_promises()
        if not output_room(case)['admit_new_work']:raise RuntimeError('Shared output headroom unavailable; no molecular job started')
        if capacity['available_RAM_bytes']-promises<(bound+32)*2**30:raise RuntimeError('Host RAM promised; no molecular job started')
        cpus=cpu_allocation(2);device=gpu_suitability(gpu)
        if (host,gpu) not in free_slots([host]):raise RuntimeError('GPU occupied; no molecular job started')
        with (output/'command.log').open('x') as log:
            child=subprocess.Popen(['taskset','--cpu-list',','.join(map(str,cpus)),*argv],env=env,
                stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,
                preexec_fn=lambda:os.sched_setaffinity(0,cpus))
        result=dict(status='molecular_native_submitted_start_unverified',host=host,GPU=gpu,pid=child.pid,
            process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            CPU_affinity=cpus,capacity=capacity,remaining_RAM_promises_bytes=promises,device=device,
            namespace=str(path),namespace_sha256=sha256_file(path),argv=argv,started=time.time(),prospective_resource_pilot=True,
            launcher_source=str(output/'launcher_source.py'),launcher_source_sha256=sha256_file(output/'launcher_source.py'))
        atomic_write_json(output/'launch.json',result)
    return child,result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','probe','launch'))
    for name in ('archive','registry','routes','assessment','out','spec'):p.add_argument('--'+name,type=Path)
    p.add_argument('--candidate');p.add_argument('--source-id');p.add_argument('--states');p.add_argument('--gpu',type=int);a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.archive,a.registry,a.routes,a.assessment,a.candidate,a.source_id,a.states.split(','),a.out)
    elif a.action=='probe':result=probe(a.spec.resolve())
    else:
        child,result=launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True)
        result.update(returncode=child.wait(),finished=time.time(),status='molecular_native_returned_independent_review_required')
        atomic_write_json(a.out/'completion.json',result);return
    print(json.dumps(result))


if __name__=='__main__':main()
