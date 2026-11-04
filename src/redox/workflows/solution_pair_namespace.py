"""Source-bound continuation of the existing molecular solvent-relaxation model."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np
from ase.io import read
from rdkit import Chem
from rdkit.Chem import AllChem,rdDetermineBonds,rdMolTransforms

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows.archived_native_namespace import remaining_host_promises


def geometry_pair(paths,smiles=None):
    if smiles is not None:
        from redox.validation.protonation import hydrogen_topology_status
        from redox.build.geometry import geometry_stereo_status
        if any(hydrogen_topology_status(smiles,path)!='ok'
               or geometry_stereo_status(smiles,path) not in ('ok','not_applicable') for path in paths):
            raise ValueError('Stored coordinates do not match the declared full graph/stereochemistry')
    molecules=[];symbols=[]
    for path in paths:
        atoms=read(path);symbols.append(atoms.get_chemical_symbols())
        molecule=Chem.MolFromXYZFile(str(path));rdDetermineBonds.DetermineConnectivity(molecule)
        molecules.append(Chem.RemoveHs(molecule,sanitize=False))
    if symbols[0]!=symbols[1]:raise ValueError('Exact atom ordering and composition required')
    edges=lambda m:{tuple(sorted((b.GetBeginAtomIdx(),b.GetEndAtomIdx()))) for b in m.GetBonds()}
    if edges(molecules[0])!=edges(molecules[1]):raise ValueError('Paired heavy connectivity differs')
    rmsd,_,mapping=AllChem.GetBestAlignmentTransform(molecules[0],molecules[1])
    if any(a!=b for a,b in mapping):
        raise ValueError('Review a nonidentity atom mapping before constructing a shared-AO pair')
    molecule=molecules[0];torsions=[]
    neighbors=lambda index:[a.GetIdx() for a in molecule.GetAtomWithIdx(index).GetNeighbors()]
    for j,k in sorted(edges(molecule)):
        left=[i for i in neighbors(j) if i!=k];right=[i for i in neighbors(k) if i!=j]
        if not left or not right:continue
        visited={j};pending=[j]
        while pending:
            n=pending.pop()
            for a in neighbors(n):
                if {a,n}=={j,k} or a in visited:continue
                visited.add(a);pending.append(a)
        if k in visited:continue
        atoms=[min(left),j,k,min(right)]
        angles=[rdMolTransforms.GetDihedralDeg(m.GetConformer(),*atoms) for m in molecules]
        torsions.append(dict(heavy_atom_indices=atoms,angles_deg=angles,difference_deg=abs((angles[0]-angles[1]+180)%360-180)))
    maximum=max((t['difference_deg'] for t in torsions),default=0.)
    if not np.isfinite([rmsd,maximum]).all() or rmsd>.3 or maximum>30:
        raise ValueError('Stored pair fails the declared 0.3-A / 30-degree same-basin diagnostic')
    return dict(status='paired_geometry_diagnostic_passed',heavy_RMSD_A=rmsd,torsions=torsions,
        maximum_bridge_torsion_difference_deg=maximum,criteria=dict(heavy_RMSD_A=.3,bridge_torsion_difference_deg=30),
        atom_ordering_identical=True,declared_full_graph_checked=smiles is not None,
        global_conformer_convergence_proven=False)


def prepare(archive,review_path,registry_path,route_path,assessment_path,output):
    if output.exists() or output.resolve().is_relative_to(ROOT) or output.resolve().is_relative_to(archive.resolve()):
        raise FileExistsError('Use a fresh durable destination outside both checkouts')
    reader=ArchivedSources(ROOT,archive,relative_to_prior=True)
    review=json.loads(review_path.read_text());reader.verify_all(review['source_sha256'])
    if review['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed':
        raise ValueError('This prospective resource pilot is restricted to the checked TEMPO pair')
    pool=json.loads(registry_path.read_text())['candidates'];candidate=next(r for r in pool if 'tempo' in r['aliases'])
    if Chem.MolToSmiles(Chem.MolFromSmiles(candidate['product_smiles']))!=review['product_smiles']:
        raise ValueError('Registered and checked product graphs differ')
    definition=importlib.util.spec_from_file_location('prospective_solution_pair_route',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(definition);definition.loader.exec_module(gate)
    decision=gate.screen(dict(id='tempo',frag=candidate['declarations'][0]['fragment']),
        json.loads(route_path.read_text())['candidates'],json.loads(assessment_path.read_text()))
    if not decision['proceed_to_full_redox'] or not decision['exact_route_binding_verified']:
        raise ValueError('Exact TEMPO route gate did not pass')
    paths=[reader.path(ROOT/'calcs/dft/tempo'/state/'opt.xyz') for state in ('ox','rad')]
    evidence=geometry_pair(paths,review['product_smiles'])
    output.mkdir(parents=True);case=output/'case';case.mkdir();sources=output/'sources';sources.mkdir()
    copied={}
    for source,name in [(ROOT/'scripts/pipeline/run_checked_solution_pair.py','worker.py'),
                        (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
                        (Path(__file__),'launcher_source.py')]:
        target=sources/name;shutil.copyfile(source,target);copied[str(target)]=sha256_file(target)
    pair=dict(states=[dict(name=state,charge=q,mult=m,geometry=str(ROOT/'calcs/dft/tempo'/state/'opt.xyz'),
        state_dir=str(ROOT/'calcs/dft/tempo'/state)) for state,q,m in [('ox',1,1),('rad',0,2)]],
        require_native_minimum_states=['ox','rad'],conformer_policy='same_basin',implementation_controls={})
    for name in ('native_smd_independent_reference.json','native_smd_gpu_reference.json'):
        path=ROOT/'results/validation/physical_reorg_2026-10-11'/name
        pair['implementation_controls'][str(path)]=sha256_file(reader.path(path))
    for name,path in [('native_review.json',review_path),('registry.json',registry_path),
                      ('route_evidence.json',route_path),('assessment.json',assessment_path)]:
        target=sources/name;shutil.copyfile(path,target);copied[str(target)]=sha256_file(target)
    source_files=['scripts/validation/solution_reorg_pair.py','scripts/validation/solution_reorg_campaign.py',
        'scripts/validation/physical_reorg.py','scripts/validation/audit_solution_reorg_pair.py',
        'scripts/pipeline/run_pi_stacking.py','src/redox/qm/nonequilibrium_smd.py',
        'src/redox/qm/nonequilibrium_pcm.py','src/redox/properties/solution_reorg.py',
        'src/redox/validation/solution_relaxation.py','src/redox/core/protocol.py','src/redox/qm/dft.py',
        'scripts/fleet/exclusive_python_worker.py','scripts/fleet/bounded_host_python_worker.py',
        'scripts/fleet/tensor_runtime_worker.py','environments/cutensor22/runtime.json']
    hashes={**copied,**review['source_sha256'],**pair['implementation_controls'],
        **{str(ROOT/p):sha256_file(reader.path(ROOT/p)) for p in source_files}}
    plan=dict(status='prepared_exact_TEMPO_solution_pair_resource_pilot',candidate_id='tempo',
        product_smiles=review['product_smiles'],product_inchikey=review['product_inchikey'],pair=pair,
        route_gate=decision,geometry_diagnostic=evidence,source_sha256=hashes,
        maximum_host_RSS_GiB=16,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
        prospective_resource_pilot=True,descriptor='directional frozen-slow SMD electronic solution relaxation',
        gas_inner_or_Born_outer_terms_added=False,material_admission=False,source_equations_changed=False)
    atomic_write_json(case/'plan.json',plan);reader.recheck()
    spec=dict(status=plan['status'],prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),
              plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=16,minimum_free_GiB=32)
    atomic_write_json(output/'namespace.json',spec);return spec


def verify(path):
    spec=json.loads(path.read_text());validate_layout(spec);plan=Path(spec['case'])/'plan.json'
    if sha256_file(plan)!=spec['plan_sha256']:raise ValueError('Prepared solution pair changed')
    reader=ArchivedSources(spec['prior_root'],spec['archive'],relative_to_prior=True)
    reader.verify_all(json.loads(plan.read_text())['source_sha256']);reader.recheck();return spec


def validate_layout(spec):
    for name in ('prior_root','archive','case','sources'):
        value=Path(spec[name])
        if not value.is_absolute() or '..' in value.parts or value==Path('/'):
            raise ValueError('Explicit absolute non-root namespace paths required')
    for name in ('case','sources'):
        value=Path(spec[name]).resolve()
        if any(value.is_relative_to(Path(spec[root]).resolve()) for root in ('prior_root','archive')):
            raise ValueError('Writable solution outputs must remain outside both source checkouts')
    if Path(spec['case']).parent!=Path(spec['sources']).parent:
        raise ValueError('One declared solution preparation root required')


def command(spec):
    validate_layout(spec)
    binary=ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    return [str(binary),'--die-with-parent','--ro-bind','/','/','--ro-bind',spec['archive'],spec['prior_root'],
        '--bind',spec['case'],spec['case'],'--proc','/proc','--dev-bind','/dev','/dev','--tmpfs','/tmp',
        '--chdir',spec['prior_root'],'--']


def probe(path):
    spec=verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code=('import json,pathlib,hashlib,importlib.metadata as m;import ase,rdkit,numpy,psutil; '
          'p=json.loads(pathlib.Path('+repr(spec['case']+'/plan.json')+').read_text());'
          'assert all(hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()==v for k,v in p["source_sha256"].items());'
          'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
          'from redox.qm.native_thermal import load_native_thermal;'
          'assert all(load_native_thermal(pathlib.Path(s["state_dir"]),s["charge"],s["mult"]) is not None for s in p["pair"]["states"]);'
          'print(json.dumps(dict(versions=wanted,both_native_minima_readable=True)))')
    result=subprocess.run(command(spec)+[runtime['python'],'-c',code],capture_output=True,text=True,check=True,timeout=180)
    proof=dict(status='solution_pair_node_sources_and_minima_checked',host=os.uname().nodename.split('.')[0],
               namespace_sha256=sha256_file(path),checks=json.loads(result.stdout),GPU_started=False)
    atomic_write_json(path.with_name('node_probe.json'),proof);return proof


def launch(path,gpu,output):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,gpu_suitability,free_slots,GATE
    spec=verify(path);proof=json.loads(path.with_name('node_probe.json').read_text());host=os.uname().nodename.split('.')[0]
    if proof['host']!=host or proof['namespace_sha256']!=sha256_file(path):raise ValueError('Exact node preflight required')
    case=Path(spec['case']);prior=Path(spec['prior_root']);bound=spec['maximum_host_RSS_GiB']
    if output.exists() or any((case/p).exists() for p in ('residency_guard','host_guard','status.json')):
        raise FileExistsError('Preserve prior placements and scientific attempts')
    output.mkdir(parents=True);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    argv=command(spec);locks=GATE.parent.parent/'.gpu_locks';argv[-1:-1]=['--bind',str(locks),str(locks)]
    argv += [str(GATE),'run',str(gpu),'--label','redox-solution-pair','--',runtime['python'],
        str(prior/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(case/'residency_guard'),
        str(Path(spec['sources'])/'host_guard.py'),'--out',str(case/'host_guard'),
        '--maximum-rss-gib',str(bound),'--minimum-free-gib','32',
        str(prior/'scripts/fleet/tensor_runtime_worker.py'),'--profile',str(prior/'environments/cutensor22/runtime.json'),
        '--out',str(case/'tensor_runtime'),str(Path(spec['sources'])/'worker.py'),'--plan',str(case/'plan.json')]
    env=dict(os.environ,PYTHONPATH=str(prior/'src'),LD_LIBRARY_PATH=runtime['library_path'],PYTHONDONTWRITEBYTECODE='1',CUPY_CACHE_DIR='/tmp/cupy_kernel_cache')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
    with allocation_guard():
        capacity=host_capacity(2,minimum_ram_gib=bound+32);promises=remaining_host_promises()
        if capacity['available_RAM_bytes']-promises<(bound+32)*2**30:raise RuntimeError('Host RAM promised; no solution calculation started')
        cpus=cpu_allocation(2);device=gpu_suitability(gpu)
        if (host,gpu) not in free_slots([host]):raise RuntimeError('GPU occupied; no solution calculation started')
        with (output/'command.log').open('x') as log:
            child=subprocess.Popen(['taskset','--cpu-list',','.join(map(str,cpus)),*argv],env=env,
                stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True,
                preexec_fn=lambda:os.sched_setaffinity(0,cpus))
        record=dict(status='solution_pair_submitted_start_unverified',host=host,GPU=gpu,pid=child.pid,
            process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            CPU_affinity=cpus,capacity=capacity,remaining_RAM_promises_bytes=promises,device=device,
            namespace=str(path),namespace_sha256=sha256_file(path),argv=argv,started=time.time(),prospective_resource_pilot=True)
        atomic_write_json(output/'launch.json',record)
    return child,record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','probe','launch'))
    for name in ('archive','review','registry','routes','assessment','out','spec'):p.add_argument('--'+name,type=Path)
    p.add_argument('--gpu',type=int);a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():
            host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.archive.resolve(),a.review.resolve(),a.registry.resolve(),a.routes.resolve(),a.assessment.resolve(),a.out.resolve())
    elif a.action=='probe':result=probe(a.spec.resolve())
    else:
        child,result=launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True)
        result.update(returncode=child.wait(),finished=time.time(),status='solution_pair_returned_independent_review_required')
        atomic_write_json(a.out/'completion.json',result);return
    print(json.dumps(result))


if __name__=='__main__':main()
