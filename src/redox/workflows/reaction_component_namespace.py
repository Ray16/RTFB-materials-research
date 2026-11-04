"""Prepare exact route reagents without adding them to the candidate registry.

Each fixed conformer has a private archive view, a one-row reagent manifest,
durable checkpoints and independent native thermochemistry. UMA supplies only
the warm start; all reported energies and derivatives remain DFT/SMD.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows import molecular_native_namespace as molecular


def validate_component(recipe, route):
    if recipe.get('schema')!=1 or recipe.get('role')!='reaction_component':
        raise ValueError('An explicit reaction-component recipe is required')
    for key in ('source_id','state'):molecular.safe_component(recipe[key])
    binding=route['route_binding'];mol=Chem.MolFromSmiles(recipe['smiles'])
    if (mol is None or len(Chem.GetMolFrags(mol))!=1
            or Chem.MolToSmiles(mol)!=Chem.MolToSmiles(Chem.MolFromSmiles(binding[recipe['route_component']]))
            or type(recipe['charge']) is not int or Chem.GetFormalCharge(mol)!=recipe['charge']
            or type(recipe['mult']) is not int or recipe['mult']<1
            or (sum(a.GetAtomicNum() for a in Chem.AddHs(mol).GetAtoms())-recipe['charge']-(recipe['mult']-1))%2):
        raise ValueError('Reagent graph, charge or spin differs from the exact route')
    if not 1<=recipe['top_k']<=20 or not 0<recipe['UMA_fmax_eV_A']<=.03:
        raise ValueError('Explicit bounded conformer and warm-start settings required')
    return mol


def prepare(archive, recipe_path, registry_path, route_path, assessment_path, output, *, conformer_indices=None, uma_reference=None):
    archive=Path(archive).resolve();output=Path(output).resolve()
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Use a fresh durable destination outside both checkouts')
    recipe=json.loads(recipe_path.read_text());routes=json.loads(route_path.read_text())['candidates']
    route=routes[recipe['route_candidate']];mol=validate_component(recipe,route)
    records=json.loads(registry_path.read_text())['candidates']
    candidates=[c for c in records if recipe['route_candidate'] in c['aliases']]
    if len(candidates)!=1:raise ValueError('One exact registered candidate is required')
    candidate=candidates[0];declaration=next(d for d in candidate['declarations'] if d['id']==recipe['route_candidate'])
    screen=molecular.load(ROOT/'scripts/pipeline/screen_grafting_routes.py','component_route_screen')
    gate=screen.screen(dict(id=recipe['route_candidate'],frag=declaration['fragment'],
        product_smiles=candidate['product_smiles']),routes,json.loads(assessment_path.read_text()))
    if not gate['proceed_to_full_redox'] or not gate.get('exact_route_binding_verified'):
        raise ValueError('The linked exact candidate route must pass prospectively')
    from redox.build.structures import conformer_ensemble,to_xyz
    from redox.build.manifest import COLS
    molecule,ids,forcefield,energies=conformer_ensemble(mol,recipe['charge'],recipe['top_k'],recipe['conformer_seed'])
    selected=list(range(len(ids))) if conformer_indices is None else conformer_indices
    if not selected or len(selected)!=len(set(selected)) or any(i not in range(len(ids)) for i in selected):
        raise ValueError('Select distinct existing conformer indices')
    if uma_reference is not None:
        reference=json.loads(uma_reference.read_text())
        if (reference.get('status')!='fixed_geometry_standard_UMA_V100_reference'
                or reference['smiles']!=recipe['smiles'] or reference['charge']!=recipe['charge'] or reference['mult']!=recipe['mult']
                or any(sha256_file(k)!=v for k,v in reference['source_sha256'].items())):
            raise ValueError('A source-bound same-reagent V100 reference is required')
    output.mkdir(parents=True);namespaces=[]
    for index,(cid,energy) in enumerate(zip(ids,energies)):
        if index not in selected:continue
        folder=output/f'conf{index:03d}';case=folder/'case';sources=folder/'sources'
        sources.mkdir(parents=True);case.mkdir()
        for name in ('dft','locks','repairs','library'):(case/name).mkdir()
        for source,name in [(ROOT/'scripts/pipeline/run_reaction_component.py','worker.py'),
                            (ROOT/'scripts/validation/review_reaction_component.py','reviewer.py'),
                            (ROOT/'scripts/validation/review_molecular_native.py','native_review.py'),
                            (ROOT/'scripts/pipeline/run_molecular_native.py','native_worker.py'),
                            (ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
                            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),
                            (Path(__file__),'preparer.py'),
                            (ROOT/'src/redox/workflows/molecular_native_namespace.py','launcher.py'),
                            (ROOT/'scripts/pipeline/screen_grafting_routes.py','route_screen.py'),
                            (ROOT/'src/redox/build/structures.py','geometry_builder.py'),
                            (recipe_path,'recipe_input.json'),(registry_path,'registry.json'),
                            (route_path,'route_evidence.json'),(assessment_path,'assessment.json')]:shutil.copyfile(source,sources/name)
        single=Chem.Mol(molecule);single.RemoveAllConformers();single.AddConformer(molecule.GetConformer(cid),assignId=True)
        geometry=case/'input.xyz';geometry.write_text(to_xyz(single,f'route reagent; seed={recipe["conformer_seed"]}; conf={index}; FF={forcefield}; FF_E={energy}'))
        row=dict(id=recipe['source_id'],state=recipe['state'],name=recipe['name'],family='reaction_component',
            charge=recipe['charge'],mult_hint=recipe['mult'],n_e=0,n_conf=1,solv_preopt=1,smiles=recipe['smiles'])
        with (case/'library/manifest.csv').open('w') as stream:
            writer=csv.DictWriter(stream,fieldnames=COLS);writer.writeheader();writer.writerow(row)
        atomic_write_json(sources/'recipe.json',recipe);atomic_write_json(sources/'route_gate.json',gate)
        if uma_reference is not None:shutil.copyfile(uma_reference,sources/'UMA_V100_reference.json')
        scientific=[*sorted((archive/'src/redox/qm').glob('*.py')),*sorted((archive/'src/redox/core').glob('*.py')),
            *sorted((archive/'config').glob('*.py')),*sorted((archive/'config').glob('*.json')),
            archive/'scripts/validation/audit/native_partition.py',archive/'scripts/pipeline/run_pi_stacking.py',
            *sorted((archive/'scripts/fleet').glob('*.py')),archive/'environments/cutensor22/runtime.json']
        input_paths=[geometry,case/'library/manifest.csv',*sorted(sources.glob('*')),*scientific]
        plan=dict(schema=1,role='reaction_component',source_id=recipe['source_id'],state=recipe['state'],
            charge=recipe['charge'],mult=recipe['mult'],smiles=recipe['smiles'],inchikey=Chem.MolToInchiKey(mol),
            route_candidate=recipe['route_candidate'],route_gate=gate,recipe=recipe,conformer_index=index,
            forcefield=forcefield,forcefield_energy_kcal_mol=energy,initial_geometry_sha256=sha256_file(geometry),
            source_sha256={str(p):sha256_file(p) for p in input_paths},
            target_directory=str(ROOT/'calcs/dft'/recipe['source_id']/recipe['state']),
            maximum_host_RSS_GiB=24,minimum_free_GiB=32,maximum_peak_GPU_MiB=31744,
            material_admission=False,global_conformer_convergence_proven=False)
        if uma_reference is not None:plan['UMA_cross_device_reference']=str(sources/'UMA_V100_reference.json')
        atomic_write_json(case/'plan.json',plan)
        spec=dict(status='prepared_reaction_component',prior_root=str(ROOT),archive=str(archive),case=str(case),
            sources=str(sources),library_overlay=str(case/'library'),plan_sha256=sha256_file(case/'plan.json'),
            maximum_host_RSS_GiB=24,minimum_free_GiB=32,source_id=recipe['source_id'],state=recipe['state'])
        atomic_write_json(folder/'namespace.json',spec);namespaces.append(str(folder/'namespace.json'))
    result=dict(status='route_reagent_conformers_prepared',role='reaction_component',namespaces=namespaces,
        actual_conformers=len(namespaces),embedded_conformers=len(ids),selected_conformer_indices=selected,
        requested_conformers=recipe['top_k'],quantum_jobs_submitted=0,
        source_id=recipe['source_id'],route_candidate=recipe['route_candidate'],material_admission=False)
    atomic_write_json(output/'prepared.json',result);return result


def probe(path):
    spec=molecular.verify(path);runtime=json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code=('import json,pathlib,importlib.metadata as m; import rdkit,ase,psutil;'
        'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
        'from redox.qm.native_thermal import implementation_controls;'
        'assert all(pathlib.Path(c["path"]).is_file() for c in implementation_controls());'
        'from redox.qm.uma_model_identity import model_identity;'
        'p=json.loads(pathlib.Path("config/uma_execution.json").read_text());i=model_identity(p["model"]);'
        'assert i["checkpoint_sha256"]==p["checkpoint_sha256"] and i["reference_sha256"]==p["reference_sha256"];'
        'assert p["model"] in json.loads(pathlib.Path(i["registry_path"]).read_text());'
        'print(json.dumps(dict(versions=wanted,UMA_identity=i,GPU_started=False)))')
    env=dict(os.environ,FAIRCHEM_CACHE_DIR=str(ROOT.parent/'.cache/fairchem'),PYTHONDONTWRITEBYTECODE='1')
    run=subprocess.run(molecular.command(spec)+[runtime['python'],'-c',code],env=env,capture_output=True,text=True,check=True,timeout=180)
    result=dict(status='reaction_component_node_sources_checked',host=os.uname().nodename.split('.')[0],
        namespace_sha256=sha256_file(path),checks=json.loads(run.stdout),GPU_started=False)
    atomic_write_json(path.with_name('node_probe.json'),result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','probe','launch','submit'))
    for name in ('archive','recipe','registry','routes','assessment','out','spec'):p.add_argument('--'+name,type=Path)
    p.add_argument('--gpu',type=int);p.add_argument('--indices');p.add_argument('--uma-reference',type=Path);a=p.parse_args()
    if a.action not in ('launch','submit'):
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':result=prepare(a.archive,a.recipe,a.registry,a.routes,a.assessment,a.out,
        conformer_indices=None if a.indices is None else [int(i) for i in a.indices.split(',')],uma_reference=a.uma_reference)
    elif a.action=='probe':result=probe(a.spec.resolve())
    elif a.action=='submit':
        molecular.verify(a.spec.resolve());a.out.parent.mkdir(parents=True,exist_ok=True)
        submission=a.out.with_name(a.out.name+'_submission.json');log=a.out.with_name(a.out.name+'_controller.log')
        if submission.exists() or log.exists() or a.out.exists():raise FileExistsError('Preserve existing submission attempts')
        argv=[sys.executable,'-m',__name__ if __name__!='__main__' else 'redox.workflows.reaction_component_namespace',
            'launch','--spec',str(a.spec.resolve()),'--gpu',str(a.gpu),'--out',str(a.out.resolve())]
        with log.open('x') as stream:
            child=subprocess.Popen(argv,stdout=stream,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
        result=dict(status='reaction_component_submission_start_unverified',host=os.uname().nodename,pid=child.pid,
            process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
            argv=argv,submitted=time.time(),log=str(log.resolve()))
        atomic_write_json(submission,result)
    else:
        child,result=molecular.launch(a.spec.resolve(),a.gpu,a.out.resolve());print(json.dumps(result),flush=True)
        result.update(returncode=child.wait(),status='reaction_component_returned_review_required')
        atomic_write_json(a.out/'completion.json',result);return
    print(json.dumps(result))


if __name__=='__main__':main()
