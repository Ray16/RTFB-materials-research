"""Prospective, route-bound charged-polymer SCF/gradient resource pilots."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.workflows import molecular_native_namespace as execution


def resting_state(case,allow_empty_reference=False):
    graph=Chem.MolFromSmiles(case['smiles']);charge=Chem.GetFormalCharge(graph)
    states=[s for s in case['proposed_global_states'] if s['charge']==case['charge']]
    empty=allow_empty_reference and not case['occupied_sites'] and charge==case['charge']==0
    if len(states)!=1 or states[0]['mult']!=1 or charge!=case['charge'] or (charge<=0 and not empty):
        raise ValueError('A closed-shell ionic graph or explicitly scoped empty reference is required')
    return states[0],empty


def prepare(archive,preparation,output,key=None,execution_profile='compact_28GiB',allow_empty_reference=False):
    if execution_profile not in ('compact_28GiB','direct_large_GPU_reference'):raise ValueError('Explicit measured execution profile required')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Use a new durable pilot outside both scientific checkouts')
    audit=execution.load(ROOT/'scripts/validation/check_charged_polymer_preparation.py','charged_pilot_initial_geometry_audit')
    checked=audit.check(preparation)
    data=json.loads(preparation.read_text());records=data['cases']
    if key is None:
        selected=max(records,key=lambda r:(Chem.AddHs(Chem.MolFromSmiles(r['case']['smiles'])).GetNumAtoms(),r['case']['key']))
    else:
        matches=[r for r in records if r['case']['key']==key]
        if len(matches)!=1:raise ValueError('One exact prepared graph/configuration required')
        selected=matches[0]
    state=selected['case'];graph=Chem.AddHs(Chem.MolFromSmiles(state['smiles']))
    selected_state,empty_reference=resting_state(state,allow_empty_reference)
    source=preparation.parent/'source_snapshot';registry=json.loads((source/'1_registry.json').read_text())['candidates']
    candidates=[r for r in registry if data['candidate_alias'] in r['aliases']]
    if len(candidates)!=1:raise ValueError('One exact candidate graph required')
    gate=execution.load(ROOT/'scripts/pipeline/screen_grafting_routes.py','charged_pilot_prospective_route')
    route=gate.screen(dict(id=data['candidate_alias'],frag=data['fragment'],product_smiles=candidates[0]['product_smiles']),
        json.loads((source/'2_route_evidence.json').read_text())['candidates'],json.loads((source/'3_assessment.json').read_text()))
    if not route['proceed_to_full_redox'] or not route.get('exact_route_binding_verified'):
        raise ValueError('The exact prospective attachment route must pass before this pilot')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir();sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    files=[(ROOT/'scripts/validation/run_charged_polymer_pilot.py','worker.py'),
        (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),(ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
        (ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
        (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(Path(__file__),'preparer.py'),
        (Path(selected['initial_geometry']),'input.xyz'),(preparation,'preparation.json'),
        (ROOT/'scripts/pipeline/screen_grafting_routes.py','route_gate.py'),
        (ROOT/'scripts/validation/check_charged_polymer_preparation.py','geometry_auditor.py')]
    for original,name in files:shutil.copyfile(original,sources/name)
    atomic_write_json(sources/'initial_audit.json',checked);atomic_write_json(sources/'prospective_route.json',route)
    spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),
        source_id='charged_polymer_resource_pilot',state='rest',maximum_host_RSS_GiB=40,minimum_free_GiB=32,
        scientific_sources_read_only=True)
    runtime=json.loads((archive/'environments/cutensor22/runtime.json').read_text())
    program=('import json,importlib.util,pathlib;from redox.core.common import ROOT;'
        's=importlib.util.spec_from_file_location("pilot_CDS_witness",ROOT/"scripts/pipeline/smooth_cds_control.py");'
        'm=importlib.util.module_from_spec(s);s.loader.exec_module(m);'
        'print(json.dumps(dict(model=m.MODEL,evidence=m.component_evidence(ROOT/"results/validation/smooth_CDS_2026-10-25/pilot"))))')
    result=subprocess.run(execution.command(spec)+[runtime['python'],'-c',program],capture_output=True,text=True,timeout=180,check=True)
    model=json.loads(result.stdout);reader=ArchivedSources(ROOT,archive);reader.verify_all(model['evidence'])
    package=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
    bound=[*sources.glob('*'),archive/'src/redox/qm/dft.py',archive/'scripts/pipeline/smooth_cds_control.py',
        archive/'scripts/pipeline/run_pi_stacking.py',archive/'src/redox/workflows/resources.py',
        archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',
        archive/'environments/cutensor22/runtime.json',package/'dft/gen_grid.py',package/'dft/numint.py',package/'lib/libgdft.so']
    bound.extend((archive/'src/redox/core').glob('*.py'))
    if execution_profile=='compact_28GiB':bound.extend((archive/'scripts/pipeline').glob('*.py'))
    hashes={**{str(reader.path(p)):h for p,h in model['evidence'].items()},**{str(p):sha256_file(p) for p in bound}}
    plan=dict(role='charged_polymer_SCF_gradient_resource_pilot',charged_polymer_pilot=True,
        matched_empty_polymer_reference=empty_reference,
        source_id=spec['source_id'],state='rest',source_sha256=hashes,exact_candidate_alias=data['candidate_alias'],
        preparation_source=dict(path=str(preparation),sha256=sha256_file(preparation)),
        prepared_case=state,selected_initial_geometry_sha256=selected['initial_geometry_sha256'],
        atoms=graph.GetNumAtoms(),charge=state['charge'],mult=selected_state['mult'],selected_global_state=selected_state,
        solvent_model=model['model'],calculation_ion_convention=data['config']['ionic_model']['calculation_convention'],
        explicit_counterions_in_coordinates=False,requested_medium='acetonitrile',
        quantum_settings=dict(xc='r2scan',basis='def2-svp',disp='d4',solvent='acetonitrile',grid_level=3,
            SCF_conv_tol=1e-12,SCF_conv_tol_grad=1e-8,SCF_max_cycle=200,gradient_grid_response=True),
        execution_profile=execution_profile,maximum_host_RSS_GiB=40,minimum_free_GiB=32,
        maximum_peak_GPU_MiB=31744 if execution_profile=='compact_28GiB' else 88000,
        minimum_GPU_memory_MiB=31*1024 if execution_profile=='compact_28GiB' else 80*1024,
        maximum_wall_seconds=7200,maximum_durable_output_bytes=512*2**20,
        minimum_initial_storage_GiB=32.5,minimum_live_storage_GiB=32,
        geometry_optimization_performed=False,native_G_qualified=False,material_admission=False)
    atomic_write_json(case/'plan.json',plan);spec['plan_sha256']=sha256_file(case/'plan.json')
    atomic_write_json(output/'namespace.json',spec);reader.recheck()
    return dict(status='exact_route_bound_charged_polymer_resource_pilot_prepared',namespace=str(output/'namespace.json'),
                selected_case=state['key'],charge=state['charge'],mult=selected_state['mult'],atoms=graph.GetNumAtoms(),quantum_jobs_submitted=0)


def prepare_optimization(archive,preparation,agreement_path,output,key=None,allow_empty_reference=False):
    agreement=json.loads(agreement_path.read_text())
    if agreement['status']!='charged_polymer_execution_profiles_compared' or not agreement['initial_state_execution_agreement_passed']:
        raise ValueError('Paired direct/compact ionic energy/force agreement required before optimization')
    for path,digest in agreement['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Ionic execution qualification changed')
    result=prepare(archive,preparation,output,key,execution_profile='compact_28GiB',allow_empty_reference=allow_empty_reference)
    spec=json.loads((output/'namespace.json').read_text());case=Path(spec['case']);sources=Path(spec['sources'])
    shutil.copyfile(ROOT/'scripts/validation/run_charged_polymer_optimization.py',sources/'worker.py')
    shutil.copyfile(agreement_path,sources/'execution_agreement.json')
    plan=json.loads((case/'plan.json').read_text());plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.glob('*')})
    plan.update(role='charged_polymer_solution_geometry_pilot',charged_polymer_optimization=True,
        maximum_wall_seconds=28800,maximum_durable_output_bytes=256*2**20,minimum_initial_storage_GiB=32.25,
        optimizer_settings=dict(maxsteps=200,convergence_energy=1e-6,convergence_grms=3e-4,convergence_gmax=4.5e-4,
            convergence_drms=1.2e-3,convergence_dmax=1.8e-3),
        checkpoint_policy='Every geometry/energy/gradient JSON retained; current and final SCF orbitals retained as working checkpoints',
        native_G_qualified=False)
    atomic_write_json(case/'plan.json',plan);spec['plan_sha256']=sha256_file(case/'plan.json');atomic_write_json(output/'namespace.json',spec)
    result['status']='exact_route_bound_ionic_solution_geometry_pilot_prepared';return result


def prepare_finalization(source_namespace,output):
    source=execution.verify(source_namespace);source_case=Path(source['case']);old=json.loads((source_case/'plan.json').read_text())
    stage=json.loads((source_case/'status.json').read_text());returned=json.loads((source_case/'optimizer_return.json').read_text())
    guard=json.loads((source_case/'residency_guard/owner.json').read_text());completions=list(source_case.parent.glob('placement*/completion.json'))
    owners=execution.load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','ionic_final_closed_owner')
    if (not old.get('charged_polymer_optimization') or stage['status']!='charged_polymer_solution_optimization_failed_preserved'
            or 'cudaErrorMemoryAllocation' not in stage.get('error','') or not returned['converged']
            or not guard.get('guard_stopped') or len(completions)!=1 or owners.alive(stage['owner']) is not False
            or any(guard.get(k)!=stage['owner'].get(k) for k in ('host','pid','process_start'))):
        raise ValueError('A proven closed resource failure after converged optimizer return is required')
    completion=json.loads(completions[0].read_text())
    if completion['returncode'] in (None,0):raise ValueError('Original failed placement must be closed')
    initial_audit=json.loads((Path(source['sources'])/'initial_audit.json').read_text())
    preparation=Path(initial_audit['source'])
    if sha256_file(preparation)!=initial_audit['source_sha256']:raise ValueError('Original preparation manifest changed')
    result=prepare(Path(source['archive']),preparation,output,old['prepared_case']['key'])
    spec=json.loads((output/'namespace.json').read_text());case=Path(spec['case']);sources=Path(spec['sources'])
    if len(returned['elements'])!=old['atoms'] or len(returned['coordinates_A'])!=old['atoms']:
        raise ValueError('Complete optimizer-returned coordinates required')
    (sources/'input.xyz').write_text(str(old['atoms'])+'\nExact optimizer-returned coordinates\n'+'\n'.join(
        symbol+' '+' '.join(format(float(v),'.17g') for v in xyz)
        for symbol,xyz in zip(returned['elements'],returned['coordinates_A']))+'\n')
    geometry=execution.load(ROOT/'scripts/validation/review_charged_polymer_optimization.py','returned_ionic_connectivity')
    atomic_write_json(sources/'optimizer_return_geometry_review.json',geometry.graph_checks(old['prepared_case']['smiles'],Path(source['sources'])/'input.xyz',sources/'input.xyz'))
    for path,name in [(source_case/'optimizer_return.json','source_optimizer_return.json'),(source_case/'status.json','source_optimization_status.json')]:shutil.copyfile(path,sources/name)
    plan=json.loads((case/'plan.json').read_text());paths=[source_namespace,source_case/'plan.json',source_case/'status.json',
        source_case/'optimizer_return.json',source_case/'returned.xyz',source_case/'residency_guard/owner.json',completions[0],
        *sorted((source_case/'geometry_evaluations').glob('*.json'))]
    source_hashes={str(p):sha256_file(p) for p in paths}
    plan['source_sha256'].update(old['source_sha256']);plan['source_sha256'].update(source_hashes)
    plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.glob('*')})
    plan.update(optimizer_returned_input=True,input_geometry_role='converged_optimizer_return_fresh_final_state_evaluation',
        original_initial_geometry_sha256=plan['selected_initial_geometry_sha256'],selected_initial_geometry_sha256=sha256_file(sources/'input.xyz'),
        source_optimization=dict(namespace=str(source_namespace),source_case=str(source_case),source_sha256=source_hashes),native_G_qualified=False)
    atomic_write_json(case/'plan.json',plan);spec['plan_sha256']=sha256_file(case/'plan.json');atomic_write_json(output/'namespace.json',spec)
    result['status']='converged_ionic_optimizer_return_prepared_for_fresh_final_evaluation';return result


def prepare_native_component(review_path,output):
    """One closed-shell ionic native pilot, at a separately checked geometry."""
    from ase.io import read
    from redox.build.geometry import _skeleton
    from redox.build.oligomer_topology import chain_topology
    from rdkit.Chem import rdDetermineBonds
    from redox.qm.hessian_restart import response_ceiling_gib
    report=json.loads(review_path.read_text())
    if report['status']!='charged_polymer_solution_stationary_geometry_independently_checked' or report['native_G_qualified']:
        raise ValueError('A separately checked stationary ionic graph is required')
    for path,digest in report['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Stationary source changed')
    source=execution.verify(Path(report['namespace']));source_case=Path(source['case']);old=json.loads((source_case/'plan.json').read_text())
    if old['mult']!=1 or old['charge']<=0:raise ValueError('This first charged native pilot is a closed-shell cation')
    evaluation=json.loads((source_case/'evaluation.json').read_text());nao=evaluation['nao']
    from pyscf import gto
    atoms=read(Path(source['sources'])/'input.xyz');symbols=atoms.get_chemical_symbols()
    mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),unit='Angstrom',basis='def2-svp',charge=old['charge'],spin=0,verbose=0)
    nocc=mol.nelectron//2
    if nao!=mol.nao or nao>1193 or nocc>235:raise ValueError('First charged native pilot must fit inside the checked neutral dimension envelope')
    initial=json.loads((Path(source['sources'])/'initial_audit.json').read_text());preparation=Path(initial['source'])
    if sha256_file(preparation)!=initial['source_sha256']:raise ValueError('Original ionic preparation changed')
    result=prepare(Path(source['archive']),preparation,output,old['prepared_case']['key'])
    spec=json.loads((output/'namespace.json').read_text());case=Path(spec['case']);sources=Path(spec['sources'])
    shutil.copyfile(Path(source['sources'])/'input.xyz',sources/'input.xyz')
    for original,name in [(ROOT/'scripts/validation/run_screening_polymer_hessian.py','worker.py'),
            (ROOT/'src/redox/qm/response_diagnostic.py','response_helper.py'),
            (ROOT/'src/redox/qm/hessian_restart.py','restart_helper.py'),
            (ROOT/'src/redox/qm/response_resume.py','response_resume_helper.py')]:shutil.copyfile(original,sources/name)
    native=execution.load(Path(source['archive'])/'scripts/pipeline/smooth_polymer_native.py','ionic_native_unchanged_settings')
    graph=Chem.AddHs(Chem.MolFromSmiles(old['prepared_case']['smiles']));observed=Chem.MolFromXYZFile(str(sources/'input.xyz'))
    rdDetermineBonds.DetermineConnectivity(observed,charge=0);mapping=_skeleton(observed).GetSubstructMatch(_skeleton(graph))
    if len(mapping)!=graph.GetNumAtoms():raise ValueError('Full ionic native input connectivity differs')
    topology=chain_topology(old['prepared_case']['smiles'],old['prepared_case']['repeat_units'],old['prepared_case']['reactive_sites'])
    site=old['prepared_case']['occupied_sites'][0];benzyl=topology['para_graft_sites'][site-1]['benzyl_atom']
    nitrogen=next(a.GetIdx() for a in graph.GetAtomWithIdx(benzyl).GetNeighbors() if a.GetAtomicNum()==7 and a.GetFormalCharge()==1)
    extras=[[mapping[nitrogen],2],[mapping[benzyl],2]]
    cap,minimum=response_ceiling_gib(nao,nocc)
    plan=json.loads((case/'plan.json').read_text());plan.update(role='charged_polymer_native_component',charged_polymer_pilot=False,
        charged_polymer_native_component=True,screening_variant='corrected',symbols=symbols,
        original_namespace=report['namespace'],original_namespace_sha256=sha256_file(Path(report['namespace'])),
        source_stationary_review=str(review_path),points=[dict(evaluation=str(source_case/'evaluation.json'))],
        numerical_settings=native.NATIVE_SETTINGS,extra_FD_columns=extras,response_available_memory_GiB=cap,
        response_nao=nao,response_nocc=nocc,minimum_aligned_response_bytes=minimum,
        maximum_host_RSS_GiB=80,maximum_wall_seconds=28800,maximum_durable_output_bytes=4*2**30,
        minimum_initial_storage_GiB=36,minimum_live_storage_GiB=32,
        selected_initial_geometry_sha256=sha256_file(sources/'input.xyz'),native_G_qualified=False,material_admission=False)
    plan['source_sha256'].update(report['source_sha256']);plan['source_sha256'][str(review_path)]=sha256_file(review_path)
    for name in ('qrrho.py','vibrations.py','cds_finite.py','hessian_memory.py','native_thermal.py'):
        path=Path(source['archive'])/'src/redox/qm'/name;plan['source_sha256'][str(path)]=sha256_file(path)
    plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.glob('*')})
    atomic_write_json(case/'plan.json',plan);spec.update(plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=80)
    atomic_write_json(output/'namespace.json',spec);result['status']='closed_shell_ionic_native_component_pilot_prepared';return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','preparation','out'):parser.add_argument('--'+name,type=Path,required=name=='out')
    parser.add_argument('--finalize-from',type=Path)
    parser.add_argument('--native-from-review',type=Path)
    parser.add_argument('--key');parser.add_argument('--execution-profile',choices=('compact_28GiB','direct_large_GPU_reference'),default='compact_28GiB');parser.add_argument('--optimization-agreement',type=Path);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.native_from_review is not None:result=prepare_native_component(args.native_from_review.resolve(),args.out.resolve())
    elif args.finalize_from is not None:result=prepare_finalization(args.finalize_from.resolve(),args.out.resolve())
    elif args.optimization_agreement is None:result=prepare(args.archive.resolve(),args.preparation.resolve(),args.out.resolve(),args.key,args.execution_profile)
    else:result=prepare_optimization(args.archive.resolve(),args.preparation.resolve(),args.optimization_agreement.resolve(),args.out.resolve(),args.key)
    print(json.dumps(result))
