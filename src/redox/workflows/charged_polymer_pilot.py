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


def prepare(archive,preparation,output,key=None,execution_profile='compact_28GiB'):
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
    resting=[s for s in state['proposed_global_states'] if s['charge']==state['charge']]
    if len(resting)!=1:raise ValueError('One declared resting global charge/spin state required for this pilot')
    selected_state=resting[0]
    if state['charge']!=Chem.GetFormalCharge(graph) or state['charge']<=0 or selected_state['mult']!=1:
        raise ValueError('This first ionic pilot requires the resting closed-shell charged graph')
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
        source_id=spec['source_id'],state='rest',source_sha256=hashes,exact_candidate_alias=data['candidate_alias'],
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


def prepare_optimization(archive,preparation,agreement_path,output,key=None):
    agreement=json.loads(agreement_path.read_text())
    if agreement['status']!='charged_polymer_execution_profiles_compared' or not agreement['initial_state_execution_agreement_passed']:
        raise ValueError('Paired direct/compact ionic energy/force agreement required before optimization')
    for path,digest in agreement['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Ionic execution qualification changed')
    result=prepare(archive,preparation,output,key,execution_profile='compact_28GiB')
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


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','preparation','out'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--key');parser.add_argument('--execution-profile',choices=('compact_28GiB','direct_large_GPU_reference'),default='compact_28GiB');parser.add_argument('--optimization-agreement',type=Path);args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.optimization_agreement is None:result=prepare(args.archive.resolve(),args.preparation.resolve(),args.out.resolve(),args.key,args.execution_profile)
    else:result=prepare_optimization(args.archive.resolve(),args.preparation.resolve(),args.optimization_agreement.resolve(),args.out.resolve(),args.key)
    print(json.dumps(result))
