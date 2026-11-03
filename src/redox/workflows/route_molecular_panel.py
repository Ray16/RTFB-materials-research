"""Start an exact route-admitted molecule from preserved source geometries.

The first stage reuses the established installed/corrected electronic worker
and its independent review. Own native thermochemistry and matched references
remain required before publishing a potential or reaction free energy.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import source_state, load, verify
from redox.build.contracts import state_declaration


def route_candidate(registry, routes, assessment, alias, source_id):
    records=json.loads(Path(registry).read_text())['candidates']
    matches=[r for r in records if alias in r['aliases'] and source_id in r['aliases']]
    if len(matches)!=1:raise ValueError('One exact graph must bind the candidate and source identity')
    candidate=matches[0];declaration=state_declaration(candidate,alias)
    gate=load(ROOT/'scripts/pipeline/screen_grafting_routes.py','route_molecular_exact_gate')
    decision=gate.screen(dict(id=alias,frag=declaration['fragment'],product_smiles=candidate['product_smiles']),
        json.loads(Path(routes).read_text())['candidates'],json.loads(Path(assessment).read_text()))
    if not decision['proceed_to_full_redox'] or not decision.get('exact_route_binding_verified'):
        raise ValueError('Exact precursor/product route must pass before a full molecular campaign')
    return candidate,declaration,decision


def prepare(archive,registry,routes,assessment,alias,source_id,output):
    from redox.workflows.screening_component_panel import prepare as prepare_energy
    archive=Path(archive).resolve();output=Path(output).resolve()
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable output outside both scientific source trees required')
    candidate,declaration,decision=route_candidate(registry,routes,assessment,alias,source_id)
    states=declaration['states']
    if len(states)<2 or len({s[0] for s in states})!=len(states):
        raise ValueError('A complete distinct declared redox-state set is required')
    checked=[];bindings={str(Path(p).resolve()):sha256_file(p) for p in (registry,routes,assessment)}
    for state,q,m,_ in states:
        row,raw,energy=source_state(archive,source_id,state,candidate['product_smiles'])
        if (raw['charge'],raw['mult'])!=(q,m):raise ValueError('Source and declared redox assignments differ')
        directory=energy.parent.parent
        checked.append((state,q,m,row,raw,directory,energy))
        bindings.update({str(p):sha256_file(p) for p in (directory/'opt.xyz',directory/'result.json',energy)})
    manifest_path=archive/'library/manifest.csv';bindings[str(manifest_path)]=sha256_file(manifest_path)
    package=output/'source_package';package.mkdir(parents=True)
    for value,name in ((candidate,'candidate.json'),(decision,'route_gate.json')):atomic_write_json(package/name,value)
    shutil.copyfile(Path(__file__),package/'preparer.py')
    rows=[]
    for state,q,m,row,raw,directory,energy in checked:
        case=package/state;case.mkdir()
        for source,name in ((directory/'opt.xyz','geometry.xyz'),(directory/'result.json','result.json'),(energy,'energy.json')):
            shutil.copyfile(source,case/name)
        rows.append(dict(id=source_id+'_'+state,source_id=source_id,source_state=state,
            smiles=candidate['product_smiles'],charge=q,mult=m,
            files={key:str((case/name).relative_to(package)) for key,name in
                (('geometry','geometry.xyz'),('result','result.json'),('energy','energy.json'))},
            native_G_qualified=False,checked_G_by_cutoff_eV=None))
    if any(sha256_file(p)!=h for p,h in bindings.items()):raise ValueError('Candidate or source changed during capture')
    manifest=dict(scope_id='route_admitted_'+alias+'_primary_electronic',source_alias_prefix=alias,
        candidate_alias=alias,source_id=source_id,product_smiles=candidate['product_smiles'],product_inchikey=candidate['inchikey'],
        states=rows,files_sha256={str(p.relative_to(package)):sha256_file(p) for p in package.rglob('*') if p.is_file()},
        source_sha256=bindings,route_gate_snapshot='route_gate.json',source_geometries_reused=True,
        geometry_scope='One preserved previously optimized geometry per state; fresh stationarity/minimum and conformer coverage are separate',
        native_G_qualified=False,redox_potential_qualified=False,material_admission=False)
    atomic_write_json(package/'manifest.json',manifest)
    result=prepare_energy(archive,package,output/'electronic')
    result.update(candidate_alias=alias,source_id=source_id,route_gate_snapshot=str(package/'route_gate.json'),
        required_next_stages=['independent actual electronic/frame/resource review','own native stationarity/Hessian/modes/partition',
            'matched corrected Fc reference and context-checked redox sum'],material_admission=False)
    atomic_write_json(output/'prepared.json',result)
    return result


def qualified_rows(manifest_path):
    data=json.loads(manifest_path.read_text());package=Path(data['source_package'])
    source=json.loads((package/'manifest.json').read_text());gate=json.loads((package/'route_gate.json').read_text())
    if (sha256_file(package/'manifest.json')!=data['source_package_manifest_sha256']
            or not gate['proceed_to_full_redox'] or not gate['exact_route_binding_verified']
            or len(data['states'])!=len(source['states']) or len(data['states'])!=data['target_states']):
        raise ValueError('The exact full state set and frozen prospective route admission are required')
    if any(sha256_file(package/p)!=h for p,h in source['files_sha256'].items()):
        raise ValueError('Frozen route, molecule or source input changed')
    states={r['id']:r for r in source['states']}
    if len(states)!=len(source['states']) or {r['key'] for r in data['states']}!=set(states):
        raise ValueError('All declared exact molecular states must remain present exactly once')
    for row in data['states']:
        spec=verify(Path(row['namespace']));plan=json.loads((Path(spec['case'])/'plan.json').read_text());state=states[row['key']]
        if (plan['role']!='screening_energy_control' or not plan['energy_control']
                or plan['source_state']['smiles']!=source['product_smiles']
                or (plan['source_state']['charge'],plan['source_state']['mult'])!=(state['charge'],state['mult'])
                or sha256_file(Path(spec['sources'])/'input.xyz')!=source['files_sha256'][state['files']['geometry']]):
            raise ValueError('Prepared electronic state differs from the exact admitted molecule')
    return data['states'],package


def dispatch(manifest_path,interpreter,storage_peers=()):
    from redox.workflows.charged_native_followup import dispatch_rows
    rows,package=qualified_rows(manifest_path)
    return dispatch_rows(rows,manifest_path.parent,interpreter,
        [manifest_path,package/'manifest.json',package/'route_gate.json',Path(__file__)],storage_peers=storage_peers)


def queue(manifest_path,output,storage_peers,interpreter,hours,resume=False):
    from redox.workflows.prepared_namespace_queue import run
    rows,package=qualified_rows(manifest_path)
    return run(rows,output,storage_peers,interpreter,hours,
        finished_status='screening_energy_control_finished_review_pending',
        validate=lambda row:qualified_rows(manifest_path),
        sources=[Path(__file__),manifest_path,package/'manifest.json',package/'route_gate.json'],
        resume=resume,label='route_molecular_electronic')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','dispatch','queue'))
    for name in ('archive','registry','routes','assessment','out','manifest','storage-peers'):p.add_argument('--'+name,type=Path)
    for name in ('candidate-alias','source-id','python'):p.add_argument('--'+name)
    p.add_argument('--hours',type=float,default=12);p.add_argument('--resume',action='store_true')
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':
        if any(v is None for v in (a.archive,a.registry,a.routes,a.assessment,a.candidate_alias,a.source_id,a.out)):
            p.error('Preparation requires archive, exact registry/route/assessment, candidate/source IDs and fresh output')
        result=prepare(a.archive,a.registry,a.routes,a.assessment,a.candidate_alias,a.source_id,a.out)
        print(json.dumps({k:result[k] for k in ('status','candidate_alias','target_states','quantum_jobs_submitted')}))
    elif a.action=='dispatch':
        if a.manifest is None or not a.python:p.error('Dispatch requires prepared manifest and node interpreter')
        peers=[] if a.storage_peers is None else json.loads(a.storage_peers.read_text())
        result=dispatch(a.manifest.resolve(),a.python,peers)
        print(json.dumps({k:result[k] for k in ('status','launches','refusals','unsubmitted')}))
    else:
        if any(v is None for v in (a.manifest,a.out,a.storage_peers,a.python)):
            p.error('Queue requires prepared manifest, storage roster, node interpreter and fresh controller output')
        result=queue(a.manifest.resolve(),a.out.resolve(),a.storage_peers.resolve(),a.python,a.hours,a.resume)
        print(json.dumps(dict(status=result['status'],scientific_retries_performed=0)))
