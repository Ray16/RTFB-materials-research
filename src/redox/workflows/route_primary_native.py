"""First native components for the route-admitted phenothiazine baseline.

Reuses the primary molecular Hessian resource envelope and quantum/partition
engine. Each state is a prospective pilot and needs its own sulfur derivatives,
wavefunction, minimum, nuclear and closed-execution review. No potential is
published by this module. Top-level imports stay standard-library-only because
the profile is copied into the isolated scientific worker namespace.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import shutil

PROFILE='route_phenothiazine_first_native_v1'
SMILES='Cc1ccc(CN2c3ccccc3Sc3ccccc32)cc1'


def validate_plan(plan):
    states={'phenothiazine_neu':(0,1,'neu'),'phenothiazine_ox':(1,2,'ox')}
    expected=dict(role='primary_molecular_Hessian_control',primary_molecular_Hessian_control=True,
        first_route_native_component=True,route_native_profile=PROFILE,candidate_alias='phenothiazine',
        smiles=SMILES,atoms=39,native_nao=397,native_memory_cap_GiB=6,
        maximum_host_RSS_GiB=40,minimum_free_GiB=32,minimum_GPU_memory_MiB=31744,
        maximum_peak_GPU_MiB=31744,maximum_wall_seconds=14400,maximum_durable_output_bytes=2**30,
        minimum_initial_storage_GiB=33,minimum_live_storage_GiB=32,
        extra_columns=[[15,0],[15,1],[15,2]],geometry_repair_permitted=False,
        native_G_qualified=False,redox_potential_qualified=False,material_admission=False)
    if (any(plan.get(k)!=v for k,v in expected.items()) or plan.get('source_id') not in states
            or (plan.get('charge'),plan.get('mult'),plan.get('native_state'))!=states.get(plan.get('source_id'))):
        raise ValueError('Only the exact admitted phenothiazine first-native pilot is supported')
    return expected


def prepare(archive,package,energy_review,output,contract):
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.core.protocol import ACTIVE_SP
    from redox.core.runtime_paths import gpu4pyscf_source_root
    from redox.workflows.molecular_native_namespace import load,verify
    from ase.io import read
    from pyscf import gto
    report=json.loads(energy_review.read_text());source_manifest=json.loads((package/'manifest.json').read_text())
    if (report['status']!='benchmark_screening_energy_control_independently_checked'
            or not report['actual_checkpoint_frames_verified'] or not report['declared_agreement_passed']
            or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('Actual independently checked corrected electronic parent required')
    gate=json.loads((package/'route_gate.json').read_text())
    if (not gate['proceed_to_full_redox'] or not gate['exact_route_binding_verified']
            or any(sha256_file(package/p)!=h for p,h in source_manifest['files_sha256'].items())):
        raise ValueError('The exact frozen prospective route and source state must remain admitted')
    row=next(r for r in source_manifest['states'] if r['id']==report['key'])
    if report['source_state']['smiles']!=SMILES or (row['charge'],row['mult'])!=(report['source_state']['charge'],report['source_state']['mult']):
        raise ValueError('Exact baseline charge/spin/graph differs from the electronic parent')
    electronic=verify(Path(report['namespace']));parent=Path(electronic['case'])
    geometry=Path(electronic['sources'])/'input.xyz';atoms=read(geometry)
    if sha256_file(geometry)!=source_manifest['files_sha256'][row['files']['geometry']]:
        raise ValueError('The native geometry must equal the checked electronic geometry exactly')
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='def2-svp',charge=row['charge'],spin=row['mult']-1,verbose=0)
    if (len(atoms),mol.nao,atoms.get_chemical_symbols()[15])!=(39,397,'S'):raise ValueError('Exact sulfur dimension/order pilot required')
    if output.exists() or any(output.is_relative_to(p) for p in (ROOT,archive)):
        raise FileExistsError('Fresh durable native outputs outside scientific source trees required')
    case=output/'case';sources=output/'sources';case.mkdir(parents=True);sources.mkdir()
    for name in ('state','library','dft','locks','repairs'):(case/name).mkdir()
    for p,name in [(ROOT/'scripts/validation/run_route_primary_native.py','worker.py'),
            (ROOT/'src/redox/core/hashing.py','hashing.py'),(ROOT/'src/redox/qm/scf_checkpoints.py','checkpoints.py'),
            (ROOT/'src/redox/qm/screening_barrier.py','helper.py'),(ROOT/'tests/fixtures/screen_index_v1_8_1.cu','kernel.cu'),
            (ROOT/'scripts/fleet/bounded_peak_host_python_worker.py','host_guard.py'),(Path(__file__),'admission.py'),
            (energy_review,'electronic_review.json'),(parent/'corrected/result.json','corrected_electronic.json'),
            (geometry,'reference_geometry.xyz'),(package/'route_gate.json','route_gate.json'),(contract,'calculation_contract.json')]:
        shutil.copyfile(p,sources/name)
    shutil.copyfile(geometry,case/'state/opt.xyz')
    raw_source=package/row['files']['result'];shutil.copyfile(raw_source,case/'state/result.json');raw=json.loads(raw_source.read_text())
    with (case/'library/manifest.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['id','state','smiles','charge','mult_hint']);writer.writeheader()
        writer.writerow(dict(id=raw['id'],state=raw['state'],smiles=SMILES,charge=row['charge'],mult_hint=row['mult']))
    native=load(archive/'src/redox/qm/native_thermal.py','route_primary_native_protocol')
    thermal=load(archive/'src/redox/qm/qrrho.py','route_primary_thermal_protocol')
    installed=gpu4pyscf_source_root()
    files=[*sources.glob('*'),case/'state/opt.xyz',case/'state/result.json',case/'library/manifest.csv',
        *(archive/'src/redox/qm').glob('*.py'),*(archive/'src/redox/core').glob('*.py'),
        *(archive/'config').glob('*.py'),*(archive/'config').glob('*.json'),archive/'src/redox/workflows/resources.py',
        archive/'scripts/pipeline/run_pi_stacking.py',archive/'scripts/validation/audit/native_partition.py',
        archive/'scripts/fleet/exclusive_python_worker.py',archive/'scripts/fleet/tensor_runtime_worker.py',
        archive/'environments/cutensor22/runtime.json',installed/'dft/gen_grid.py',installed/'dft/numint.py',installed/'lib/libgdft.so']
    plan=dict(role='primary_molecular_Hessian_control',primary_molecular_Hessian_control=True,
        first_route_native_component=True,route_native_profile=PROFILE,candidate_alias='phenothiazine',
        source_id=row['id'],state='rest',native_id=raw['id'],native_state=raw['state'],smiles=SMILES,
        charge=row['charge'],mult=row['mult'],atoms=39,native_nao=397,geometry_sha256=sha256_file(geometry),
        native_protocol=native.NATIVE_THERMAL,thermal_protocol=thermal.PROTOCOL,electronic_protocol=ACTIVE_SP,
        native_memory_cap_GiB=6,maximum_host_RSS_GiB=40,minimum_free_GiB=32,minimum_GPU_memory_MiB=31744,
        maximum_peak_GPU_MiB=31744,maximum_wall_seconds=14400,maximum_durable_output_bytes=2**30,
        minimum_initial_storage_GiB=33,minimum_live_storage_GiB=32,extra_columns=[[15,0],[15,1],[15,2]],
        geometry_repair_permitted=False,prospective_resource_pilot=True,native_G_qualified=False,
        redox_potential_qualified=False,material_admission=False,
        source_sha256={**report['source_sha256'],**{str(p):sha256_file(p) for p in files}})
    validate_plan(plan);atomic_write_json(case/'plan.json',plan)
    spec=dict(prior_root=str(ROOT),archive=str(archive),case=str(case),sources=str(sources),source_id=plan['source_id'],state='rest',
        plan_sha256=sha256_file(case/'plan.json'),maximum_host_RSS_GiB=40,minimum_free_GiB=32,
        scientific_sources_read_only=True,library_overlay=str(case/'library'))
    atomic_write_json(output/'namespace.json',spec);verify(output/'namespace.json')
    return dict(key=plan['source_id'],namespace=str(output/'namespace.json'),namespace_sha256=sha256_file(output/'namespace.json'))


def queued_rows(manifest_path):
    from redox.workflows.molecular_native_namespace import verify
    data=json.loads(manifest_path.read_text());rows=data['states']
    if len(rows)!=2 or {r['key'] for r in rows}!={'phenothiazine_neu','phenothiazine_ox'}:
        raise ValueError('Both exact declared native states must remain present')
    for row in rows:
        spec=verify(Path(row['namespace']));validate_plan(json.loads((Path(spec['case'])/'plan.json').read_text()))
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','queue'))
    for name in ('archive','package','energy-review','out','contract','manifest','storage-peers'):p.add_argument('--'+name,type=Path)
    p.add_argument('--python');p.add_argument('--hours',type=float,default=12);p.add_argument('--resume',action='store_true');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':
        for k in ('archive','package','energy_review','out','contract'):
            if getattr(a,k) is None:p.error('Preparation requires archive, package, checked energy, contract and fresh output')
        print(json.dumps(prepare(a.archive.resolve(),a.package.resolve(),a.energy_review.resolve(),a.out.resolve(),a.contract.resolve())))
    else:
        if any(v is None for v in (a.manifest,a.out,a.storage_peers,a.python)):p.error('Queue requires manifest, storage roster, interpreter and fresh output')
        from redox.workflows.prepared_namespace_queue import run
        manifest=a.manifest.resolve();rows=queued_rows(manifest)
        result=run(rows,a.out.resolve(),a.storage_peers.resolve(),a.python,a.hours,
            finished_status='route_primary_native_component_finished_review_pending',validate=lambda row:queued_rows(manifest),
            sources=[Path(__file__),manifest],resume=a.resume,label='route_primary_native')
        print(json.dumps(dict(status=result['status'],native_G_qualified=False)))
