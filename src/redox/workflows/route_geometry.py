"""Refine both signed branches of an independently diagnosed molecular saddle.

Shares the strict primary geometry worker/reviewer with chloride controls,
while retaining its own exact chemistry/state and 32-GB resource admission.
No energy at the old geometry or native free energy transfers to a new pose.
Top-level imports stay standard-library-only for the copied worker profile.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

PROFILE='phenothiazine_signed_neutral_geometry_v1'
SMILES='Cc1ccc(CN2c3ccccc3Sc3ccccc32)cc1'
METHOD=dict(xc='r2scan',basis='def2-svp',disp='d4',solvent='acetonitrile',
    solvent_model='primary SMD with legacy CDS',SCF_conv_tol=1e-12,SCF_conv_tol_grad=1e-8,
    SCF_max_cycle=200,grid_level=3,gradient_grid_response=True)
CONVERGENCE=dict(convergence_energy=1e-9,convergence_grms=5e-6,convergence_gmax=1e-5,
    convergence_drms=1e-4,convergence_dmax=2e-4)


def validate_plan(plan):
    expected=dict(role='route_primary_geometry_control',first_route_geometry_repair=True,
        route_geometry_profile=PROFILE,candidate_alias='phenothiazine',smiles=SMILES,charge=0,mult=1,
        atoms=39,nao=397,method=METHOD,optimizer_convergence=CONVERGENCE,
        maximum_host_RSS_GiB=32,minimum_free_GiB=32,minimum_GPU_memory_MiB=31744,
        maximum_peak_GPU_MiB=31744,maximum_wall_seconds=14400,maximum_durable_output_bytes=2*2**30,
        minimum_initial_storage_GiB=64,minimum_live_storage_GiB=32,geometry_qualified=False,
        Hessian_computed=False,native_G_qualified=False,material_admission=False)
    if (any(plan.get(k)!=v for k,v in expected.items())
            or plan.get('source_id') not in ('phenothiazine_neu_mode_minus','phenothiazine_neu_mode_plus')
            or plan.get('seed',{}).get('sign') not in (-1,1)
            or plan.get('seed',{}).get('maximum_atom_displacement_A')!=.1
            or plan.get('geometry_memory_pilot')):
        raise ValueError('Only the exact neutral phenothiazine signed geometry repairs are admitted')
    return expected


def prepare(diagnosis,output):
    import numpy as np
    from ase.io import read,write
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.workflows.molecular_native_namespace import verify
    from redox.validation.protonation import hydrogen_topology_status
    report=json.loads(diagnosis.read_text())
    if (report['status']!='route_primary_unstable_native_mode_independently_reconstructed'
            or report['native_G_qualified'] or report['charge']!=0 or report['mult']!=1
            or report['minimum_frequency_cm']>=-5
            or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('An actual independently reconstructed closed failed minimum is required')
    original=verify(Path(report['namespace']));old=json.loads((Path(original['case'])/'plan.json').read_text())
    if old['smiles']!=SMILES:raise ValueError('Exact route-admitted neutral graft required')
    if output.exists() or any(output.is_relative_to(Path(original[k])) for k in ('prior_root','archive','case','sources')):
        raise FileExistsError('Use fresh durable repair outputs outside scientific source trees')
    atoms=read(report['geometry']);mode=np.load(report['signed_mode'],allow_pickle=False)
    if mode.shape!=(39,3) or not np.isfinite(mode).all():raise ValueError('Actual Cartesian signed mode required')
    rows=[]
    for sign,label in ((-1,'minus'),(1,'plus')):
        root=output/label;case=root/'case';sources=root/'sources';case.mkdir(parents=True);sources.mkdir()
        for name in ('dft','locks','repairs'):(case/name).mkdir()
        for name in ('checkpoints.py','helper.py','kernel.cu','host_guard.py'):
            shutil.copyfile(Path(original['sources'])/name,sources/name)
        for source,name in ((ROOT/'scripts/validation/run_counterion_geometry.py','worker.py'),
                (Path(__file__),'route_geometry_profile.py'),(diagnosis,'failed_mode_diagnosis.json'),
                (Path(report['signed_mode']),'failed_mode.npy'),(Path(original['sources'])/'route_gate.json','route_gate.json')):
            shutil.copyfile(source,sources/name)
        seed=atoms.copy();seed.positions+=sign*.1*mode;write(case/'input.xyz',seed,format='xyz')
        if hydrogen_topology_status(SMILES,case/'input.xyz')!='ok':raise ValueError('Signed seed changed the intended graph')
        plan=dict(role='route_primary_geometry_control',first_route_geometry_repair=True,route_geometry_profile=PROFILE,
            candidate_alias='phenothiazine',source_id='phenothiazine_neu_mode_'+label,state='neutral_geometry_repair',
            smiles=SMILES,charge=0,mult=1,atoms=39,nao=397,method=METHOD,optimizer_convergence=CONVERGENCE,
            seed=dict(sign=sign,maximum_atom_displacement_A=.1,source_Hessian_sha256=report['Hessian_sha256'],source_diagnosis=str(diagnosis)),
            maximum_host_RSS_GiB=32,minimum_free_GiB=32,minimum_GPU_memory_MiB=31744,maximum_peak_GPU_MiB=31744,
            maximum_wall_seconds=14400,maximum_durable_output_bytes=2*2**30,minimum_initial_storage_GiB=64,minimum_live_storage_GiB=32,
            geometry_qualified=False,Hessian_computed=False,native_G_qualified=False,material_admission=False,
            source_sha256={**old['source_sha256'],**report['source_sha256'],
                **{str(p):sha256_file(p) for p in [*sources.iterdir(),case/'input.xyz'] if p.is_file()}})
        validate_plan(plan);atomic_write_json(case/'plan.json',plan)
        spec=dict(prior_root=original['prior_root'],archive=original['archive'],case=str(case),sources=str(sources),
            source_id=plan['source_id'],state=plan['state'],plan_sha256=sha256_file(case/'plan.json'),
            maximum_host_RSS_GiB=32,minimum_free_GiB=32,scientific_sources_read_only=True)
        atomic_write_json(root/'namespace.json',spec);verify(root/'namespace.json')
        rows.append(dict(key=plan['source_id'],namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json'),launcher_module='redox.workflows.route_geometry'))
    result=dict(status='two_signed_neutral_geometry_controls_prepared',states=rows,source_diagnosis=str(diagnosis),
        quantum_jobs_submitted=0,native_G_qualified=False,material_admission=False)
    atomic_write_json(output/'manifest.json',result);return result


def qualified_rows(manifest):
    from redox.core.hashing import sha256_file
    from redox.workflows.molecular_native_namespace import verify
    rows=json.loads(manifest.read_text())['states']
    if len(rows)!=2 or {r['key'] for r in rows}!={'phenothiazine_neu_mode_minus','phenothiazine_neu_mode_plus'}:
        raise ValueError('Both signed branches must remain explicit')
    for row in rows:
        if row.get('launcher_module')!='redox.workflows.route_geometry' or sha256_file(row['namespace'])!=row['namespace_sha256']:
            raise ValueError('Exact namespace and maintained launcher required')
        spec=verify(Path(row['namespace']));validate_plan(json.loads((Path(spec['case'])/'plan.json').read_text()))
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','launch','queue'))
    for name in ('diagnosis','out','spec','manifest','storage-peers'):p.add_argument('--'+name,type=Path)
    p.add_argument('--gpu',type=int);p.add_argument('--python');p.add_argument('--hours',type=float,default=12);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if a.action!='launch':
        from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.action=='prepare':
        if a.diagnosis is None or a.out is None:p.error('Preparation requires a checked diagnosis and fresh output')
        result=prepare(a.diagnosis.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],states=len(result['states']))))
    elif a.action=='launch':
        from redox.workflows.molecular_native_namespace import verify
        from redox.workflows.bounded_molecular import launch
        if a.spec is None or a.gpu is None or a.out is None:p.error('Launch requires exact spec, GPU and fresh placement')
        spec=verify(a.spec.resolve());validate_plan(json.loads((Path(spec['case'])/'plan.json').read_text()))
        result=launch(a.spec.resolve(),a.gpu,a.out.resolve(),caller_source=Path(__file__))
        print(json.dumps({k:result[k] for k in ('status','returncode')}))
    else:
        if any(v is None for v in (a.manifest,a.out,a.storage_peers,a.python)):p.error('Queue requires manifest, output, storage peers and interpreter')
        from redox.workflows.prepared_namespace_queue import run
        manifest=a.manifest.resolve();rows=qualified_rows(manifest)
        result=run(rows,a.out.resolve(),a.storage_peers.resolve(),a.python,a.hours,
            finished_status='route_primary_geometry_control_finished_independent_review_pending',
            validate=lambda row:qualified_rows(manifest),sources=[Path(__file__),manifest],resume=a.resume,label='route_geometry')
        print(json.dumps(dict(status=result['status'],native_G_qualified=False)))
