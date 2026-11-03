"""Separate T4 fixed-gradient validation; full molecular/Hessian limits stay intact."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil


PROFILE='primary_chloride_T4_fixed_gradient_pilot_v1'
ROLE='primary_counterion_fixed_gradient_memory_pilot'
METHOD=dict(xc='r2scan',basis='def2-svp',disp='d4',solvent='acetonitrile',
    solvent_model='primary SMD with legacy CDS',SCF_conv_tol=1e-12,
    SCF_conv_tol_grad=1e-8,SCF_max_cycle=200,grid_level=3,gradient_grid_response=True)


def validate_plan(plan):
    expected=dict(role=ROLE,geometry_resource_profile=PROFILE,fixed_gradient_memory_pilot=True,
        atoms=42,nao=412,charge=1,mult=1,maximum_host_RSS_GiB=8,minimum_free_GiB=32,
        minimum_GPU_memory_MiB=15360,maximum_peak_GPU_MiB=12288,
        maximum_wall_seconds=3600,maximum_durable_output_bytes=256*2**20,
        Hessian_computed=False,geometry_optimization_performed=False,native_G_qualified=False)
    if any(plan.get(k)!=v for k,v in expected.items()) or plan.get('method')!=METHOD:
        raise ValueError('Only the exact prospective 42-atom fixed-gradient T4 control is admitted')
    return expected


def device_for_plan(plan,gpu):
    validate_plan(plan)
    from redox.workflows.resources import inspect_gpu
    device=inspect_gpu(gpu)
    if (device['memory_MiB']<15360 or float(device['compute_capability'])!=7.5
            or 'T4' not in device['name'].split()):
        raise RuntimeError('This independent fixed-gradient resource pilot specifically tests a free T4')
    return device


def prepare(geometry_review,output):
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.workflows.molecular_native_namespace import verify,load
    review=json.loads(geometry_review.read_text())
    if (review['status']!='primary_counterion_geometry_independently_checked'
            or (review['total_charge'],review['multiplicity'],review['formal_core_charge'])!=(1,1,2)
            or review['method']!=METHOD or review['nao']!=412
            or review['native_G_qualified'] or review['material_admission']
            or any(sha256_file(p)!=v for p,v in review['source_sha256'].items())):
        raise ValueError('Exact independently checked primary chloride geometry is required')
    source=verify(Path(review['namespace']));origin=Path(source['case']);old=json.loads((origin/'plan.json').read_text())
    if output.exists() or any(output.is_relative_to(Path(source[k])) for k in ('prior_root','archive')):
        raise FileExistsError('Use a fresh durable control outside scientific source trees')
    case=output/'case';sources=output/'sources';case.mkdir(parents=True);sources.mkdir()
    for name in ('dft','locks','repairs'):(case/name).mkdir()
    for path in Path(source['sources']).iterdir():
        if path.is_file():shutil.copyfile(path,sources/path.name)
    for path,name in [(ROOT/'scripts/validation/run_counterion_gradient_pilot.py','worker.py'),
        (Path(__file__),'resource_profile.py'),(geometry_review,'reference_geometry_review.json'),
        (origin/'final_evaluation.json','reference_final_evaluation.json'),
        (origin/'optimizer_return.xyz','input.xyz')]:shutil.copyfile(path,sources/name)
    if sha256_file(sources/'input.xyz')!=review['geometry_sha256']:raise ValueError('Exact source geometry required')
    plan=dict(old);plan['source_sha256']=dict(old['source_sha256'])
    plan.update(role=ROLE,geometry_resource_profile=PROFILE,fixed_gradient_memory_pilot=True,
        maximum_host_RSS_GiB=8,minimum_free_GiB=32,minimum_GPU_memory_MiB=15360,
        maximum_peak_GPU_MiB=12288,maximum_wall_seconds=3600,
        maximum_durable_output_bytes=256*2**20,Hessian_computed=False,
        geometry_optimization_performed=False,native_G_qualified=False,
        primary_candidate_admission=False,geometry_followups_admitted=False,
        purpose='Independent T4 SCF/gradient memory and architecture comparison only; no Hessian or native G')
    plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.iterdir() if p.is_file()})
    validate_plan(plan);atomic_write_json(case/'plan.json',plan)
    spec=dict(source,case=str(case),sources=str(sources),maximum_host_RSS_GiB=8,
              plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json',spec);verify(output/'namespace.json')
    return dict(status='T4_fixed_gradient_pilot_prepared_no_followup_admission',namespace=str(output/'namespace.json'),quantum_jobs_submitted=0)


def review(namespace):
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    from redox.core.io import atomic_write_json
    from redox.workflows.molecular_native_namespace import verify,load
    import numpy as np
    from ase.io import read as read_atoms
    from redox.validation.execution_review import review_closure,freeze_source
    from redox.validation.primary_checkpoints import checkpoint_series
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(path):used.append(Path(path));return json.loads(Path(path).read_text())
    plan=read(case/'plan.json');validate_plan(plan)
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','T4_gradient_closed_owners')
    stage=review_closure(case,read,owners.alive,expected_status='counterion_fixed_gradient_finished_review_pending',
        host_GiB=8,GPU_MiB=12288,worker=sources/'worker.py')
    result=read(case/'evaluation.json');reference=read(sources/'reference_final_evaluation.json')
    geometry=sources/'input.xyz';used.append(geometry);atoms=read_atoms(geometry)
    series=checkpoint_series(case/'SCF_evaluations',geometry,
        dict(charge=1,mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,
             conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),read,used,
             expected_count=1,fixed_geometry=True)
    energy=abs(result['E_Ha']-reference['E_Ha']);gradient=np.asarray(result['gradient_Ha_Bohr'])
    force=float(abs(gradient-np.asarray(reference['gradient_Ha_Bohr'])).max())
    gradient_paths=sorted((case/'SCF_evaluations').glob('*/*/gradient.npy'))
    if len(gradient_paths)!=1:raise ValueError('One actual checkpointed gradient required')
    used.append(gradient_paths[0]);actual_gradient=np.load(gradient_paths[0],allow_pickle=False)
    reference_review=read(sources/'reference_geometry_review.json')
    if any(sha256_file(p)!=h for p,h in reference_review['source_sha256'].items()):
        raise ValueError('The independent reference geometry backing changed')
    origin=verify(Path(reference_review['namespace']));old_case=Path(origin['case'])
    old_gradients=sorted((old_case/'SCF_evaluations').glob('*/*/gradient.npy'))
    if not old_gradients:raise ValueError('Actual reference gradient evidence required')
    used.append(old_gradients[-1]);reference_gradient=np.load(old_gradients[-1],allow_pickle=False)
    runtime=read(case/'tensor_runtime/status.json')
    if (not result['SCF_converged'] or gradient.shape!=(42,3) or not np.isfinite(gradient).all()
            or actual_gradient.shape!=(42,3) or not np.isfinite(actual_gradient).all()
            or np.max(abs(actual_gradient-gradient))>1e-12
            or reference_gradient.shape!=(42,3) or not np.isfinite(reference_gradient).all()
            or np.max(abs(reference_gradient-np.asarray(reference['gradient_Ha_Bohr'])))>1e-12
            or abs(series['energies_Ha'][0]-result['E_Ha'])>1e-10
            or runtime['compute_capability']!=[7,5] or 'T4' not in runtime['device_name'].split()
            or not math.isfinite(energy) or energy>1e-8 or not math.isfinite(force) or force>1e-6
            or stage['observed_peak_GPU_MiB']*1.25>15360):
        raise ValueError('T4 resource headroom and exact fixed-state energy/force agreement required')
    used.extend([freeze_source(case,Path(__file__)),freeze_source(case,ROOT/'src/redox/validation/primary_checkpoints.py')])
    return dict(status='T4_primary_chloride_fixed_gradient_independently_checked',namespace=str(namespace),
        geometry_sha256=sha256_file(geometry),atoms=len(atoms),charge=1,mult=1,method=METHOD,
        energy_difference_Ha=energy,maximum_gradient_difference_Ha_Bohr=force,
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],actual_SCF_frames_checked=series['actual_SCF_frames_checked'],
        fixed_gradient_resource_control_passed=True,full_geometry_resource_class_qualified=False,
        Hessian_class_qualified=False,native_G_qualified=False,material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used})


if __name__=='__main__':
    from redox.core.io import atomic_write_json
    from redox.workflows.molecular_native_namespace import verify
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','launch','review'))
    for name in ('geometry-review','spec','out'):parser.add_argument('--'+name,type=Path)
    parser.add_argument('--gpu',type=int);args=parser.parse_args()
    if args.action!='launch':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.action=='prepare':result=prepare(args.geometry_review.resolve(),args.out.resolve())
    elif args.action=='review':
        if args.out.exists():raise FileExistsError('Preserve existing resource reviews')
        result=review(args.spec.resolve());atomic_write_json(args.out,result)
    else:
        spec=verify(args.spec.resolve());validate_plan(json.loads((Path(spec['case'])/'plan.json').read_text()))
        from redox.workflows.bounded_molecular import launch
        result=launch(args.spec.resolve(),args.gpu,args.out.resolve(),caller_source=Path(__file__))
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','argv')}))
