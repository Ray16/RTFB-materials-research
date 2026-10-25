#!/usr/bin/env python
"""Independently gated smooth-CDS geometry pilot with a distinct record schema.

Original structural and energy/force gates are retained. The original failed
geometry is only a warm start. No primary result.json or electronic/native
certificate is written for the alternate solvent model.
"""
import argparse
from contextlib import contextmanager
import hashlib
import inspect
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

GEOMETRY=ROOT/'scripts/pipeline/optimize_polymer_geometry.py'
MODEL=ROOT/'scripts/pipeline/smooth_cds_control.py'
CHECKER=ROOT/'scripts/validation/review_smooth_model_panel.py'


def alternate_run_source(body):
    if body.count("sd/'result.json'")!=2:
        raise ValueError('Geometry writer template changed; do not risk a primary-schema write')
    body=body.replace("sd/'result.json'","sd/'experimental_geometry_result.json'")
    old='forcecheck=forces.final_energy_check(sd/'
    if body.count(old)!=1:raise ValueError('Final force-check template changed')
    body=body.replace(old,'forcecheck=_expanded_final_check(forces,sd/',1)
    old='convergence_gmax=1e-4,convergence_grms=6e-5,convergence_energy=1e-7,\n                        convergence_dmax=.002,convergence_drms=.001'
    if body.count(old)!=1:raise ValueError('Geometry convergence template changed')
    return body.replace(old,'convergence_gmax=1e-5,convergence_grms=6e-6,convergence_energy=1e-8,\n                        convergence_dmax=.0005,convergence_drms=.0003',1)


def expanded_force_source(body,axes):
    old='indices = np.argsort(abs(force).ravel())[-2:]'
    if body.count(old)!=1:raise ValueError('Final force selector template changed')
    return body.replace(old,old+'\n    indices = np.unique(np.r_[indices, '+repr([3*a+b for a,b in axes])+']).astype(int)',1)


def prepare(source,controls,out,axes):
    plan=api.read_plan(source)
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    status_path=source/'geometry/status.json';state=json.loads(status_path.read_text())
    owner=api.load_module('smooth_seed_stopped_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if (not state.get('finished') or owner.owner_alive(state['owner'])
            or 'force fails energy differences' not in state.get('error','')):
        raise ValueError('A completed, preserved force-check failure is required for this pilot')
    if len(plan['tasks'])!=1 or (plan['tasks'][0]['charge'],plan['tasks'][0]['mult'])!=(0,1):
        raise ValueError('Pilot is scoped to one neutral singlet graph')
    task=plan['tasks'][0]
    if any(not 0<=a<task['atoms'] or not 0<=b<3 for a,b in axes):raise ValueError('Invalid critical axis')
    model=api.load_module('smooth_geometry_declared_model',MODEL)
    checked=api.load_module('smooth_geometry_quantum_controls',CHECKER).checked_controls(controls)
    sd=Path(state['state_dir']);seed=sd/'opt.xyz'
    raw_result=json.loads((sd/'result.json').read_text())
    if not raw_result['geometry_optimization_converged'] or not raw_result['optimization_SCF_converged']:
        raise ValueError('Source must have an optimizer-returned, SCF-converged geometry')
    out.mkdir(parents=True,exist_ok=False)
    shutil.copytree(source/'live_inputs',out/'live_inputs')
    ready=out/'live_inputs'/task['key'];xyz=ready/'input.xyz';shutil.copyfile(seed,xyz)
    audit=api.load_module('smooth_initial_graph_review',ROOT/'scripts/validation/audit_polymer_geometries.py')
    review=audit.audit(task,xyz)
    if review['issues']:raise ValueError('Source graph/stereochemistry is invalid')
    atomic_write_json(ready/'initial_graph_review.json',review)
    row=json.loads((ready/'task.json').read_text());row['input_xyz']=str(xyz);row['input_sha256'][str(xyz)]=api.digest(xyz)
    atomic_write_json(ready/'task.json',row)
    for stage in ('energy','thermal'):
        for name in ('tasks','done','claims','failures','workers','pilot','launches'):
            (out/stage/name).mkdir(parents=True)
    inputs=dict(plan['input_sha256'])
    for p in (source/'plan.json',status_path,sd/'result.json',seed,xyz,ready/'task.json',ready/'initial_graph_review.json'):
        inputs[str(p)]=api.digest(p)
    runtime=ROOT/'environments/cutensor22/runtime.json';inputs[str(runtime)]=api.digest(runtime)
    producers={**plan['producer_sha256'],**checked}
    for p in (Path(__file__).resolve(),GEOMETRY,MODEL,CHECKER,ROOT/'scripts/pipeline/smd_gc_workspace.py',
              ROOT/'scripts/fleet/tensor_runtime_worker.py'):
        producers[str(p)]=api.digest(p)
    updated=dict(plan,input_sha256=inputs,producer_sha256=producers,
        raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),execution_driver=str(Path(__file__).resolve()),
        runtime_profile=str(runtime),maximum_host_RSS_GiB=32,
        quantum_settings=dict(plan['quantum_settings'],geometry_SCF_conv_tol=1e-12),
        geometry_resource_policy=dict(minimum_GPU_MiB=32768,maximum_peak_GPU_MiB=31744,minimum_host_RAM_GiB=64),
        experimental_model=model.MODEL,experimental_method_id=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest(),
        full_quantum_controls=str(controls),extra_force_axes=axes,primary_candidate_admission=False,
        purpose='New solvent-model geometry pilot; original failed coordinates only a warm start, all physical gates repeated',created=time.time())
    atomic_write_json(out/'plan.json',updated)
    print(json.dumps(dict(status='experimental_smooth_geometry_pilot_prepared',primary_admission=False)))


def run(campaign):
    import os
    from redox.workflows.resources import host_capacity
    plan=api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    guard=json.loads((campaign/'host_guard/status.json').read_text())
    if guard['pid']!=os.getpid() or guard['maximum_RSS_GiB']!=32 or guard['minimum_free_GiB']!=32:
        raise ValueError('Declared same-PID host-memory guard required')
    model=api.load_module('smooth_geometry_profile_runtime',MODEL)
    checked=api.load_module('smooth_geometry_full_quantum_review',CHECKER).checked_controls(Path(plan['full_quantum_controls']))
    if plan['experimental_model']!=model.MODEL:raise ValueError('Experimental model declaration changed')
    geometry=api.load_module('retained_geometry_physical_gates',GEOMETRY)
    raw=Path(plan['raw_root'])
    geometry.host_capacity=lambda threads,**kwargs:host_capacity(threads,minimum_ram_gib=64)
    engine_body=inspect.getsource(geometry.geometry_engine)
    old=",forces.adaptive_forces(directory/'CDS')"
    if engine_body.count(old)!=1:raise ValueError('Geometry-force context template changed')
    engine_body=engine_body.replace(old,'',1)
    namespace=dict(vars(geometry));local={}
    exec(compile(engine_body,str(Path(__file__).resolve()),'exec'),namespace,local)
    geometry.geometry_engine=local['geometry_engine']
    original_writer=geometry.atomic_write_json
    def write(path,value):
        if path.name=='experimental_geometry_result.json' or path==campaign/'geometry/status.json':
            value=dict(value,experimental_model=plan['experimental_model'],experimental_method_id=plan['experimental_method_id'],
                primary_candidate_admission=False)
            if value.get('status')=='stationary_geometry_checked_active_parent_and_native_minimum_pending':
                value['status']='experimental_geometry_base_checks_passed_model_review_pending'
        original_writer(path,value)
    geometry.atomic_write_json=write
    def final_check(module,xyz,out):
        body=expanded_force_source(inspect.getsource(module.final_energy_check),plan['extra_force_axes'])
        local={};scope=dict(vars(module))
        exec(compile(body,str(Path(__file__).resolve()),'exec'),scope,local)
        return local['final_energy_check'](xyz,out)
    namespace=dict(vars(geometry),_expanded_final_check=final_check);local={}
    exec(compile(alternate_run_source(inspect.getsource(geometry.run)),str(Path(__file__).resolve()),'exec'),namespace,local)
    gc=api.load_module('smooth_geometry_GC_guard',ROOT/'scripts/pipeline/smd_gc_workspace.py')
    metrics=dict(collections=0,unreachable_objects_collected=0)
    with gc.collect_before_build(metrics),model.smooth_cds(campaign/'model',ROOT/'results/validation/smooth_CDS_2026-10-25/pilot'):
        local['run'](campaign)
    path=campaign/'geometry/status.json';record=json.loads(path.read_text());proof=campaign/'model/model.json'
    if record['status']!='experimental_geometry_base_checks_passed_model_review_pending':
        raise ValueError('Experimental model or base geometry checks did not finish normally')
    record.update(status='experimental_smooth_geometry_checked_energy_and_native_QC_pending',
        model_evidence=str(proof),full_quantum_control_evidence=checked,finished=time.time())
    record['output_sha256'][str(proof)]=api.digest(proof)
    atomic_write_json(path,record)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','run'))
    parser.add_argument('--source',type=Path);parser.add_argument('--controls',type=Path)
    parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--critical-axis',nargs=2,type=int,action='append',default=[])
    args=parser.parse_args()
    if args.action=='prepare':
        if not args.source or not args.controls:parser.error('prepare requires source and controls')
        prepare(args.source.resolve(),args.controls.resolve(),args.campaign.resolve(),args.critical_axis)
    else:run(args.campaign.resolve())


if __name__=='__main__':main()
