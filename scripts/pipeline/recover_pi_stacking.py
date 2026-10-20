#!/usr/bin/env python
"""Bounded recovery of the ten unconverged neutral oligomer cases.

Reuse the production queue, molecular protocol, SP evaluator and independent
gradient gate. Increase only the optimization step budget, tighten geometry-SCF
convergence, and save every geometry evaluation on NFS. Original attempts remain
immutable. No Hessian, population or electron-transfer rate is inferred.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json

SCRIPT=Path(__file__).resolve()
UPSTREAM=ROOT/'results/polymer_models/campaigns/2026-10-19_pi_stacking'
BASE=ROOT/'results/polymer_models/campaigns/2026-10-20_geometry_recovery'
spec=importlib.util.spec_from_file_location('polymer_production',ROOT/'scripts/pipeline/run_pi_stacking.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
R.BASE=BASE;R.SCRIPT=SCRIPT


from redox.workflows.geometry_checkpoints import write_geometry, checkpoint_callback


def calculate_with_checkpoints(task,plan,pilot=False):
    """Process-local hooks leave shared numerical producer files unchanged."""
    if pilot:raise RuntimeError('Recovery inherits the checked memory pilot; no repeated pilot')
    import redox.qm.dft as dft
    original_opt=dft._geomopt_fn;original_dft=dft.dft_smd;original_mf=dft._build_smd_mf
    checkpoint=Path(plan['raw_root'])/task['key']/'rest/geometry_checkpoints'
    def driver(backend):
        kernel=original_opt(backend)
        def run(method,**kwargs):
            kwargs['callback']=checkpoint_callback(checkpoint)
            result=kernel(method,**kwargs)
            write_geometry(result[1],checkpoint/'optimizer_returned.xyz',
                           'Optimizer returned coordinates; convergence='+str(bool(result[0])))
            atomic_write_json(checkpoint/'optimizer_returned.json',dict(
                optimizer_converged=bool(result[0]),solution_geometry_qualified=False,
                note='Full production SCF, topology, SP and independent gradient checks still required'))
            return result
        return run
    def energy(*args,**kwargs):
        if kwargs.get('do_opt',True):kwargs['max_opt_steps']=plan['recovery_settings']['max_opt_steps']
        return original_dft(*args,**kwargs)
    def mean_field(*args,**kwargs):
        mf=original_mf(*args,**kwargs)
        xc=args[1] if len(args)>1 else kwargs.get('xc')
        if xc=='r2scan':mf.conv_tol=plan['recovery_settings']['geometry_SCF_conv_tol']
        return mf
    try:
        dft._geomopt_fn=driver;dft.dft_smd=energy;dft._build_smd_mf=mean_field
        result=ORIGINAL_CALCULATE(task,plan,pilot=False)
        # Apply the same geometric sanity review used to admit the restart seeds.
        audit_spec=importlib.util.spec_from_file_location('polymer_geometry_audit',ROOT/'scripts/validation/audit_polymer_geometries.py')
        module=importlib.util.module_from_spec(audit_spec);audit_spec.loader.exec_module(module)
        review=module.audit(task,Path(plan['raw_root'])/task['key']/'rest/opt.xyz')
        atomic_write_json(BASE/'tasks'/task['key']/'geometry_review.json',review)
        if review['issues']:
            path=BASE/'tasks'/task['key']/'status.json';record=json.loads(path.read_text())
            record.update(status='failed_preserved',error='Geometry sanity review failed',geometry_review=review)
            atomic_write_json(path,record)
            raise RuntimeError('Geometry sanity review failed: '+str(review['issues']))
        path=BASE/'tasks'/task['key']/'status.json';record=json.loads(path.read_text())
        record['geometry_sanity_review']='passed'
        record['output_sha256'][str(BASE/'tasks'/task['key']/'geometry_review.json')]=R.sha(BASE/'tasks'/task['key']/'geometry_review.json')
        atomic_write_json(path,record)
        return result
    finally:
        dft._geomopt_fn=original_opt;dft.dft_smd=original_dft;dft._build_smd_mf=original_mf


ORIGINAL_CALCULATE=R.calculate
R.calculate=calculate_with_checkpoints


def prepare():
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    assert not (BASE/'plan.json').exists(),'Preserve existing recovery campaign'
    old=json.loads((UPSTREAM/'plan.json').read_text())
    audit_path=ROOT/'results/reviews/2026-10-20/polymer_geometry_recovery/geometry_audit.json'
    audit=json.loads(audit_path.read_text());assert not audit['flagged']
    failed={p.stem for p in (UPSTREAM/'failures').glob('*.json')}
    assert len(failed)==10 and not list((UPSTREAM/'claims').iterdir())
    pilot=json.loads((UPSTREAM/'pilot/status.json').read_text())
    assert pilot['status']=='memory_pilot_passed';R.verify(pilot['output_sha256'])
    for name in ('claims','tasks','done','failures','workers','launches','pilot','live_inputs'):(BASE/name).mkdir(parents=True,exist_ok=True)
    tasks=[]
    for task in old['tasks']:
        if task['key'] not in failed:continue
        ready_path=UPSTREAM/'live_inputs'/task['key']/'task.json'
        ready=json.loads(ready_path.read_text());R.verify(ready['input_sha256'])
        destination=BASE/'live_inputs'/task['key'];destination.mkdir()
        shutil.copyfile(ready['input_xyz'],destination/'input.xyz')
        item=dict(ready,input_xyz=str(destination/'input.xyz'),recovery_of=str(UPSTREAM/'tasks'/task['key']),
                  input_sha256={str(destination/'input.xyz'):R.sha(destination/'input.xyz'),str(ready_path):R.sha(ready_path)})
        atomic_write_json(destination/'task.json',item)
        tasks.append(task)
    shutil.copyfile(UPSTREAM/'pilot/status.json',BASE/'pilot/status.json')
    files=list((ROOT/'src/redox').rglob('*.py'))+list((ROOT/'config').glob('*.py'))+list((ROOT/'config').glob('*.json'))
    files += [SCRIPT,R.__file__,ROOT/'scripts/pipeline/collect_pi_stacking.py',ROOT/'scripts/validation/audit_polymer_geometries.py']
    inputs=[UPSTREAM/'plan.json',audit_path,BASE/'pilot/status.json']
    for task in tasks:inputs.append(BASE/'live_inputs'/task['key']/'task.json')
    plan=dict(old,created=time.time(),tasks=tasks,parent_campaign=str(UPSTREAM.relative_to(ROOT)),
              raw_root=str(ROOT/'calcs/polymer_models'/BASE.name),
              input_sha256={str(p):R.sha(p) for p in inputs},
              producer_sha256={str(p):R.sha(p) for p in sorted(set(map(Path,files)))},
              recovery_settings=dict(max_opt_steps=400,geometry_SCF_conv_tol=1e-11,
                                     original_max_opt_steps=100,acceptance_criteria='Unchanged production convergence and structural/gradient gates'),
              recovery_reason='All ten original cases exhausted 100 optimization steps; rerun checked seeds with durable checkpoints',
              original_successes_reused=17,Hessians=False)
    atomic_write_json(BASE/'plan.json',plan)
    pointer=ROOT/'results/polymer_models/recovery_campaigns.json'
    assert not pointer.exists(),'Review existing recovery overlay'
    atomic_write_json(pointer,dict(parent_campaign=str(UPSTREAM.relative_to(ROOT)),
                                  campaigns=[dict(path=str(BASE.relative_to(ROOT)),plan_sha256=R.sha(BASE/'plan.json'))]))
    print(json.dumps(dict(recovery_cases=len(tasks),preserved_original_successes=17,raw_root=plan['raw_root'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','launch','worker'])
    parser.add_argument('--role',choices=['queue_large'],default='queue_large')
    args=parser.parse_args()
    if args.action=='prepare':prepare()
    else:getattr(R,args.action)(args.role)
