#!/usr/bin/env python
"""Refine an independently reviewed counterpoise checkpoint in a new attempt.

The frozen geometry, common integration grid, auxiliary basis and scientific
acceptance threshold are preserved. Only the SCF stopping targets are tightened.
"""
import argparse
import inspect
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

ORIGINAL=ROOT/'scripts/pipeline/counterpoise_cpu_remaining.py'
TARGETS=dict(SCF_conv_tol_Ha=1e-13,SCF_conv_tol_grad_au=1e-9,
             independent_acceptance_threshold_au=1e-7,scientific_model_changed=False)


def prior_checkpoint(prior,checkpoint_review,stationarity_review):
    raw=json.loads((prior/'status.json').read_text());launch=json.loads((prior/'launch.json').read_text())
    review=json.loads(checkpoint_review.read_text());stationarity=json.loads(stationarity_review.read_text())
    if (raw.get('status')!='CPU_remaining_term_converged_independent_review_pending'
            or raw.get('role') not in ('AB','B_own','B_ghost') or not raw.get('converged')
            or review.get('status')!='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending'
            or review.get('role')!=raw['role'] or stationarity.get('role')!=raw['role']
            or stationarity.get('status')!='remaining_stationarity_failed_preserved'
            or stationarity.get('passed') is not False
            or stationarity.get('orbital_gradient_threshold_au')!=1e-7):
        raise ValueError('Require a completed checkpoint with an independent stationarity failure')
    error=stationarity.get('orbital_gradient_norm_au')
    if (not isinstance(error,(int,float)) or not math.isfinite(error) or error<=1e-7
            or stationarity['energy_error_Ha']>1e-8 or stationarity['checkpoint_energy_error_Ha']>1e-8
            or abs(stationarity['reference_energy_Ha']-raw['E_Ha'])>1e-9
            or raw['source_sha256']!=launch['source_sha256']):
        raise ValueError('Prior failure is not an isolated orbital-stationarity failure')
    owner=raw['owner']
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot verify remote checkpoint owner')
    proc=Path('/proc')/str(owner['pid'])/'stat'
    if proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19]==owner['process_start']:
        raise ValueError('Prior producer is still live')
    api.verify_sources(raw['source_sha256']);api.verify_sources(review['source_sha256']);api.verify_sources(stationarity['source_sha256'])
    checkpoint=prior/raw['role']/'SCF.h5'
    if api.digest(checkpoint)!=raw['checkpoint_sha256']:raise ValueError('Prior checkpoint changed')
    for p in (prior/'launch.json',prior/'status.json',prior/raw['role']/'energy.json',checkpoint,
              prior/'grid_coordinates_AU.npy',prior/'grid_weights.npy'):
        if (review['source_sha256'].get(str(p))!=api.digest(p)
                or stationarity['source_sha256'].get(str(p))!=api.digest(p)):
            raise ValueError('Independent reviews do not bind the exact prior term')
    return raw,launch,dict(checkpoint=str(checkpoint),checkpoint_sha256=api.digest(checkpoint),
                          prior_attempt=str(prior),prior_status_sha256=api.digest(prior/'status.json'))


def bind_worker(prior,checkpoint_review,stationarity_review):
    raw,launch,previous=prior_checkpoint(prior,checkpoint_review,stationarity_review)
    module=api.load_module('same_model_counterpoise_polish',ORIGINAL)
    original_inputs=module.inputs
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),ORIGINAL,checkpoint_review,stationarity_review,
              prior/'launch.json',prior/'status.json',prior/raw['role']/'SCF.h5')}
    def inputs(*args):
        plan,sources,state,versions=original_inputs(*args)
        api.verify_sources(proof)
        return plan,{**sources,**proof},state,versions
    module.inputs=inputs
    body=inspect.getsource(module.run)
    old="mf.conv_tol=core.SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=core.SETTINGS['SCF_conv_tol_grad'];"
    new="mf.conv_tol=_targets['SCF_conv_tol_Ha'];mf.conv_tol_grad=_targets['SCF_conv_tol_grad_au'];"
    if body.count(old)!=1:raise ValueError('Counterpoise SCF template changed')
    body=body.replace(old,new,1)
    old="primary_candidate_admission=False,full_counterpoise_computed=False,qualified=False)"
    new="execution_refinement=dict(_targets),primary_candidate_admission=False,full_counterpoise_computed=False,qualified=False)"
    if body.count(old)!=1:raise ValueError('Counterpoise result schema changed')
    namespace=dict(vars(module),_targets=TARGETS);local={}
    exec(compile(body.replace(old,new,1),__file__,'exec'),namespace,local)
    module.run=local['run']
    return module,raw,launch,previous


def launch(prior,checkpoint_review,stationarity_review,out,threads):
    if out.exists():raise FileExistsError('Preserve previous refinement attempts')
    module,state,old,previous=bind_worker(prior,checkpoint_review,stationarity_review)
    binder=api.load_module('polish_original_argument_reader',ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py').argument
    command=old['command']
    names=('case','geometry-review','own','own-stationarity','ghost','ghost-stationarity')
    paths=[binder(command,'--'+name).resolve() for name in names]
    _,sources,_,versions=module.inputs(*paths)
    with allocation_guard():
        capacity=host_capacity(threads,minimum_ram_gib=112);cpus=cpu_allocation(threads)
        out.mkdir(parents=True);(out/'scratch').mkdir()
        argv=[sys.executable,str(Path(__file__).resolve()),'run','--prior',str(prior),
              '--checkpoint-review',str(checkpoint_review),'--stationarity-review',str(stationarity_review),
              '--out',str(out),'--threads',str(threads)]
        worker_command=[sys.executable,str(Path(__file__).resolve()),'run','--role',state['role']]
        for name,value in zip(names,paths):worker_command+=['--'+name,str(value)]
        worker_command+=['--raw',str(out)]
        record=dict(status='CPU_checkpoint_refinement_launch_prepared',role=state['role'],case=str(paths[0]),
                    threads=threads,cpu_affinity=cpus,capacity=capacity,source_sha256=sources,interpreter=sys.executable,
                    command=worker_command,actual_command=argv,versions=versions,resume_from=previous,
                    execution_refinement=TARGETS,maximum_host_RSS_GiB=32,full_counterpoise_computed=False,qualified=False)
        atomic_write_json(out/'launch.json',record)
        env=dict(os.environ,PYSCF_TMPDIR=str(out/'scratch'))
        for key in module.CAPS:env[key]=str(threads)
        with (out/'worker.log').open('xb') as log:
            process=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,
                                     start_new_session=True,env=env,preexec_fn=lambda:os.sched_setaffinity(0,set(cpus)))
        owner=dict(host=os.uname().nodename,pid=process.pid,
                   process_start=Path('/proc/'+str(process.pid)+'/stat').read_text().rsplit(')',1)[1].split()[19])
        atomic_write_json(out/'launcher.json',dict(owner=owner,CPU_affinity=cpus,started=time.time()))
    return owner


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','run'))
    for name in ('prior','checkpoint-review','stationarity-review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    values=[getattr(a,key).resolve() for key in ('prior','checkpoint_review','stationarity_review','out')]
    if not 1<=a.threads<=80 or any(not v.is_relative_to(ROOT/'calcs') for v in values):
        raise ValueError('Use durable project paths and a valid thread count')
    if a.action=='launch':print(json.dumps(launch(*values,a.threads)))
    else:
        module,state,old,_=bind_worker(*values[:3])
        binder=api.load_module('polish_run_arguments',ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py').argument
        names=('case','geometry-review','own','own-stationarity','ghost','ghost-stationarity')
        paths=[binder(old['command'],'--'+name).resolve() for name in names]
        module.run(state['role'],*paths,values[3],a.threads)
