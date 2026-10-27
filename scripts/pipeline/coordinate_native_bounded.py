#!/usr/bin/env python
"""Resume native thermochemistry with bounded placement and durable errors."""
import inspect
import json
import os
from pathlib import Path
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources

SCRIPT=Path(__file__).resolve()
ORIGINAL=ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py'


def adapt(module):
    guards=api.load_module('native_bounded_host_lock',ROOT/'scripts/pipeline/coordinate_smooth_geometry_bounded.py')
    factory=module.native_placement
    def placement(profile):
        base=factory(profile);original=base.launch_local
        def launch_local(batch,campaign,gpu):
            host=os.uname().nodename.split('.')[0];directory=campaign/'placement_checks';directory.mkdir(exist_ok=True)
            destination=directory/(str(time.time_ns())+'.json')
            row=dict(phase='remote_native_placement_entered',host=host,pid=os.getpid(),GPU=gpu,
                     started=time.time(),scientific_start_inferred=False)
            atomic_write_json(destination,row)
            previous=resources.allocation_guard
            resources.allocation_guard=lambda:guards.bounded_allocation(ROOT/'logs/resource_locks'/f'{host}.lock',row,destination)
            try:
                try:result=original(batch,campaign,gpu)
                except guards.HostAllocationBusy:result=dict(status='resource_capacity_unavailable',host=host,reason='host_allocation_guard_busy')
                row.update(phase='remote_native_placement_returned',observation=result,finished=time.time());atomic_write_json(destination,row)
                return result
            except BaseException as error:
                row.update(phase='remote_native_placement_exception_preserved',error=repr(error),finished=time.time());atomic_write_json(destination,row)
                raise
            finally:resources.allocation_guard=previous
        base.launch_local=launch_local
        return base
    module.native_placement=placement;module.SCRIPT=SCRIPT
    body=inspect.getsource(module.coordinate)
    changes=[
        ('        failures, processed, observations = {}, set(), {}','        failures, processed, observations = {}, set(), {}; placement_errors = {}'),
        ('if not pending or failures.get(host,0)>=3: continue','if not pending or host in placement_errors or failures.get(host,0)>=3: continue'),
        ("                atomic_write_json(path,record)","                record['unverified_placement_hosts'] = placement_errors\n                atomic_write_json(path,record)"),
        ("                            observations[entry['key']]=dict(owner=dict(host=host,pid=-1,process_start='unknown'),error=repr(error))\n                            pending.pop(0);continue",
         "                            detail=dict(owner=dict(host=host,pid=-1,process_start='unknown'),error=repr(error),stderr=str(getattr(error,'stderr','') or '')[-4000:],checked=time.time())\n                            observations[entry['key']]=detail; placement_errors[host]=detail\n                            atomic_write_json(case/'placement_error.json',detail)\n                            pending.pop(0);break")]
    for old,new in changes:
        if body.count(old)!=1:raise ValueError('Native placement template changed: '+old[:55])
        body=body.replace(old,new,1)
    local={};exec(compile(body,__file__,'exec'),dict(vars(module)),local);module.coordinate=local['coordinate']


def main():
    os.chdir(ROOT)
    inventory=api.load_module('native_bounded_reservation_inventory',ROOT/'scripts/pipeline/coordinate_native_inventory_resilient.py')
    hosts=sys.argv[sys.argv.index('--hosts')+1].split(',') if '--hosts' in sys.argv else [os.uname().nodename.split('.')[0]]
    if len(hosts)!=len(set(hosts)) or not set(hosts)<=inventory.permitted_hosts(inventory.POLICY):raise ValueError('Distinct current-policy hosts required')
    proof={str(p):api.digest(p) for p in (SCRIPT,ORIGINAL,Path(inventory.__file__),inventory.GATE,inventory.POLICY,
        ROOT/'scripts/pipeline/coordinate_smooth_geometry_bounded.py')}
    previous=resources.free_slots
    def free_slots(requested=None):
        api.verify_sources(proof)
        chosen=hosts if requested is None else list(requested)
        if not set(chosen)<=inventory.permitted_hosts(inventory.POLICY):raise ValueError('Reservation host policy changed')
        slots,_=inventory.query_all(chosen)
        return slots
    resources.free_slots=free_slots
    module=api.load_module('bounded_native_unchanged_scientific_queue',ORIGINAL);adapt(module)
    try:module.main()
    finally:resources.free_slots=previous


if __name__=='__main__':main()
