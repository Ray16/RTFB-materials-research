#!/usr/bin/env python
"""Keep remote geometry placement within its CPU allocation/RPC budget.

Retains the source-bound inventory and scientific producer of the prior adapter.
A busy host allocation lock is a resource wait, and launch phase receipts expose
where any later placement fails. Unverified RPC hosts stop receiving new inputs.
"""
from contextlib import contextmanager
import fcntl
import inspect
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources

SCRIPT=Path(__file__).resolve()
BASE=ROOT/'scripts/pipeline/coordinate_smooth_geometry_resilient.py'


class HostAllocationBusy(RuntimeError):
    pass


@contextmanager
def bounded_allocation(path, record, destination):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            record.update(phase='host_allocation_guard_busy',checked=time.time())
            atomic_write_json(destination,record)
            raise HostAllocationBusy('Host CPU allocation guard is busy')
        record.update(phase='host_allocation_guard_acquired',checked=time.time())
        atomic_write_json(destination,record)
        yield


def adapt(module,profile,base_adapt):
    base_adapt(module,profile)
    module.SCRIPT=SCRIPT
    original_launch=module.launch_local

    def launch_local(batch,campaign,gpu):
        directory=campaign/'placement_checks';directory.mkdir(exist_ok=True)
        destination=directory/(str(time.time_ns())+'.json')
        host=os.uname().nodename.split('.')[0]
        record=dict(phase='remote_placement_entered',host=host,pid=os.getpid(),GPU=gpu,
                    process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                    started=time.time(),scientific_start_inferred=False,
                    source_sha256={str(p):api.digest(p) for p in (SCRIPT,BASE)})
        atomic_write_json(destination,record)
        original_guard=resources.allocation_guard
        def guard():
            return bounded_allocation(ROOT/'logs/resource_locks'/f'{host}.lock',record,destination)
        resources.allocation_guard=guard
        try:
            try:result=original_launch(batch,campaign,gpu)
            except HostAllocationBusy:
                result=dict(status='resource_capacity_unavailable',host=host,reason='host_allocation_guard_busy')
            record.update(phase='remote_placement_returned',observation=result,finished=time.time())
            atomic_write_json(destination,record)
            return result
        except BaseException as error:
            record.update(phase='remote_placement_exception_preserved',error=repr(error),finished=time.time())
            atomic_write_json(destination,record)
            raise
        finally:
            resources.allocation_guard=original_guard
    module.launch_local=launch_local
    original=api.load_module('bounded_original_geometry_queue',ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    body=inspect.getsource(original.coordinate)
    replacements=[
      ('stationary_geometry_checked_active_parent_and_native_minimum_pending',
       'experimental_smooth_geometry_checked_energy_and_native_QC_pending'),
      ('        failures_by_host = {}; processed = set(); assignments = {}',
       '        failures_by_host = {}; placement_failures_by_host = {}; processed = set(); assignments = {}'),
      ("                atomic_write_json(status_path, record)",
       "                record['unverified_placement_hosts'] = placement_failures_by_host\n                atomic_write_json(status_path, record)"),
      ('if not pending or failures_by_host.get(host, 0) >= 3 or not preflight(batch, host):',
       'if not pending or host in placement_failures_by_host or failures_by_host.get(host, 0) >= 3 or not preflight(batch, host):'),
      ("                            assignments[entry['key']] = dict(owner=dict(host=host, pid=-1, process_start='unknown'), error=repr(error))\n                            pending.pop(0); continue",
       "                            detail = dict(owner=dict(host=host, pid=-1, process_start='unknown'), error=repr(error), stderr=str(getattr(error, 'stderr', '') or '')[-4000:], checked=time.time(), status='placement_response_unverified_preserved')\n                            assignments[entry['key']] = detail\n                            placement_failures_by_host[host] = detail\n                            atomic_write_json(Path(entry['campaign']) / 'placement_error.json', detail)\n                            pending.pop(0); break")]
    for old,new in replacements:
        if body.count(old)!=1:raise ValueError('Geometry placement template changed: '+old[:60])
        body=body.replace(old,new,1)
    local={};namespace=dict(vars(module),SCRIPT=SCRIPT)
    exec(compile(body,__file__,'exec'),namespace,local)
    module.coordinate=local['coordinate']


def main():
    base=api.load_module('source_bound_geometry_inventory',BASE)
    original_adapt=base.adapt
    base.adapt=lambda module,profile:adapt(module,profile,original_adapt)
    base.main()


if __name__=='__main__':main()
