"""One reserved GPU serves an atomic finite queue of distinct state calculations."""
import gc
import json
import os
from pathlib import Path
import runpy
import signal
import socket
import sys
import time

from redox.core.io import atomic_write_json
from redox.qm.uma import file_sha256
from redox.workflows.resources import verify_gate_owner

base=Path(__file__).resolve().parent
plan=json.loads((base/'plan.json').read_text())
verify_gate_owner()
assert all(file_sha256(p)==h for p,h in plan['producer_sha256'].items())
from redox.qm.native_thermal import implementation_controls
implementation_controls()  # Fail environment/certificate problems before claiming any state.
owner=dict(host=socket.gethostname(),pid=os.getpid(),
    process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
    cpus=sorted(os.sched_getaffinity(0)),gpu=int(os.environ['CUDA_VISIBLE_DEVICES']))
path=base/'workers'/f"{owner['host']}_gpu{owner['gpu']}.json"
report=dict(owner,status='claiming',items=[],started=time.time())
atomic_write_json(path,report)
consecutive=0


def stop(signum,frame):
    raise KeyboardInterrupt(f'signal {signum}')


signal.signal(signal.SIGTERM,stop)
try:
    for index,task in enumerate(plan['tasks']):
        key=task['key']
        if (base/'done'/f'{key}.json').exists() or (base/'failures'/f'{key}.json').exists():
            continue
        claim=base/'claims'/key
        try:claim.mkdir()
        except FileExistsError:continue
        atomic_write_json(claim/'owner.json',owner)
        try:
            # A prior attempt cannot be silently recycled when a claim disappears.
            if (base/'tasks'/key).exists():
                raise RuntimeError('Existing attempt has no completion receipt; inspect original owner')
            report.update(status='computing',active_key=key)
            atomic_write_json(path,report)
            sys.argv=[str(base/'worker.py'),'--task',str(index)]
            try:
                # Same process keeps verify_gate_owner's direct parent invariant.
                runpy.run_path(str(base/'worker.py'),run_name='__main__')
            except Exception as exc:
                status_path=base/'tasks'/key/'status.json'
                state=json.loads(status_path.read_text()) if status_path.exists() else {}
                if state.get('stage') not in ('dft','active_sp','native_thermal'):
                    raise RuntimeError('Pre-execution failure; stop before consuming further queue items') from exc
                failure=dict(key=key,owner=owner,error=repr(exc),state_status=str(status_path),attempts=1,finished=time.time())
                atomic_write_json(base/'failures'/f'{key}.json',failure)
                report['items'].append(dict(key=key,status='failed_preserved'))
                consecutive+=1
            else:
                status_path=base/'tasks'/key/'status.json'
                state=json.loads(status_path.read_text())
                assert state['status']=='state_records_verified'
                atomic_write_json(base/'done'/f'{key}.json',dict(key=key,owner=owner,status_path=str(status_path),status_sha256=file_sha256(status_path),finished=time.time()))
                report['items'].append(dict(key=key,status='state_records_verified'))
                consecutive=0
        finally:
            gc.collect()
            from redox.qm.dft import _free_gpu_pool
            _free_gpu_pool()
            (claim/'owner.json').unlink()
            claim.rmdir()
            atomic_write_json(path,report)
        if consecutive>=3:
            raise RuntimeError('Three consecutive actual failures; stop GPU worker for review')
    report.update(status='no_unclaimed_unattempted_states',active_key=None)
except BaseException as exc:
    report.update(status='stopped_preserving_outputs',error=repr(exc))
    raise
finally:
    report['finished']=time.time()
    atomic_write_json(path,report)
