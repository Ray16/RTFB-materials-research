"""Verify original process identities are terminal on the execution host."""
import json
import socket
import time
from pathlib import Path

import psutil

from redox.core.io import atomic_write_json

base = Path(__file__).resolve().parent
launch = json.loads((base / 'launch.json').read_text())
plan = json.loads((base / 'plan.json').read_text())
assert socket.gethostname().split('.')[0] == 'lambda13'
processes = []
controllers = []


def inspect(pid, start, role):
    stat = Path(f'/proc/{pid}/stat')
    try:
        fields = stat.read_text().rsplit(')', 1)[1].split()
        same_start = fields[19] == str(start)
        state = fields[0]
    except FileNotFoundError:
        same_start, state = False, None
    terminal = not same_start or state == 'Z'
    processes.append(dict(pid=pid, role=role, same_start=same_start,
                          state=state, terminal=terminal))
    assert terminal, processes[-1]


for worker in launch['workers']:
    inspect(worker['controller_pid'], worker['process_start'], 'controller')
    execution_path = base / f"reserved_{worker['gpu']}" / 'execution.json'
    execution = json.loads(execution_path.read_text())
    assert execution['status'] == 'process_completed' and execution['returncode'] == 0
    gate_pid = execution['gate_pid']
    if psutil.pid_exists(gate_pid):
        process = psutil.Process(gate_pid)
        # A reused PID is not the original gate. A zombie is already terminal.
        assert process.status() == psutil.STATUS_ZOMBIE or str(base / 'worker.py') not in process.cmdline()
    controllers.append(dict(gpu=worker['gpu'], execution=str(execution_path),
                            gate_pid=gate_pid, returncode=execution['returncode']))
    for index in worker['tasks']:
        task = plan['tasks'][index]
        record = json.loads((base / 'tasks' / task['key'] / 'status.json').read_text())
        assert record['status'] == 'fixed_geometry_electronic_method_point_verified'
        inspect(record['pid'], record['process_start'], task['key'])
    owner = Path('/nfs/lambda_stor_01/homes/rzhu/.gpu_locks') / f"lambda13_gpu{worker['gpu']}" / 'owner'
    if owner.exists():
        assert owner.read_text().split()[0] != str(gate_pid), 'Original reservation remains'

report = dict(status='all_nine_method_points_and_eight_controllers_verified_terminal',
              host=socket.gethostname(), processes=processes, controllers=controllers,
              checked=time.time(), scientific_review='independent_review.json')
atomic_write_json(base / 'terminal_verification.json', report)
print(report['status'])
