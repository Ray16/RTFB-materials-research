"""Run an archived namespace review on its source host with one disjoint CPU."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


def command_for_host(command, host):
    if host.split('.')[0] == os.uname().nodename.split('.')[0]:
        return command
    from redox.core.common import ROOT
    # Allocate before entering the RO archive namespace, where the resource-lock
    # directory intentionally cannot be changed. No GPU or physical gate changes.
    code = ('import os,subprocess; '
            'from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity; '
            'guard=allocation_guard(); guard.__enter__(); host_capacity(1); '
            'os.sched_setaffinity(0,cpu_allocation(1)); guard.__exit__(None,None,None); '
            'raise SystemExit(subprocess.run('+repr(command)+').returncode)')
    prefix = ['env', 'OMP_NUM_THREADS=1', 'OMP_THREAD_LIMIT=1', 'MKL_NUM_THREADS=1',
              'OPENBLAS_NUM_THREADS=1', 'NUMEXPR_NUM_THREADS=1', 'XTB_THREADS=1',
              'PYTHONPATH='+str(ROOT/'src'), sys.executable, '-c', code]
    return ['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host.split('.')[0], shlex.join(prefix)]
