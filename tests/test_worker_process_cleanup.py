"""A queue must survive a state boundary that killed its cached CPU kernel."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


def test_next_state_can_use_thermodynamics_after_owned_child_cleanup():
    from redox.qm import qrrho
    if not qrrho.LIBRARY.exists():pytest.skip('Qualified local xTB library unavailable')
    script=r'''
import json,subprocess,sys
from redox.qm import qrrho
from redox.workflows.processes import stop_owned_children
first=qrrho.atomic_masses([1,6,8]).tolist()
old=qrrho._WORKER['process'];old.kill();old.wait(timeout=5)
child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(120)'])
stop_owned_children(timeout=5)
assert child.poll() is not None
second=qrrho.atomic_masses([1,6,8]).tolist()
assert qrrho._WORKER['process'].pid!=old.pid and second==first
stop_owned_children(timeout=5)
assert qrrho._WORKER is None
print(json.dumps(dict(reused_client_recovered=True,owned_child_joined=True)))
'''
    env=dict(os.environ)
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    result=subprocess.run([sys.executable,'-c',script],env=env,text=True,capture_output=True,timeout=30,check=True)
    assert json.loads(result.stdout)==dict(reused_client_recovered=True,owned_child_joined=True)
