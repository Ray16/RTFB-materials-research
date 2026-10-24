"""CPU resource configuration must leave the shared Python API intact on exit."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize('fail', [False,True])
def test_CPU_context_restores_functions_after_success_or_failure(tmp_path, fail):
    root=Path(__file__).resolve().parents[1]
    code=r'''
import importlib.util,sys
from pathlib import Path
from pyscf import gto
from redox.qm import dft
from redox.core import protocol
path=Path(sys.argv[1]);folder=Path(sys.argv[2]);fail=sys.argv[3]=='True'
spec=importlib.util.spec_from_file_location('CPU_execution_test',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
before=(gto.M,dft._kernel_robust,dft._free_gpu_pool,protocol.software_versions)
try:
 with module.cpu_execution(folder,1024):
  mol=gto.M(atom='H 0 0 0; H 0 0 .74',basis='sto-3g',verbose=0)
  assert mol.max_memory==1024
  if fail:raise RuntimeError('intentional CPU scope exit')
except RuntimeError:
 assert fail
assert before==(gto.M,dft._kernel_robust,dft._free_gpu_pool,protocol.software_versions)
assert 'cupy' not in sys.modules and 'gpu4pyscf' not in sys.modules
'''
    result=subprocess.run([sys.executable,'-c',code,str(root/'scripts/pipeline/cpu_quantum_execution.py'),
        str(tmp_path/'execution'),str(fail)],cwd=root,env=dict(os.environ,PYTHONPATH=str(root/'src')),
        capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
