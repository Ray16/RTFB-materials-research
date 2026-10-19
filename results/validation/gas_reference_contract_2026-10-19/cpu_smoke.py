"""Small actual CPU solver control; no GPU, reference fitting or bulk recomputation."""
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np
from pyscf import dft,gto
from redox.core.io import atomic_write_json
from redox.qm.gas_reference import energy_Ha,optimized_molecule
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

out=Path(__file__).resolve().parent/'cpu_smoke'
out.mkdir(exist_ok=False)
with allocation_guard():
    resources=host_capacity(1);cpus=cpu_allocation(1);os.sched_setaffinity(0,cpus)
assert all(os.environ[k]=='1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'))
def build(distance):
    mol=gto.M(atom=f'H 0 0 0; H 0 0 {distance}',basis='sto-3g',verbose=0)
    mf=dft.RKS(mol).density_fit();mf.xc='pbe';mf.conv_tol=1e-10;mf.max_cycle=100
    return mf
start=time.time();mf=build(1.0);optimized=optimized_molecule(mf,maxsteps=50)
final=dft.RKS(optimized).density_fit();final.xc='pbe';final.conv_tol=1e-10
energy=energy_Ha(final);gradient=final.nuc_grad_method();gradient.grid_response=True;g=np.asarray(gradient.kernel())
assert np.isfinite(energy) and final.converged and abs(g).max()<.00045
failure_refused=False
try:optimized_molecule(build(2.0),maxsteps=1)
except RuntimeError:failure_refused=True
assert failure_refused,'Stopped geometry escaped the convergence gate'
report=dict(status='actual_CPU_solver_controls_passed',system='H2 / PBE / STO-3G / density fitting / gas',energy_Ha=energy,
    optimized_bond_A=float(np.linalg.norm(optimized.atom_coords(unit='Angstrom')[0]-optimized.atom_coords(unit='Angstrom')[1])),
    maximum_gradient_Ha_Bohr=float(abs(g).max()),SCF_converged=bool(final.converged),stopped_optimizer_refused=failure_refused,
    CPU_affinity=cpus,resources=resources,GPU_used=False,experimental_accuracy_validated=False,
    elapsed_seconds=time.time()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
      for p in (Path(__file__),Path(__import__('redox.qm.gas_reference',fromlist=['']).__file__))})
atomic_write_json(out/'report.json',report)
print(json.dumps(report,indent=2))
