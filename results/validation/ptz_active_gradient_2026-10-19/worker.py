"""Active-level force/finite-energy diagnostic; frozen input, no promotion."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import time
import traceback

import numpy as np
from ase.io import read
from pyscf import gto
from redox.core.io import atomic_write_json
from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
from redox.workflows.resources import verify_gate_owner

ap=argparse.ArgumentParser();ap.add_argument('--task',type=int,required=True);a=ap.parse_args()
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());task=plan['tasks'][a.task]
out=base/task['state'];out.mkdir(exist_ok=False)
report=dict(status='running',task=task,host=socket.gethostname(),pid=os.getpid(),
 process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
 CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),production_admission=False,
 geometry_optimized=False,thermochemistry_recomputed=False)
def save():atomic_write_json(out/'report.json',report)
save()
try:
 verify_gate_owner();assert len(os.sched_getaffinity(0))==2
 for path,h in {**plan['producer_sha256'],**task['input_sha256']}.items():
  assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h,path
 atoms=read(Path(task['source'])/'geometry.xyz');symbols=atoms.get_chemical_symbols()
 coords=atoms.positions/.52917721067;p=plan['protocol'];solver=plan['validation_solver']
 def solve(x):
  mol=gto.M(atom=list(zip(symbols,x.tolist())),unit='Bohr',charge=task['charge'],spin=task['mult']-1,
   basis=p['basis'],ecp=_ecp_for(symbols,p['basis']),verbose=0,max_memory=6000)
  mf=_build_smd_mf(mol,p['xc'],p['solvent'],'gpu',disp=p['disp'],nlc=p['nlc'])
  mf.grids.level=3;mf.init_guess=task['guess'];mf.conv_tol=solver['conv_tol']
  mf.conv_tol_grad=solver['conv_tol_grad'];mf.max_cycle=solver['max_cycle']
  e=float(mf.kernel())
  if not mf.converged or not np.isfinite(e):raise RuntimeError('Gradient-control SCF failed')
  return mf,e
 mf,e=solve(coords);delta=e*27.211386245988-task['expected_eV']
 report.update(energy_Ha=e,energy_reproduction_difference_eV=delta);save()
 if abs(delta)>.001:raise RuntimeError('Active energy not reproduced within 1 meV')
 grad=mf.nuc_grad_method();grad.grid_response=True;raw=grad.kernel()
 g=np.asarray(raw.get() if hasattr(raw,'get') else raw).copy();assert np.isfinite(g).all()
 np.save(out/'gradient_Ha_Bohr.npy',g)
 column=int(np.argmax(abs(g)));step=plan['FD_step_Bohr']
 report.update(max_gradient_Ha_Bohr=float(abs(g).max()),RMS_gradient_Ha_Bohr=float(np.sqrt(np.mean(g*g))),
               checked_column=column,analytic_component_Ha_Bohr=float(g.flat[column]),status='checking_finite_energy')
 save();del mf,grad,raw;_free_gpu_pool();energies=[]
 for sign in (-1,1):
  x=coords.copy();x.flat[column]+=sign*step;mf,value=solve(x);energies.append(value)
  del mf;_free_gpu_pool()
 numerical=(energies[1]-energies[0])/(2*step);error=abs(numerical-g.flat[column])
 report.update(finite_energy_components_Ha=energies,finite_energy_derivative_Ha_Bohr=numerical,
               finite_energy_absolute_error_Ha_Bohr=float(error))
 save()
 if error>plan['FD_absolute_acceptance_Ha_Bohr']:raise RuntimeError('NLC/SMD force finite-energy check failed')
 report.update(status='active_gradient_and_finite_energy_checked',gradient_sha256=hashlib.sha256((out/'gradient_Ha_Bohr.npy').read_bytes()).hexdigest())
except BaseException as exc:
 report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
finally:
 report['finished']=time.time();save()
