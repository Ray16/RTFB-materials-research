"""Private RKS reproduction / UKS zero-spin orbital-seed diagnostics; never promote."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import time
import traceback

import numpy as np
from ase.io import read
from pyscf import gto
from redox.core.io import atomic_write_json
from redox.qm.dft import _build_smd_mf, _ecp_for, _orbital_diag
from redox.workflows.resources import verify_gate_owner

ap=argparse.ArgumentParser();ap.add_argument('--baseline',action='store_true');ap.add_argument('--seed',type=int)
a=ap.parse_args();base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text())
assert a.baseline != (a.seed is not None)
out=base/('baseline' if a.baseline else f'seed_{a.seed}');out.mkdir(exist_ok=False)
report=dict(status='running',host=socket.gethostname(),pid=os.getpid(),
 process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
 CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),production_admission=False,
 quantum_geometry_optimized=False,thermochemistry_recomputed=False)
def save():atomic_write_json(out/'report.json',report)
save()
try:
 verify_gate_owner();assert len(os.sched_getaffinity(0))==2
 for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
  assert os.environ[key]=='2'
 for path,digest in plan['producer_sha256'].items():
  assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
 assert hashlib.sha256((base/'geometry.xyz').read_bytes()).hexdigest()==plan['source_geometry_sha256']
 assert hashlib.sha256((base/'selected_SP.json').read_bytes()).hexdigest()==plan['source_SP_sha256']
 atoms=read(base/'geometry.xyz');symbols=atoms.get_chemical_symbols();p=plan['protocol']
 mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),basis=p['basis'],charge=2,spin=0,
           ecp=_ecp_for(symbols,p['basis']),verbose=0,max_memory=6000)
 import cupy as cp
 if a.baseline:
  mf=_build_smd_mf(mol,p['xc'],p['solvent'],'gpu',disp=p['disp'],nlc=p['nlc'])
  dm=None
 else:
  control=json.loads((base/'baseline/report.json').read_text())
  assert control['status']=='restricted_baseline_reproduced'
  frozen=base/'baseline/orbitals.npz'
  assert hashlib.sha256(frozen.read_bytes()).hexdigest()==control['orbitals_sha256']
  data=np.load(frozen);coeff=data['mo_coeff'];occupied=np.flatnonzero(data['mo_occ']>0)
  virtual=np.flatnonzero(data['mo_occ']==0);h,l=int(occupied[-1]),int(virtual[0])
  theta=plan['seed_angles_radians'][a.seed];densities=[]
  for sign in (1,-1):
   c=coeff[:,occupied].copy()
   c[:,-1]=math.cos(theta)*coeff[:,h]+sign*math.sin(theta)*coeff[:,l]
   densities.append(c@c.conj().T)
  electrons=[float(np.trace(d@data['overlap']).real) for d in densities]
  assert all(abs(n-mol.nelectron/2)<1e-6 for n in electrons)
  from gpu4pyscf import dft
  mf=dft.UKS(mol).density_fit().SMD();mf.xc=p['xc'];mf.nlc=p['nlc'];mf.with_solvent.solvent=p['solvent']
  dm=cp.asarray(np.asarray(densities));report.update(seed_angle_radians=theta,initial_spin_electron_counts=electrons)
 mf.conv_tol=p['conv_tol'];mf.max_cycle=p['max_cycle'];mf.grids.level=3
 energy=float(mf.kernel() if dm is None else mf.kernel(dm0=dm))
 if not mf.converged or not math.isfinite(energy):raise RuntimeError('Probe SCF did not converge')
 energy_eV=energy*27.211386245988;diag=_orbital_diag(mf,0)
 report.update(energy_Ha=energy,energy_eV=energy_eV,SCF_converged=True,orbital_diagnostics=diag,
               model='RKS' if a.baseline else 'UKS with zero net spin; diagnostic only',protocol=p)
 if a.baseline:
  delta=energy_eV-plan['expected_RKS_eV'];report['difference_from_selected_SP_eV']=delta
  if abs(delta)>.001:raise RuntimeError('Restricted baseline does not reproduce selected active SP within 1 meV')
  np.savez(out/'orbitals.npz',mo_coeff=mf.mo_coeff.get(),mo_occ=mf.mo_occ.get(),overlap=mf.get_ovlp().get())
  report.update(status='restricted_baseline_reproduced',orbitals_sha256=hashlib.sha256((out/'orbitals.npz').read_bytes()).hexdigest())
 else:
  report.update(status='unrestricted_seed_probe_complete',difference_from_RKS_eV=energy_eV-control['energy_eV'],
                pure_singlet_not_assumed=True)
except BaseException as exc:
 report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
finally:
 report['finished']=time.time();save()
