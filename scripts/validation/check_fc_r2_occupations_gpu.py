#!/usr/bin/env python
"""Bounded GPU occupation check of the native Fc+ solution surface."""
import hashlib,json,os,socket,time,sys
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto,lib
from redox.core.common import ROOT,DFT,HARTREE_EV
from redox.qm.dft import _build_smd_mf,_orbital_diag
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
def main():
 assert os.getenv("CUDA_VISIBLE_DEVICES") and sorted(os.sched_getaffinity(0))==[10,11] and os.getloadavg()[0]+2<os.cpu_count()
 lib.num_threads(2);xyz=DFT/'ferrocene/ox/opt.xyz';atoms=read(xyz);sha=hashlib.sha256(xyz.read_bytes()).hexdigest()
 status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],cpus=[10,11],gpu=os.getenv("CUDA_VISIBLE_DEVICES"),geometry_sha256=sha,started=time.time(),guesses=[],source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),protocol=dict(xc='r2scan',disp='d4',basis='def2-svp',solvent='acetonitrile',solvent_model='SMD',density_fit=True,grid_level=3,conv_tol=1e-12,conv_tol_grad=1e-8))
 def save():
  p=BASE/'fc_r2_occupation_gpu_status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
 save()
 atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist()))
 parent=gto.M(atom=atom,charge=0,spin=0,basis='def2-svp',verbose=0,max_memory=6000)
 ref=_build_smd_mf(parent,'r2scan','acetonitrile','gpu',disp='d4');ref.grids.level=3;ref.conv_tol=1e-12;ref.conv_tol_grad=1e-8;ref.max_cycle=200;ref.kernel()
 assert ref.converged
 import cupy
 C=ref.mo_coeff;occ=ref.mo_occ.get()/2;nocc=int(np.count_nonzero(occ));mol=gto.M(atom=atom,charge=1,spin=1,basis='def2-svp',verbose=0,max_memory=6000)
 seeds={}
 for k in range(2):
  beta=occ.copy();beta[nocc-1-k]=0
  mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');seeds['rm'+str(k)]=mf.make_rdm1((C,C),(cupy.asarray(occ),cupy.asarray(beta)))
 del ref
 for name in ['minao','atom','huckel','rm0','rm1']:
  mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
  if name in seeds:mf.kernel(dm0=seeds[name])
  else:mf.init_guess=name;mf.kernel()
  record=dict(guess=name,converged=bool(mf.converged),e_Ha=float(mf.e_tot),e_eV=float(mf.e_tot)*HARTREE_EV,spin_diagnostics=_orbital_diag(mf,1))
  status['guesses'].append(record);save();print('R2_GPU_GUESS',record,flush=True)
  del mf
 verified=[r for r in status['guesses'] if r['converged'] and abs(r['spin_diagnostics'].get('spin_contam',999))<.1]
 assert verified and sha==hashlib.sha256(xyz.read_bytes()).hexdigest()
 canonical=json.loads((BASE/'fc_conservative_minimum_refinement/status.json').read_text())['selected']['electronic_opt_energy_Ha']
 best=min(verified,key=lambda r:r['e_Ha']);delta=(best['e_Ha']-canonical)*HARTREE_EV
 status.update(best=best,native_source_opt_energy_Ha=canonical,best_minus_native_surface_eV=delta,
  status='native_surface_occupation_confirmed' if abs(delta)<1e-5 else 'native_surface_lower_occupation_found' if delta<0 else 'native_surface_energy_disagreement',finished=time.time(),gpu_imported=any(k in sys.modules for k in ('cupy','gpu4pyscf')))
 save()
if __name__=='__main__':main()
