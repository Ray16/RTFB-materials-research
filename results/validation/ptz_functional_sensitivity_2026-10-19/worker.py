"""Predeclared vertical electronic-method sensitivity at fixed private geometries."""
import argparse,fcntl,hashlib,json,os,socket,time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from redox.qm.dft import _build_smd_mf,_ecp_for,_orbital_diag,_free_gpu_pool
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner
base=Path(__file__).resolve().parent;p=argparse.ArgumentParser();p.add_argument('--tasks',required=True);args=p.parse_args();plan=json.loads((base/'plan.json').read_text());verify_gate_owner()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for index in map(int,args.tasks.split(',')):
 task=plan['tasks'][index];out=base/'tasks'/task['key'];out.mkdir(parents=True,exist_ok=True);lock=(out/'writer.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 report=dict(status='starting',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),production_admission=False,thermal_computed=False)
 def save():atomic_write_json(out/'status.json',report)
 save()
 try:
  for path,h in {**plan['producer_sha256'],**task['input_sha256']}.items():assert sha(path)==h,path
  atoms=read(task['geometry']);coords=atoms.positions/.52917721067;symbols=atoms.get_chemical_symbols();method=task['method'];solver=plan['SCF']
  mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=task['charge'],spin=task['mult']-1,basis=plan['basis'],ecp=_ecp_for(symbols,plan['basis']),verbose=0,max_memory=6000)
  mf=_build_smd_mf(mol,method['xc'],'acetonitrile','gpu',disp=method['disp'],nlc=method['nlc']);mf.grids.level=3;mf.init_guess=task['guess'];mf.conv_tol=solver['conv_tol'];mf.conv_tol_grad=solver['conv_tol_grad'];mf.max_cycle=solver['max_cycle']
  report.update(status='computing',physical_gpu=os.environ['CUDA_VISIBLE_DEVICES']);save()
  energy=float(mf.kernel());assert mf.converged and np.isfinite(energy)
  diag=_orbital_diag(mf,task['mult']-1);assert all(isinstance(diag.get(k),(int,float)) and np.isfinite(diag[k]) for k in ['s_squared','spin_contam','homo_eV']);assert abs(diag['spin_contam'])<=.1
  difference=energy*27.211386245988-task['WB_control_eV'] if method['id']=='wb97mv' else None
  if difference is not None:assert abs(difference)<.001
  for path,h in task['input_sha256'].items():assert sha(path)==h,path
  report.update(status='fixed_geometry_electronic_method_point_verified',energy_Ha=energy,energy_eV=energy*27.211386245988,SCF_converged=True,spin_diagnostics=diag,WB_control_difference_eV=difference)
  del mf;_free_gpu_pool()
 except BaseException as exc:report.update(status='failed_preserving_outputs',error=repr(exc));raise
 finally:report['finished']=time.time();save();lock.close()
