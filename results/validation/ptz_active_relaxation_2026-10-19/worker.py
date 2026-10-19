"""Private active-level geometry sensitivity; no canonical or thermal adoption."""
import argparse,hashlib,json,os,socket,time,traceback,signal
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read,write
from pyscf import gto
from pyscf.geomopt.geometric_solver import kernel
from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
from redox.qm.uma import geometry_checks
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner
base=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--task',type=int,required=True);args=parser.parse_args()
plan=json.loads((base/'plan.json').read_text());task=plan['tasks'][args.task];out=base/task['state'];out.mkdir(exist_ok=False)
report=dict(status='starting',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),production_admission=False,thermal_recomputed=False,Hessian_minimum_proven=False,canonical_values_changed=False)
def save():atomic_write_json(out/'status.json',report)
def verify():
 for p,h in {**plan['producer_sha256'],**task['input_sha256']}.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
signal.signal(signal.SIGTERM,lambda signum,frame:(_ for _ in ()).throw(KeyboardInterrupt(f'signal {signum}')))
save()
try:
 verify_gate_owner();verify();assert len(os.sched_getaffinity(0))==2
 atoms=read(Path(task['source'])/'geometry.xyz');symbols=atoms.get_chemical_symbols();p=plan['electronic_protocol'];solver=plan['SCF'];BOHR=.52917721067;HA=27.211386245988
 def solve(coords):
  mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=task['charge'],spin=task['mult']-1,basis=p['basis'],ecp=_ecp_for(symbols,p['basis']),verbose=0,max_memory=6000)
  mf=_build_smd_mf(mol,p['xc'],p['solvent'],'gpu',disp=p['disp'],nlc=p['nlc']);mf.grids.level=3;mf.init_guess=task['guess'];mf.conv_tol=solver['conv_tol'];mf.conv_tol_grad=solver['conv_tol_grad'];mf.max_cycle=solver['max_cycle']
  e=float(mf.kernel());assert mf.converged and np.isfinite(e)
  return mf,e
 coords=atoms.positions/BOHR;mf,start=solve(coords);assert abs(start*HA-task['expected_eV'])<.001
 report.update(status='optimizing',initial_energy_Ha=start,initial_reproduction_difference_eV=start*HA-task['expected_eV'],physical_gpu=os.environ['CUDA_VISIBLE_DEVICES']);save()
 grad=mf.nuc_grad_method();grad.grid_response=True;scanner=grad.as_scanner()
 def checkpoint(env):
  xyz=np.asarray(env['coords']);g=np.asarray(env['gradients'].get() if hasattr(env['gradients'],'get') else env['gradients'])
  write(out/'checkpoint.xyz',Atoms(symbols,positions=xyz*BOHR),format='xyz')
  report.update(cycle=int(env['self'].cycle),last_energy_Ha=float(env['energy']),last_max_gradient_Ha_Bohr=float(abs(g).max()),last_progress=time.time());save()
 converged,mol=kernel(scanner,assert_convergence=True,callback=checkpoint,**plan['optimizer'])
 assert converged and scanner.converged
 final_coords=mol.atom_coords(unit='Bohr').copy();optimizer_energy=report['last_energy_Ha'];assert np.max(abs(final_coords-read(out/'checkpoint.xyz').positions/BOHR))<1e-5
 write(out/'relaxed.xyz',Atoms(symbols,positions=final_coords*BOHR),format='xyz')
 assert not geometry_checks(task['smiles'],out/'relaxed.xyz',task['charge'],task['mult'])
 del mf,grad,scanner;_free_gpu_pool();mf,e=solve(final_coords);assert abs(e-optimizer_energy)*HA<.001;report.update(optimizer_final_energy_Ha=optimizer_energy,cold_reproduction_difference_eV=(e-optimizer_energy)*HA);gradient=mf.nuc_grad_method();gradient.grid_response=True;raw=gradient.kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw).copy();np.save(out/'final_gradient_Ha_Bohr.npy',g)
 assert np.isfinite(g).all() and abs(g).max()<=plan['stationarity_max_Ha_Bohr']
 column=int(abs(g).argmax());step=plan['FD_step_Bohr'];energies=[];del mf,gradient,raw;_free_gpu_pool()
 report.update(status='checking_cold_final_gradient',optimized=True,geometry_optimization_converged=True,final_energy_Ha=e,electronic_relaxation_eV=(e-start)*HA,final_max_gradient_Ha_Bohr=float(abs(g).max()),FD_column=column);save()
 for sign in (-1,1):
  x=final_coords.copy();x.flat[column]+=sign*step;mf,value=solve(x);energies.append(value);del mf;_free_gpu_pool()
 fd=(energies[1]-energies[0])/(2*step);error=abs(fd-g.flat[column]);assert error<=plan['FD_absolute_acceptance_Ha_Bohr']
 assert (e-start)*HA<=.001;verify()
 report.update(status='active_level_stationary_geometry_verified',finite_energy_Ha=energies,FD_derivative_Ha_Bohr=fd,FD_absolute_error_Ha_Bohr=float(error),output_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [out/'relaxed.xyz',out/'final_gradient_Ha_Bohr.npy']})
except BaseException as exc:
 report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
finally:
 report['finished']=time.time();save()
