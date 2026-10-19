"""Isolated 1e-6-Bohr CDS derivative qualification; never admit production records.

The physical protocol and acceptance thresholds remain unchanged. Two complete
baseline replays and a difficult target need independent review before any loader
or production execution path can use this derivative setting.
"""
import argparse,hashlib,json,os,shutil,socket,time,traceback
from pathlib import Path
import numpy as np
from redox.core.common import ROOT,DFT
from redox.qm import native_thermal as N
from redox.workflows.resources import verify_gate_owner
from redox.validation.spins import atomic_json

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--manifest',required=True,type=Path);ap.add_argument('--task',required=True,type=int);a=ap.parse_args()
 verify_gate_owner();plan=json.loads(a.manifest.read_text());task=plan['tasks'][a.task]
 for path,h in plan['producer_sha256'].items():assert sha(path)==h,'Producer changed'
 assert len(os.sched_getaffinity(0))==2
 sd=Path(task['source']);raw=json.loads((sd/'result.json').read_text());q,m=raw['charge'],raw['mult']
 for path,h in task['input_sha256'].items():assert sha(path)==h,'Input changed'
 out=a.manifest.parent/task['key'];out.mkdir(exist_ok=False);work=out/'work';work.mkdir()
 for name in ('opt.xyz','result.json'):shutil.copy2(sd/name,work/name)
 (out/'producer_source.py').write_bytes(Path(__file__).read_bytes())
 report=dict(status='running',task=task,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],cpu_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),protocol=N.NATIVE_THERMAL,cds_step_Bohr=1e-6,production_admission=False,raw_inputs_modified=False)
 def save():atomic_json(out/'status.json',report)
 save()
 import gpu4pyscf.solvent.hessian.smd as H
 from gpu4pyscf.solvent.grad import smd as G
 original=H.get_cds
 def fine_cds(obj):
  mol=obj.mol.copy();tmp=obj.copy();coords=mol.atom_coords(unit='B');matrix=np.empty((mol.natm,mol.natm,3,3));step=1e-6
  for atom in range(mol.natm):
   for axis in range(3):
    gradients=[]
    for sign in (-1,1):
     x=coords.copy();x[atom,axis]+=sign*step;mol.set_geom_(x,unit='B');tmp.reset(mol);gradients.append(np.asarray(G.get_cds(tmp)).copy())
    matrix[atom,:,axis]=(gradients[1]-gradients[0])/(2*step)
  return matrix
 H.get_cds=fine_cds
 try:
  baseline=N.load_native_thermal(sd,q,m) if task['control'] else None
  if task['control']:assert baseline is not None
  _,thermal=N.ensure_native_thermal(task['id'],task['state'],state_dir=work,memory_cap_gib=24)
  path,_=N.cache_path(work,q,m);matrix=np.load(path.with_suffix('')/'hessian_cartesian.npy')
  record=json.loads(path.read_text());record.update(status='isolated_CDS_1e6_qualification_pending',execution_derivatives=dict(cds_step_Bohr=1e-6,qualification_only=True,production_admission=False));atomic_json(path,record)
  assert N.load_native_thermal(work,q,m) is None
  report.update(status='checking_worst_contact',native_standard_gates_passed=True,maximum_asymmetry=thermal['maximum_asymmetry'],minimum_frequency_cm=min(thermal['frequencies_cm']),derivative_checks=thermal['derivative_checks'],thermal_eV=thermal['g_thermal_total_eV']);save()
  if task['control']:
   bp,_=N.cache_path(sd,q,m);B=np.load(bp.with_suffix('')/'hessian_cartesian.npy')
   report.update(relative_matrix_error=float(np.linalg.norm(matrix-B)/np.linalg.norm(B)),maximum_matrix_difference=float(abs(matrix-B).max()),maximum_thermal_difference_eV=max(abs(thermal['nuclear_by_cutoff_eV'][str(c)]-baseline['nuclear_by_cutoff_eV'][str(c)]) for c in (25,50,100)))
   assert all(report[k]<=1e-6 for k in ('relative_matrix_error','maximum_matrix_difference','maximum_thermal_difference_eV'))
  from ase.io import read
  from pyscf import gto
  from redox.qm.dft import _build_smd_mf,_ecp_for,_free_gpu_pool
  from gpu4pyscf.solvent.smd import SMD
  atoms=read(work/'opt.xyz');symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
  basis=N.NATIVE_THERMAL['anion_basis'] if q<0 else N.NATIVE_THERMAL['basis']
  def build(x):
   mol=gto.M(atom=list(zip(symbols,x.tolist())),unit='Bohr',charge=q,spin=m-1,basis=basis,ecp=_ecp_for(symbols,basis),verbose=0,max_memory=6000)
   mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
   mf.kernel();assert mf.converged;return mol,mf
  mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=q,spin=m-1,basis=basis,verbose=0)
  coarse=np.asarray(original(SMD(mol,solvent='acetonitrile'))).transpose(0,2,1,3).reshape(3*len(atoms),3*len(atoms));np.save(out/'coarse_CDS_hessian.npy',coarse)
  i,j=np.unravel_index(np.argmax(abs(coarse-coarse.T)),coarse.shape);assert i!=j
  report.update(worst_contact_indices=[int(i),int(j)],maximum_coarse_CDS_asymmetry=float(abs(coarse-coarse.T).max()),worst_contact_checks=[]);save()
  for column in sorted({int(i),int(j)}):
   for step in (.0005,.001):
    gradients=[]
    for sign in (-1,1):
     x=coords.copy();x.flat[column]+=sign*step;_,mf=build(x);grad=mf.nuc_grad_method();grad.grid_response=True;v=grad.kernel();g=np.asarray(v.get() if hasattr(v,'get') else v).ravel().copy();gradients.append(g);np.save(out/f'gradient_{column}_{step}_{sign}.npy',g);del mf,grad,v;_free_gpu_pool()
    numerical=(gradients[1]-gradients[0])/(2*step);error=float(np.linalg.norm(numerical-matrix[column])/np.linalg.norm(numerical));report['worst_contact_checks'].append(dict(column=int(column),step_Bohr=step,relative_error=error));save();assert error<=.01
  assert all(sha(p)==h for p,h in task['input_sha256'].items()),'Source changed during qualification'
  report.update(status='isolated_control_verified' if task['control'] else 'isolated_target_verified',finished=time.time(),evidence_sha256={str(p):sha(p) for p in out.rglob('*') if p.is_file() and p.name not in ('status.json',)})
  save()
 except BaseException as exc:
  # No failure record at a private path may resemble an admitted native record.
  try:
   path,_=N.cache_path(work,q,m)
   if path.exists():
    rec=json.loads(path.read_text());rec['status']='isolated_CDS_1e6_qualification_failed';atomic_json(path,rec)
  except Exception:pass
  report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc(),finished=time.time());save();raise
 finally:H.get_cds=original
if __name__=='__main__':main()
