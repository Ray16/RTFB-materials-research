"""Reciprocal charge-state optimization control; never overwrite production endpoints."""
import argparse,hashlib,json,os,socket,time,traceback,gc
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from redox.core.common import ROOT,smiles_to_xyz_map,read_manifest
from redox.core.io import atomic_write_json
from redox.qm.dft import _build_smd_mf,_free_gpu_pool
from redox.workflows.resources import verify_gate_owner
base=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def cpu(x):return np.asarray(x.get() if hasattr(x,'get') else x)
def main():
 p=argparse.ArgumentParser();p.add_argument('--target',required=True,choices=['neu','ox1']);arg=p.parse_args();verify_gate_owner()
 plan=json.loads((base/'plan.json').read_text());task=next(t for t in plan['tasks'] if t['target']==arg.target)
 for path,h in plan['source_sha256'].items():assert sha(path)==h,path
 for path,h in task['input_sha256'].items():assert sha(path)==h,path
 work=base/arg.target;work.mkdir(exist_ok=False);snapshots=work/'producers';snapshots.mkdir()
 for path,h in plan['source_sha256'].items():(snapshots/(Path(path).stem+'_'+h[:12]+'.py')).write_bytes(Path(path).read_bytes())
 status=dict(status='baseline_SCF',target=arg.target,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],cpu_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),production_records_modified=False,task=task)
 def save():atomic_write_json(work/'status.json',status)
 save();q,m=task['charge'],task['mult'];smiles=task['smiles'];reference=read(task['selected_geometry']);seed=read(task['opposite_geometry']);assert reference.get_chemical_symbols()==seed.get_chemical_symbols()
 def molecule(atoms):return gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='def2-svp',charge=q,spin=m-1,verbose=0,max_memory=6000)
 def build(mol):
  mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;return mf
 def lowest(atoms,label):
  mol=molecule(atoms);best=None;records=[]
  for guess in (['minao','atom','huckel'] if m>1 else ['minao']):
   mf=build(mol);mf.init_guess=guess;e=float(mf.kernel());ok=bool(mf.converged);contamination=None
   if ok and m>1:contamination=float(mf.spin_square()[0])-(m-1)/2*((m-1)/2+1)
   records.append(dict(guess=guess,e_Ha=e,converged=ok,spin_contamination=contamination))
   if ok and (contamination is None or abs(contamination)<=.1) and (best is None or e<best[0]):best=e,cpu(mf.make_rdm1()).copy()
   del mf;_free_gpu_pool()
  if best is None:raise RuntimeError('No qualified SCF branch')
  atomic_write_json(work/(label+'_SCF.json'),records)
  import cupy
  mf=build(mol);restored=float(mf.kernel(dm0=cupy.asarray(best[1])))
  if not mf.converged or abs(restored-best[0])>1e-8:raise RuntimeError('SCF restoration changed selected branch')
  return mf,restored
 try:
  baseline,e_reference=lowest(reference,'selected_endpoint');del baseline;_free_gpu_pool();mf,e_start=lowest(seed,'reciprocal_start')
  status.update(status='optimizing',selected_endpoint_e_Ha=e_reference,starting_e_Ha=e_start);save()
  from pyscf.geomopt import geometric_solver as GS
  import geometric.optimize as GO
  original=GO.Optimizer;capture={};cycles=[]
  class Captured(original):
   def optimizeGeometry(self):
    try:return super().optimizeGeometry()
    finally:
     grms,gmax=self.calcGradNorm();drms,dmax=GO.calc_drms_dmax(self.X,self.Xprev)
     values=dict(energy_change_Ha=float(abs(self.E-self.Eprev)),gradient_rms_Ha_Bohr=float(grms),gradient_max_Ha_Bohr=float(gmax),displacement_rms_A=float(drms),displacement_max_A=float(dmax))
     limits=dict(energy_change_Ha=float(self.params.Convergence_energy),gradient_rms_Ha_Bohr=float(self.params.Convergence_grms),gradient_max_Ha_Bohr=float(self.params.Convergence_gmax),displacement_rms_A=float(self.params.Convergence_drms),displacement_max_A=float(self.params.Convergence_dmax))
     capture.update(state=int(self.state),CONVERGED_state=int(GO.OPT_STATE.CONVERGED),iterations=int(self.Iteration),values=values,limits=limits,all_five_criteria_passed=all(values[k]<limits[k] for k in values),final_coordinates_Bohr=self.X.reshape(-1,3).tolist());atomic_write_json(work/'geometric_convergence.json',capture)
  def callback(info):
   cycles.append(dict(cycle=len(cycles)+1,e_Ha=float(info['energy']),gradient_max_atomic_norm_Ha_Bohr=float(np.linalg.norm(cpu(info['gradients']),axis=1).max()),SCF_converged=bool(info['g_scanner'].converged)));atomic_write_json(work/'optimization_progress.json',dict(last=cycles[-1],cycles=cycles))
  grad=mf.nuc_grad_method();grad.grid_response=True;GO.Optimizer=Captured
  try:conv,final=GS.kernel(grad.as_scanner(),maxsteps=120,assert_convergence=True,callback=callback,prefix=str(work/'geometric'),convergence_energy=1e-9,convergence_grms=1e-6,convergence_gmax=3e-6,convergence_drms=3e-5,convergence_dmax=5e-5)
  finally:GO.Optimizer=original
  final_atoms=reference.copy();final_atoms.positions=final.atom_coords(unit='Angstrom');write(work/'opt.xyz',final_atoms)
  aligned=float(np.max(abs(final.atom_coords(unit='Bohr')-np.asarray(capture['final_coordinates_Bohr']))))<1e-6
  capture.update(pyscf_kernel_converged=bool(conv),final_geometry_matches_converged_state=aligned);atomic_write_json(work/'geometric_convergence.json',capture)
  if not(conv and aligned and capture['state']==capture['CONVERGED_state'] and capture['all_five_criteria_passed']):raise RuntimeError('Explicit optimizer gate failed')
  if smiles_to_xyz_map(smiles,work/'opt.xyz') is None:raise RuntimeError('Reciprocal path changes exact graph')
  del mf,grad,final;gc.collect();_free_gpu_pool()
  final_mf,e_final=lowest(final_atoms,'reciprocal_final');final_grad=final_mf.nuc_grad_method();final_grad.grid_response=True;g=cpu(final_grad.kernel());np.save(work/'final_gradient.npy',g)
  if not final_mf.converged:raise RuntimeError('Final gradient SCF failed')
  for path,h in task['input_sha256'].items():assert sha(path)==h,'Selected endpoints changed'
  status.update(status='reciprocal_optimization_verified_pending_basin_review',final_e_Ha=e_final,delta_from_selected_eV=(e_final-e_reference)*27.211386245988,final_gradient_max_Ha_Bohr=float(abs(g).max()),finished=time.time(),global_conformer_ensemble_proven=False);save()
 except BaseException as exc:status.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc(),finished=time.time());save();raise
if __name__=='__main__':main()
