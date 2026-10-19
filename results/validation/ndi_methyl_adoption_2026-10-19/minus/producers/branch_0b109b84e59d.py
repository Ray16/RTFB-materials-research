"""Bounded isolated signed-mode search, fine-CDS qualification and active SP.

Canonical data and the original failed qualification case are never modified.
No branch is admitted before a separate certificate/source/promotion review.
"""
import argparse,fcntl,gc,hashlib,json,os,runpy,socket,sys,time,traceback
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from redox.core.common import ROOT,read_manifest,smiles_to_xyz_map
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.qm import native_thermal as N,dft as D
from redox.qm.sp import _ENERGY_KEYS,ensure
from redox.workflows.resources import verify_gate_owner
BASE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def host(v):return np.asarray(v.get() if hasattr(v,'get') else v)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--branch',choices=('minus','plus'),required=True);a=ap.parse_args();verify_gate_owner()
 source=json.loads((BASE/'source_plus.json').read_text())
 for path,h in source['inputs_sha256'].items():assert sha(path)==h,'Source changed before mode search'
 for path,h in source['producer_sha256'].items():assert sha(path)==h,'Producer changed'
 work=BASE/a.branch;writer=(work/'writer.lock').open('a+');fcntl.flock(writer,fcntl.LOCK_EX|fcntl.LOCK_NB);assert (work/'seed.xyz').exists() and not (work/'result.json').exists()
 status=dict(status='optimizing',branch=a.branch,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],cpu_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),source=source,production_admission=False,canonical_modified=False)
 def save(name,value):atomic_write_json(work/name,value)
 def progress():save('status.json',status)
 progress();atoms=read(work/'seed.xyz');symbols=atoms.get_chemical_symbols();raw=json.loads((BASE/'source_archive/result.json').read_text());q,m=source['charge'],source['mult']
 try:
  from pyscf.geomopt import geometric_solver as GS
  import geometric.optimize as GO
  mol=gto.M(atom=list(zip(symbols,atoms.positions.tolist())),charge=q,spin=m-1,basis=N.NATIVE_THERMAL['anion_basis'] if q<0 else N.NATIVE_THERMAL['basis'],verbose=0,max_memory=6000)
  mf=D._build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.kernel();assert mf.converged
  grad=mf.nuc_grad_method();grad.grid_response=True;capture={};cycles=[];original=GO.Optimizer
  class CapturedOptimizer(original):
   def optimizeGeometry(self):
    try:return super().optimizeGeometry()
    finally:
     rms,maxg=self.calcGradNorm();drms,dmax=GO.calc_drms_dmax(self.X,self.Xprev);params=self.params
     values=dict(energy_change_Ha=float(abs(self.E-self.Eprev)),gradient_rms_Ha_Bohr=float(rms),gradient_max_Ha_Bohr=float(maxg),displacement_rms_A=float(drms),displacement_max_A=float(dmax))
     limits=dict(energy_change_Ha=float(params.Convergence_energy),gradient_rms_Ha_Bohr=float(params.Convergence_grms),gradient_max_Ha_Bohr=float(params.Convergence_gmax),displacement_rms_A=float(params.Convergence_drms),displacement_max_A=float(params.Convergence_dmax))
     capture.update(state=int(self.state),CONVERGED_state=int(GO.OPT_STATE.CONVERGED),iterations=int(self.Iteration),values=values,limits=limits,all_five_criteria_passed=all(values[k]<limits[k] for k in values),final_coordinates_Bohr=self.X.reshape(-1,3).tolist());save('geometric_convergence.json',capture)
  GO.Optimizer=CapturedOptimizer
  def callback(info):
   g=host(info['gradients']);cycles.append(dict(cycle=len(cycles)+1,e_Ha=float(info['energy']),gradient_max_atomic_norm_Ha_Bohr=float(np.linalg.norm(g,axis=1).max()),SCF_converged=bool(info['g_scanner'].converged)));save('optimization_progress.json',dict(cycles=cycles,last=cycles[-1]))
  try:conv,final=GS.kernel(grad.as_scanner(),maxsteps=120,assert_convergence=True,callback=callback,prefix=str(work/'geometric'),convergence_energy=1e-9,convergence_grms=1e-6,convergence_gmax=3e-6,convergence_drms=3e-5,convergence_dmax=5e-5)
  finally:GO.Optimizer=original
  marker=('Converged! =D' in (work/'geometric.log').read_text()) if (work/'geometric.log').exists() else None;aligned=np.max(abs(final.atom_coords(unit='Bohr')-np.asarray(capture['final_coordinates_Bohr'])))<1e-6
  capture.update(pyscf_kernel_converged=bool(conv),geometric_log_convergence_marker=marker,final_geometry_matches_converged_state=bool(aligned));save('geometric_convergence.json',capture)
  if not(conv and aligned and capture['state']==capture['CONVERGED_state'] and capture['all_five_criteria_passed']):raise RuntimeError('Explicit optimizer convergence gate failed')
  result_atoms=atoms.copy();result_atoms.positions=final.atom_coords(unit='Angstrom');write(work/'opt.xyz',result_atoms)
  if smiles_to_xyz_map(source['smiles'],work/'opt.xyz') is None:raise RuntimeError('Mode-follow branch changes exact graph')
  trial=dict(raw,optimized=True,opt_out=str(work/'opt.xyz'),g_thermal_eV=None,n_imag=None,geometry_optimization_explicitly_converged=True,signed_mode_source=str(BASE/'source_plus.json'))
  for k in _ENERGY_KEYS:trial[k]=None
  save('result.json',trial);del mf,grad,final;gc.collect();D._free_gpu_pool()
  status['status']='qualifying_finer_native_derivatives';progress()
  directory=work/'qualification';directory.mkdir();qualifier=BASE/'qualify.py'
  plan=dict(protocol=N.NATIVE_THERMAL,producer_sha256=source['producer_sha256'],tasks=[dict(key='candidate',id=source['id'],state=source['state'],source=str(work),charge=q,mult=m,control=False,input_sha256={str(work/k):sha(work/k) for k in ('opt.xyz','result.json')})])
  atomic_write_json(directory/'plan.json',plan);sys.argv=[str(qualifier),'--manifest',str(directory/'plan.json'),'--task','0'];runpy.run_path(str(qualifier),run_name='__main__')
  qualification=json.loads((directory/'candidate/status.json').read_text());assert qualification['status']=='isolated_target_verified'
  status['status']='computing_active_energy';progress();kind,sp=ensure(source['id'],source['state'],work,q,m,backend='gpu');assert kind=='computed' and load_record(work,q,m) is not None
  for path,h in source['inputs_sha256'].items():assert sha(path)==h,'Canonical source changed during isolated search'
  status.update(status='private_minimum_and_energy_verified_pending_certificate',active_SP=str(record_path(work,q,m)),active_SP_sha256=sha(record_path(work,q,m)),e_smd_eV=sp['e_smd_eV'],qualification_status=str(directory/'candidate/status.json'),minimum_frequency_cm=qualification['minimum_frequency_cm'],finished=time.time(),production_admission=False);progress()
 except BaseException as exc:status.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc(),finished=time.time());progress();raise
if __name__=='__main__':main()
