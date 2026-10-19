"""Lower-level-curvature sensitivity estimate only; no geometry or energy adoption."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from ase.io import read
from redox.core.common import DFT
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import cache_path
from redox.qm.vibrations import project_rigid_modes
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
base=Path(__file__).resolve().parent
with allocation_guard():
 resources=host_capacity(1);cpus=cpu_allocation(1);os.sched_setaffinity(0,cpus)
plan=json.loads((base/'plan.json').read_text());result=[]
for task in plan['tasks']:
 state=task['state'];geometry=Path(task['source'])/'geometry.xyz';atoms=read(geometry)
 directory=DFT/'bm_n_ethyl_3_7_dimethoxyphenothiazine'/state
 assert hashlib.sha256((directory/'opt.xyz').read_bytes()).hexdigest()==task['input_sha256'][str(geometry)]
 native,_=cache_path(directory,task['charge'],task['mult']);hessian=native.with_suffix('')/'hessian_cartesian.npy'
 H=np.load(hessian);g=np.load(base/state/'gradient_Ha_Bohr.npy').ravel();xyz=atoms.positions/.52917721067
 projector,rank=project_rigid_modes(np.eye(g.size),xyz,np.ones(len(atoms)))
 vals,vecs=np.linalg.eigh(projector);B=vecs[:,vals>.5]
 internal=B.T@((H+H.T)/2)@B;forces=B.T@g
 assert np.linalg.eigvalsh(internal).min()>0
 delta=-B@np.linalg.solve(internal,forces);max_step=float(np.linalg.norm(delta.reshape(-1,3),axis=1).max())
 scale=min(1.,(.1/.52917721067)/max_step)
 local=scale*delta;estimated=float((g@local+.5*local@H@local)*27.211386245988)
 result.append(dict(state=state,rigid_rank=rank,unscaled_max_step_A=max_step*.52917721067,
  local_trust_radius_A=.1,step_scale=scale,estimated_electronic_change_eV=estimated,
  curvature_model='r2SCAN-D4/SVP(D) native Hessian; active gradient; unit-mass Euclidean rigid exclusion',
  evidence_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (geometry,hessian,base/state/'gradient_Ha_Bohr.npy')}))
changes={r['state']:r['estimated_electronic_change_eV'] for r in result}
report=dict(status='local_quadratic_diagnostic_completed',states=result,
 estimated_charging_gap_change_eV=changes['ox2']+changes['neu']-2*changes['ox1'],
 actual_active_minimum_found=False,energy_correction_applied=False,production_admission=False,
 caveat='Different-level curvature and a local trust step are not an active-level minimization or a validated free-energy correction',
 CPU_affinity=cpus,resources=resources,completed=time.time())
atomic_write_json(base/'quadratic_sensitivity.json',report);print(json.dumps(report,indent=2))
