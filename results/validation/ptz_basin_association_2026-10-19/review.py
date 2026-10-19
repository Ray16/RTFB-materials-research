"""Independently assess reciprocal return to selected endpoint geometry/energy."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from ase.io import read
from rdkit import Chem
from redox.core.common import ROOT,smiles_to_xyz_map
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());rows=[];checks=[];evidence={}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rmsd(a,b,reflection=False):
 a=a-a.mean(axis=0);b=b-b.mean(axis=0);u,_,vt=np.linalg.svd(a.T@b);sign=np.linalg.det(u@vt);rotation=u@np.diag([1.,1.,-sign if reflection else sign])@vt
 return float(np.sqrt(np.mean(np.sum((a@rotation-b)**2,axis=1))))
for task in plan['tasks']:
 work=base/task['target'];status=json.loads((work/'status.json').read_text())
 if status['status']!='reciprocal_optimization_verified_pending_basin_review':
  rows.append(dict(target=task['target'],status=status['status'],basin_associated=False));continue
 selected=read(task['selected_geometry']);trial=read(work/'opt.xyz');smiles=task['smiles'];mol=Chem.MolFromSmiles(smiles)
 match_a=smiles_to_xyz_map(smiles,task['selected_geometry']);match_b=smiles_to_xyz_map(smiles,work/'opt.xyz');assert match_a and match_b
 xa=selected.positions[[a.symbol!='H' for a in selected]][list(match_a)];xb=trial.positions[[a.symbol!='H' for a in trial]][list(match_b)]
 permutations=mol.GetSubstructMatches(mol,uniquify=False,useChirality=True,maxMatches=10000)
 proper=min(rmsd(xa,xb[list(p)]) for p in permutations);mirror=min(rmsd(xa,xb[list(p)],True) for p in permutations)
 convergence=json.loads((work/'geometric_convergence.json').read_text());g=np.load(work/'final_gradient.npy')
 local_checks=dict(explicit_optimizer_termination=convergence['state']==convergence['CONVERGED_state'] and convergence['pyscf_kernel_converged'],all_five_criteria=all(abs(v)<convergence['limits'][k] for k,v in convergence['values'].items()),final_gradient=float(abs(g).max())<=4.5e-4,local_energy_match=abs(status['delta_from_selected_eV'])<=.001,proper_rotation_heavy_geometry_match=proper<=.05)
 for k,v in local_checks.items():checks.append(dict(target=task['target'],check=k,passed=bool(v)))
 for path,h in task['input_sha256'].items():assert sha(path)==h;evidence[path]=h
 for p in [work/'status.json',work/'opt.xyz',work/'geometric_convergence.json',work/'final_gradient.npy',work/'selected_endpoint_SCF.json',work/'reciprocal_final_SCF.json']:evidence[str(p)]=sha(p)
 rows.append(dict(target=task['target'],status=status['status'],proper_symmetry_RMSD_A=proper,reflection_symmetry_RMSD_A=mirror,energy_difference_eV=status['delta_from_selected_eV'],graph_automorphisms=len(permutations),basin_associated=all(local_checks.values())))
complete=all(r['status']=='reciprocal_optimization_verified_pending_basin_review' for r in rows)
report=dict(status='reciprocal_local_branch_review' if complete else 'waiting_for_reciprocal_endpoint',endpoints=rows,checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),both_selected_local_basins_associated=complete and all(r['basin_associated'] for r in rows),criteria=dict(energy_eV=.001,heavy_RMSD_A=.05,gradient_Ha_Bohr=.00045),global_conformer_ensemble_proven=False,active_SP_stationary_endpoints_proven=False,empirical_kinetics_validated=False,production_records_modified=False,evidence_sha256=evidence,source_sha256=sha(Path(__file__)),reviewed=time.time())
atomic_write_json(base/'review.json',report);atomic_write_json(base/'checks.json',checks);print({k:v for k,v in report.items() if k!='evidence_sha256'})
