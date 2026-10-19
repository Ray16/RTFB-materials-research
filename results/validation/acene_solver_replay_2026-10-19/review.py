"""Read four raw energies, units, geometry identity and numerical reproducibility."""
import hashlib,json,math,os,time
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
base=Path(__file__).resolve().parent
with allocation_guard():
 resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
plan=json.loads((base/'plan.json').read_text());launch=json.loads((base/'launch.json').read_text());checks=[];cases=[]
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def topology(mol):
 rw=Chem.RWMol()
 for a in mol.GetAtoms():rw.AddAtom(Chem.Atom(a.GetAtomicNum()))
 for b in mol.GetBonds():rw.AddBond(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),Chem.BondType.SINGLE)
 return rw.GetMol()
for p,h in plan['producer_sha256'].items():check('producer '+p,sha(p)==h)
for task,w in zip(plan['tasks'],launch['workers']):
 row=task['row'];gid=row['id'];path=base/'calc'/(gid+'.json');r=json.loads(path.read_text());proof=r['four_point_evidence'];e=proof['energies_Ha']
 for p,h in task['input_sha256'].items():check(gid+' original bytes '+p,sha(p)==h)
 report=json.loads((base/'tasks'/gid/'status.json').read_text());execution=json.loads((base/f"reserved_{w['task']}"/'execution.json').read_text())
 check(gid+' completed solver execution',r['status']=='ok' and report['status']=='fresh_solver_qualified_record_written' and execution['status']=='process_completed' and execution['returncode']==0)
 check(gid+' result hash',sha(path)==report['result_sha256'])
 check(gid+' finite four-point energies',len(e)==4 and all(type(v) in (float,int) and math.isfinite(v) for v in e.values()))
 contract=r['solver_verification']['contract'];check(gid+' explicit convergence/full-response contract',all(contract[k] is True for k in ['finite_SCF_energy','SCF_converged','optimizer_converged','optimization_SCF_converged','full_grid_response']))
 ident=r['solver_verification']['identity'];check(gid+' exact method and material',ident['basis']=='6-31g(d,p)' and ident['xc']=='b3lyp' and ident['phase']=='gas' and ident['smiles']==row['smiles'])
 ha=27.211386245988
 rc=(e['ion_at_neutral']-e['ion_at_ion'])*ha;rn=(e['neutral_at_ion']-e['neutral_at_neutral'])*ha
 total=math.fsum([e['ion_at_neutral'],e['neutral_at_ion'],-e['ion_at_ion'],-e['neutral_at_neutral']])*ha
 check(gid+' unrounded energy identity',abs(total-proof['lambda_unrounded_eV'])<1e-8)
 check(gid+' serialized total',round(total,4)==r['our_hole'])
 check(gid+' serialized charged relaxation',round(rc,4)==r['relax_cat'])
 check(gid+' serialized neutral relaxation',round(rn,4)==r['relax_neu'])
 check(gid+' nonnegative directional relaxations',rc>=0 and rn>=0)
 expected=Chem.AddHs(Chem.MolFromSmiles(row['smiles']));desired=topology(expected)
 for key in ['neutral_geometry_A','ion_geometry_A']:
  atoms=proof[key];block=str(len(atoms))+'\nreview\n'+'\n'.join(symbol+' '+' '.join(str(v) for v in xyz) for symbol,xyz in atoms)+'\n'
  observed=Chem.MolFromXYZBlock(block);check(gid+' '+key+' composition',observed is not None and sorted(a.GetAtomicNum() for a in observed.GetAtoms())==sorted(a.GetAtomicNum() for a in expected.GetAtoms()))
  rdDetermineBonds.DetermineConnectivity(observed,charge=0);actual=topology(observed)
  check(gid+' '+key+' graph',actual.HasSubstructMatch(desired) and desired.HasSubstructMatch(actual))
 check(gid+' charge/spin convention',proof['ion_charge']==1 and proof['ion_spin']==1)
 old=json.loads(Path(task['legacy_record']).read_text())
 cases.append(dict(id=gid,fresh_lambda_eV=r['our_hole'],fresh_unrounded_lambda_eV=total,legacy_serialized_lambda_eV=old['our_hole'],change_at_saved_precision_eV=r['our_hole']-old['our_hole'],curated_computational_literature_eV=float(row['lit_hole_eV']),difference_from_computational_literature_eV=r['our_hole']-float(row['lit_hole_eV']),literature_is_experimental=False,legacy_record=task['legacy_record'],fresh_record=str(path),fresh_record_sha256=sha(path),true_minimum_Hessian_checked=False,preferred_spin_beyond_stated_model_verified=False))
terminal=json.loads((base/'terminal_verification.json').read_text());check('actual process termination',terminal['status']=='both_reference_workers_and_controllers_verified_terminal')
result=dict(status='two_actual_reference_replays_independently_checked' if all(c['passed'] for c in checks) else 'failed',cases=cases,checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),details=checks,actual_quantum_recomputed_cases=2,legacy_records_modified=False,other_legacy_records_not_retroactively_certified=True,production_promotion=False,experimental_solution_kinetics_validated=False,resources=resources,CPU_affinity=sorted(os.sched_getaffinity(0)),reviewed=time.time())
atomic_write_json(base/'review.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['details','resources']}))
assert not result['failed_checks']
