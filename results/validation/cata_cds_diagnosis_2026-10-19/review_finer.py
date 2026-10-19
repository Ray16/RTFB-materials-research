"""Independently reconstruct spectrum and all saved worst-contact finite differences."""
import hashlib,json,math,os,time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import NATIVE_THERMAL
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
base=Path(__file__).resolve().parent;replay=base/'finer_replay';case=replay/'cata_neutral'
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
plan=json.loads((replay/'plan.json').read_text());task=plan['tasks'][0];report=json.loads((case/'status.json').read_text());checks=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
check('qualified isolated target',report['status']=='isolated_target_verified' and report['production_admission'] is False)
for p,h in {**plan['producer_sha256'],**task['input_sha256'],**report['evidence_sha256']}.items():check('input/output source '+p,sha(p)==h)
work=case/'work';path=next((work/'native_thermal').glob('*.json'));r=json.loads(path.read_text());matrix=np.load(path.with_suffix('')/'hessian_cartesian.npy');atoms=read(work/'opt.xyz');nat=len(atoms)
check('isolated record remains unadmitted',r['status']=='isolated_CDS_1e6_qualification_pending' and r['execution_derivatives']['qualification_only'] is True)
check('geometry and protocol',r['identity']['geometry_sha256']==sha(work/'opt.xyz')==sha(Path(task['source'])/'opt.xyz') and r['identity']['protocol']==NATIVE_THERMAL)
check('finite symmetric matrix',matrix.shape==(nat*3,nat*3) and np.isfinite(matrix).all() and abs(matrix-matrix.T).max()<=1e-6)
check('SCF and stationary gradient',r['stationary_minimum'] and r['converged'] and r['grid_response'] and math.isfinite(r['gradient_max']) and 0<=r['gradient_max']<=.00045)
check('stated spin diagnostics',abs(r['smd_spin_diagnostics']['spin_contam'])<=.1)
mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),(atoms.positions/.52917721067).tolist())),unit='Bohr',basis='def2-svp',charge=0,spin=0,verbose=0)
freq=thermo.harmonic_analysis(mol,matrix.reshape(nat,3,nat,3).transpose(0,2,1,3),imaginary_freq=False)['freq_wavenumber']
check('independent vibrational spectrum',len(freq)==3*nat-6 and np.isfinite(freq).all() and min(freq)>=-5 and np.max(abs(freq-np.asarray(r['frequencies_cm'])))<1e-7)
for item in r['derivative_checks']:check('standard derivative '+str(item),math.isfinite(item['relative_error']) and 0<=item['relative_error']<=.01)
coarse=np.load(case/'coarse_CDS_hessian.npy');i,j=np.unravel_index(np.argmax(abs(coarse-coarse.T)),coarse.shape)
check('worst-contact identity',[int(i),int(j)]==report['worst_contact_indices'])
errors=[]
for column in sorted({int(i),int(j)}):
 for step in (.0005,.001):
  plus=np.load(case/f'gradient_{column}_{step}_1.npy');minus=np.load(case/f'gradient_{column}_{step}_-1.npy');numerical=(plus-minus)/(2*step);error=float(np.linalg.norm(numerical-matrix[column])/np.linalg.norm(numerical));recorded=next(c['relative_error'] for c in report['worst_contact_checks'] if c['column']==column and c['step_Bohr']==step)
  check(f'independent full-gradient FD {column} {step}',np.isfinite(plus).all() and np.isfinite(minus).all() and error<=.01 and abs(error-recorded)<1e-12);errors.append(error)
old=np.load(base/'frozen_cata/full_hessian.npy');cpucoarse=np.load(base/'cpu_cds_diagnosis/cata_coarse_CDS.npy');cpufine=np.load(base/'cpu_cds_diagnosis/cata_fine_CDS.npy');estimate=old-cpucoarse+cpufine
relative=float(np.linalg.norm(matrix-estimate)/np.linalg.norm(matrix));check('fresh quantum matrix vs independent component reconstruction',relative<=1e-6)
evidence={str(case/'status.json'):sha(case/'status.json'),**report['evidence_sha256'],**task['input_sha256']}
result=dict(status='finer_target_independently_verified' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),details=checks,minimum_frequency_cm=float(min(freq)),maximum_asymmetry=float(abs(matrix-matrix.T).max()),maximum_contact_relative_error=max(errors),relative_matrix_difference_from_component_reconstruction=relative,evidence_sha256=evidence,canonical_modified=False,production_admission=False,reviewed=time.time())
atomic_write_json(base/'finer_independent_review.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ('details','evidence_sha256')}))
assert not result['failed_checks']
