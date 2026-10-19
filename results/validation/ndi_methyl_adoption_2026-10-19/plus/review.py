"""Independent physical/energy review before promoting the preserved dianion."""
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():
    resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,smiles_to_xyz_map
from redox.core.protocol import ACTIVE_SP,load_record
from redox.qm.native_thermal import load_native_thermal
from redox.core.io import atomic_write_json

BASE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    receipt=json.loads((BASE/'preparation.json').read_text());stage=Path(receipt['prepared'])
    case=Path(receipt['qualified_case']);qual=json.loads((case/'status.json').read_text())
    native=json.loads(Path(receipt['native_record']).read_text());q,m=receipt['charge'],receipt['mult']
    checks=[];evidence={str(BASE/'preparation.json'):sha(BASE/'preparation.json')}
    def check(name,value):checks.append(dict(check=name,passed=bool(value)))
    for path,digest in receipt['evidence_sha256'].items():
        check('frozen:'+path,sha(path)==digest);evidence[path]=digest
    atoms=read(stage/'opt.xyz');nat=len(atoms)
    matrix=np.load(Path(receipt['native_record']).with_suffix('')/'hessian_cartesian.npy')
    check('native_matrix',matrix.shape==(3*nat,3*nat) and np.isfinite(matrix).all()
          and abs(matrix-matrix.T).max()<=1e-6)
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),(atoms.positions/.52917721067).tolist())),
              unit='Bohr',basis='def2-svpd',charge=q,spin=m-1,verbose=0)
    freq=np.asarray(thermo.harmonic_analysis(mol,matrix.reshape(nat,3,nat,3).transpose(0,2,1,3),
                                          imaginary_freq=False)['freq_wavenumber'])
    check('native_minimum',len(freq)==3*nat-6 and min(freq)>=-5)
    check('native_spectrum_replay',max(abs(freq-native['frequencies_cm']))<1e-5)
    check('exact_graph',smiles_to_xyz_map(receipt['smiles'],stage/'opt.xyz') is not None)
    geometry=json.loads((Path(receipt['branch'])/'geometric_convergence.json').read_text())
    check('optimizer_termination',geometry['state']==geometry['CONVERGED_state']
          and geometry['pyscf_kernel_converged'] is True)
    check('all_five_optimizer_criteria',set(geometry['values'])==set(geometry['limits'])
          and all(math.isfinite(v) and abs(v)<geometry['limits'][k] for k,v in geometry['values'].items()))
    check('optimizer_final_geometry',np.max(abs(atoms.positions/.52917721067-
          np.asarray(geometry['final_coordinates_Bohr'])))<1e-6)
    estimates=[]
    for column in sorted(set(qual['worst_contact_indices'])):
        midpoint={}
        for step in (.0005,.001):
            plus=np.load(case/f'gradient_{column}_{step}_1.npy')
            minus=np.load(case/f'gradient_{column}_{step}_-1.npy')
            derivative=(plus-minus)/(2*step)
            error=float(np.linalg.norm(derivative-matrix[column])/np.linalg.norm(derivative))
            saved=next(v['relative_error'] for v in qual['worst_contact_checks']
                       if v['column']==column and v['step_Bohr']==step)
            check(f'worst_contact:{column}:{step}',np.isfinite(plus).all() and np.isfinite(minus).all()
                  and error<=.01 and abs(error-saved)<1e-12)
            midpoint[step]=(plus+minus)/2
        estimates.append((4*midpoint[.0005]-midpoint[.001])/3)
    spread=float(np.max(np.ptp(np.asarray(estimates),axis=0)))
    estimated_max=float(np.max(abs(np.asarray(estimates))))
    check('independent_midpoint_stationarity_diagnostic',estimated_max+spread<=.00045)
    check('producer_stationarity',math.isfinite(native['gradient_max']) and 0<=native['gradient_max']<=.00045)
    energy=load_record(stage,q,m)
    check('addressed_energy',energy is not None and energy['protocol']==ACTIVE_SP)
    if energy:
        for phase in ('smd','gas'):
            check('energy_units:'+phase,energy['converged_'+phase] is True
                  and abs(energy['e_'+phase+'_eV']-energy['e_'+phase+'_Ha']*27.211386245988)<1e-8
                  and math.isfinite(energy[phase+'_spin_contam']) and abs(energy[phase+'_spin_contam'])<=.1)
    spec=importlib.util.spec_from_file_location('independent_adopted_state',ROOT/'scripts/validation/audit/recompute_axes.py')
    audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
    raw=json.loads((stage/'result.json').read_text());correction,imaginary,definition=audit.thermal_inputs(stage,raw)
    primary=load_native_thermal(stage,q,m)
    check('both_primary_readers',primary is not None and correction is not None and imaginary==0
          and definition=='native_solution_dft_qrrho_spin'
          and abs(correction-primary['g_thermal_total_eV'])<1e-8)
    for path in stage.rglob('*'):
        if path.is_file() and path.suffix!='.lock':evidence[str(path)]=sha(path)
    for path,digest in evidence.items():check('final_hash:'+path,sha(path)==digest)
    report=dict(status='independent_review_passed' if all(c['passed'] for c in checks) else 'failed',
                checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),evidence_sha256=evidence,
                minimum_frequency_cm=float(min(freq)),midpoint_stationarity_diagnostic_max_Ha_Bohr=estimated_max,
                midpoint_spread_Ha_Bohr=spread,corrected_thermal_eV=correction,
                source_optimization_path_independently_checked=True,resources=resources,
                source_sha256=sha(Path(__file__)),canonical_modified=False,finished=time.time())
    atomic_write_json(BASE/'independent_review.json',report);atomic_write_json(BASE/'checks.json',checks)
    print({k:v for k,v in report.items() if k not in ('evidence_sha256','resources')})
    assert report['failed_checks']==0


if __name__=='__main__':main()
