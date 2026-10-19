#!/usr/bin/env python
"""Test full-grid-response analytical RKS Hessian against saved conservative derivatives.

Prepared validation only; a numerical reference is required before this can pass.
Run through the shared GPU reservation gate with two pinned CPU threads.
"""
import argparse,hashlib,json,os,socket,time,subprocess
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,DFT
from redox.qm.dft import _build_smd_mf,_free_gpu_pool
from redox.core.protocol import software_versions
from redox.qm.vibrations import write_xtb_hessian,parse_xtb_thermal_eV
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    ref_path=BASE/'numerical_solution_hessian_gridresponse_verified/result.json'
    if not ref_path.exists():raise RuntimeError('Conservative numerical reference not complete')
    reference=json.loads(ref_path.read_text())
    if reference['status']!='numerical_reference_verified':raise RuntimeError('Numerical reference not verified')
    if a.out.exists():raise ValueError('Preserve prior analytical experiments; use a new output')
    if len(os.sched_getaffinity(0))!=2:raise RuntimeError('Two pinned CPU threads required')
    xyz=DFT/'aq_benzylamino/red2/opt.xyz'
    if hashlib.sha256(xyz.read_bytes()).hexdigest()!=reference['identity']['geometry_sha256']:raise RuntimeError('Reference geometry stale')
    a.out.mkdir(parents=True);atoms=read(xyz)
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
        cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),source_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
        software=software_versions(),grid_level=3,grid_response=True,production_modified=False)
    def save():
        p=a.out/'status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    save()
    try:
        # Retain the demonstrated safe DF batching control from the numerical-reference
        # campaign's analytical diagnostic; only tiling changes, not the Hessian terms.
        import gpu4pyscf.df.hessian.rhf as dfh
        original_avail=dfh.get_avail_mem
        dfh.get_avail_mem=lambda *args,**kwargs:min(original_avail(*args,**kwargs),6*1024**3)
        status['df_batch_memory_limit_GiB']=6;save()
        mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),charge=-2,spin=0,basis='def2-svpd',max_memory=6000,verbose=4)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.kernel()
        if not mf.converged:raise RuntimeError('Analytical Hessian source SCF failed')
        gradient=mf.nuc_grad_method();gradient.grid_response=True
        raw=gradient.kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw)
        np.save(a.out/'gradient_full_response.npy',g)
        hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
        raw=hobj.kernel();H=np.asarray(raw.get() if hasattr(raw,'get') else raw)
        np.save(a.out/'hessian.npy',H)
        mat=H.transpose(0,2,1,3).reshape(3*len(atoms),3*len(atoms))
        exact=np.load(BASE/'numerical_solution_hessian_gridresponse_verified/hessian_matrix.npy')
        diff=mat-exact
        columns=[4,50,70,73]
        checks=[dict(column=i,relative_norm_error=float(np.linalg.norm(diff[:,i])/np.linalg.norm(exact[:,i]))) for i in columns]
        spectrum=thermo.harmonic_analysis(mol,H,exclude_trans=True,exclude_rot=True,imaginary_freq=False)
        freq=np.asarray(spectrum['freq_wavenumber'])
        status.update(max_abs_error_Ha_Bohr2=float(abs(diff).max()),relative_full_matrix_error=float(np.linalg.norm(diff)/np.linalg.norm(exact)),
            maximum_asymmetry_Ha_Bohr2=float(abs(mat-mat.T).max()),column_checks=checks,
            gradient_max=float(abs(g).max()),frequencies_cm=freq.tolist(),minimum_frequency_cm=float(freq.min()))
        write_xtb_hessian(a.out/'hessian',(mat+mat.T)/2);(a.out/'mol.xyz').write_bytes(xyz.read_bytes())
        thermal={};delta={}
        for cutoff in (25,50,100):
            log=a.out/f'thermo_{cutoff}.log'
            with log.open('w') as f:
                p=subprocess.run(['/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb','thermo',
                    '--sthr',str(cutoff),'--temp','298.15','mol.xyz','hessian'],cwd=a.out,stdout=f,stderr=subprocess.STDOUT,timeout=120)
            if p.returncode:raise RuntimeError('Analytical Hessian thermochemistry processing failed')
            thermal[str(cutoff)]=parse_xtb_thermal_eV(log.read_text())
            delta[str(cutoff)]=thermal[str(cutoff)]-reference['g_thermal_qrrho_by_cutoff_eV'][str(cutoff)]
        status.update(thermal_by_cutoff_eV=thermal,thermal_error_from_numerical_eV=delta)
        passed=(status['relative_full_matrix_error']<.01 and all(r['relative_norm_error']<.01 for r in checks)
            and status['maximum_asymmetry_Ha_Bohr2']<1e-6 and freq.min()>=-5
            and max(abs(x) for x in delta.values())<.005)
        status.update(status='analytic_full_response_verified' if passed else 'analytic_full_response_failed',finished=time.time());save()
        _free_gpu_pool()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
