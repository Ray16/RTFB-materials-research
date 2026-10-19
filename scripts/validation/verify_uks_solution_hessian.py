#!/usr/bin/env python
"""Independent full numerical-vs-analytical UKS/SMD Hessian control for BQ radical."""
import argparse,hashlib,json,os,socket,subprocess,time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,DFT,read_result
from redox.qm.dft import _build_smd_mf,_free_gpu_pool,_orbital_diag
from redox.qm.vibrations import write_xtb_hessian,parse_xtb_thermal_eV
from redox.core.protocol import software_versions


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--id',default='bm_1_4_benzoquinone');ap.add_argument('--state',default='red1');a=ap.parse_args()
    if a.out.exists():raise ValueError('use a fresh output directory')
    if len(os.sched_getaffinity(0))!=2:raise RuntimeError('two pinned CPU threads required')
    gid=a.id;state=a.state;xyz=DFT/gid/state/'opt.xyz';raw=read_result(gid,state,raw=True)
    if raw is None or raw['mult']!=2:raise RuntimeError('control must be a declared doublet')
    q=raw['charge']
    a.out.mkdir(parents=True);atoms=read(xyz);symbols=atoms.get_chemical_symbols();coordinates=atoms.positions/.52917721067
    basis=raw['opt_basis'];n=3*len(atoms)
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
        process_start=Path('/proc/self/stat').read_text().split()[21],cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),
        id=gid,state=state,charge=q,mult=2,geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
        software=software_versions(),xc='r2scan',basis=basis,disp='d4',solvent='SMD(acetonitrile)',
        density_fit=True,grid_level=3,grid_response=True,conv_tol=1e-12,conv_tol_grad=1e-8,
        step_Bohr=.0005,expected_columns=n,columns=[],production_modified=False)
    def save():
        p=a.out/'status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    save()
    import gpu4pyscf.df.hessian.rhf as df_rhf
    import gpu4pyscf.df.hessian.uhf as df_uhf
    for module in (df_rhf,df_uhf):
        if hasattr(module,'get_avail_mem'):
            old=module.get_avail_mem
            module.get_avail_mem=lambda *args,_old=old,**kwargs:min(_old(*args,**kwargs),6*1024**3)
    def cpu(value):return np.asarray(value.get() if hasattr(value,'get') else value)
    baseline=None
    def solve(coords,use_warm=True):
        mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=q,spin=1,basis=basis,verbose=0,max_memory=6000)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.grids.level=3;mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
        mf.kernel(dm0=baseline if use_warm and baseline is not None else None)
        if not mf.converged:raise RuntimeError('control SCF did not converge')
        grad=mf.nuc_grad_method();grad.grid_response=True
        g=cpu(grad.kernel());return mol,mf,g
    try:
        mol,mf,g0=solve(coordinates)
        status['smd_spin_diagnostics']=_orbital_diag(mf,1)
        if abs(status['smd_spin_diagnostics'].get('spin_contam',99))>.1:raise RuntimeError('UKS control spin contaminated')
        baseline=cpu(mf.make_rdm1());status['gradient_max']=float(abs(g0).max())
        hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
        HA=cpu(hobj.kernel());np.save(a.out/'hessian_analytical_tensor.npy',HA)
        A=HA.transpose(0,2,1,3).reshape(n,n);np.save(a.out/'hessian_analytical_matrix.npy',A)
        del mf,hobj;_free_gpu_pool();save()
        columns=[]
        for k in range(n):
            ia,axis=divmod(k,3);values=[]
            for sign in (-1,1):
                c=coordinates.copy();c[ia,axis]+=sign*.0005
                _,mf,g=solve(c);values.append(g);del mf;_free_gpu_pool()
            column=(values[1]-values[0])/.001;np.save(a.out/f'column_{k:03d}.npy',column)
            columns.append(column.ravel());status['columns'].append(k);save()
            print('UKS_COLUMN_COMPLETE',k+1,n,flush=True)
        J=np.stack(columns,axis=1);H=(J+J.T)/2;np.save(a.out/'hessian_numerical_matrix.npy',H)
        checks=[]
        probe=next((i for i,s in enumerate(symbols) if s=='O'),next((i for i,s in enumerate(symbols) if s=='Fe'),0))
        for k in (1,3*probe+2):
            ia,axis=divmod(k,3);grads=[]
            for sign in (-1,1):
                c=coordinates.copy();c[ia,axis]+=sign*.001
                _,fm,g=solve(c,use_warm=False);grads.append(g);del fm;_free_gpu_pool()
            rough=(grads[1]-grads[0])/.002
            checks.append(dict(column=k,step_and_cold_relative_error=float(np.linalg.norm(rough.ravel()-J[:,k])/np.linalg.norm(J[:,k]))))
        numerical_tensor=H.reshape(len(atoms),3,len(atoms),3).transpose(0,2,1,3)
        # Recreate the actual center molecule; do not project at the last displacement.
        mol=gto.M(atom=list(zip(symbols,coordinates.tolist())),unit='Bohr',charge=q,spin=1,basis=basis,verbose=0)
        data_n=thermo.harmonic_analysis(mol,numerical_tensor,imaginary_freq=False)
        data_a=thermo.harmonic_analysis(mol,HA,imaginary_freq=False)
        thermal={};deltas={}
        for kind,mat in (('numerical',H),('analytical',(A+A.T)/2)):
            work=a.out/kind;work.mkdir();(work/'mol.xyz').write_bytes(xyz.read_bytes());write_xtb_hessian(work/'hessian',mat)
            by_cutoff={}
            for cutoff in (25,50,100):
                with (work/f'thermo_{cutoff}.log').open('w') as f:
                    proc=subprocess.run(['/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb','thermo','--sthr',str(cutoff),'--temp','298.15','mol.xyz','hessian'],cwd=work,stdout=f,stderr=subprocess.STDOUT,timeout=120)
                if proc.returncode:raise RuntimeError('control qRRHO failed')
                by_cutoff[str(cutoff)]=parse_xtb_thermal_eV((work/f'thermo_{cutoff}.log').read_text())
            thermal[kind]=by_cutoff
        deltas={k:thermal['analytical'][k]-thermal['numerical'][k] for k in thermal['numerical']}
        error=float(np.linalg.norm(A-H)/np.linalg.norm(H));asym=float(abs(J-J.T).max())
        freq_n=np.asarray(data_n['freq_wavenumber']);freq_a=np.asarray(data_a['freq_wavenumber'])
        status.update(relative_full_matrix_error=error,numerical_max_asymmetry=asym,step_cold_checks=checks,
            numerical_frequencies_cm=freq_n.tolist(),analytical_frequencies_cm=freq_a.tolist(),
            qrrho=thermal,thermal_difference_eV=deltas)
        passed=(error<.01 and asym<1e-4 and all(r['step_and_cold_relative_error']<.01 for r in checks)
            and min(freq_n)>=-5 and min(freq_a)>=-5 and max(abs(x) for x in deltas.values())<.005)
        status.update(status='uks_full_response_verified' if passed else 'uks_full_response_failed',finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
