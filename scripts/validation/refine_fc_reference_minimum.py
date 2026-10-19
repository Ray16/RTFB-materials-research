#!/usr/bin/env python
"""Resolve the actual Fc+ solution-surface instability without replacing production."""
import hashlib,json,os,socket,subprocess,time,signal
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from pyscf.hessian import thermo
from pyscf.geomopt.geometric_solver import optimize
from redox.core.common import ROOT,DFT
from redox.qm.dft import _build_smd_mf,_free_gpu_pool,_orbital_diag,dft_smd
from redox.qm.vibrations import write_xtb_hessian,parse_xtb_thermal_eV
from redox.core.protocol import load_record,software_versions
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    out=BASE/'fc_conservative_minimum_refinement';out.mkdir(exist_ok=True)
    xyz=DFT/'ferrocene/ox/opt.xyz';atoms=read(xyz);n=len(atoms)
    source=BASE/'uks_fc_full_response_control';r=json.loads((source/'status.json').read_text())
    if r['geometry_sha256']!=hashlib.sha256(xyz.read_bytes()).hexdigest():raise RuntimeError('source reference changed')
    H=np.load(source/'hessian_numerical_matrix.npy')
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),charge=1,spin=1,basis='def2-svp',verbose=0)
    modes=thermo.harmonic_analysis(mol,H.reshape(n,3,n,3).transpose(0,2,1,3),imaginary_freq=False)
    k=int(np.argmin(modes['freq_wavenumber']));mode=np.asarray(modes['norm_mode'][k]);mode/=np.max(np.linalg.norm(mode,axis=1))
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
        cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),
        source_sha256=r['geometry_sha256'],source_frequency_cm=float(modes['freq_wavenumber'][k]),
        max_displacement_A=.1,production_modified=False,branches=[])
    def save():
        p=out/'status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    save()
    def stop(signum,frame):raise KeyboardInterrupt(f'termination signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    def cpu(v):return np.asarray(v.get() if hasattr(v,'get') else v)
    import gpu4pyscf.df.hessian.rhf as dr
    import gpu4pyscf.df.hessian.uhf as du
    for module in (dr,du):
        if hasattr(module,'get_avail_mem'):
            previous=module.get_avail_mem
            module.get_avail_mem=lambda *args,_old=previous,**kwargs:min(_old(*args,**kwargs),6*1024**3)
    def build(coords):
        molecule=gto.M(atom=list(zip(atoms.get_chemical_symbols(),coords.tolist())),charge=1,spin=1,basis='def2-svp',verbose=0,max_memory=6000)
        mf=_build_smd_mf(molecule,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=3
        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
        return molecule,mf
    try:
        # Isolate response-grid accuracy from the separate geometry instability.
        _,testmf=build(atoms.positions);testmf.kernel()
        if not testmf.converged:raise RuntimeError('metal response-grid SCF failed')
        testmf.cphf_grids=testmf.grids.copy();hh=testmf.Hessian();hh.grid_response=True;hh.max_cycle=100
        fine=cpu(hh.kernel()).transpose(0,2,1,3).reshape(3*n,3*n)
        np.save(out/'fine_response_hessian_original.npy',fine)
        status['fine_response_relative_matrix_error_original']=float(np.linalg.norm(fine-H)/np.linalg.norm(H))
        status['fine_response_implementation_verified']=status['fine_response_relative_matrix_error_original']<.01
        del testmf,hh;_free_gpu_pool();save()
        for sign in (-1,1):
            work=out/('minus' if sign<0 else 'plus');work.mkdir()
            seed=atoms.copy();seed.positions+=sign*.1*mode;write(work/'seed.xyz',seed)
            _,mf=build(seed.positions);mf.kernel()
            if not mf.converged:raise RuntimeError('source SCF failed')
            grad=mf.nuc_grad_method();grad.grid_response=True
            scanner=grad.as_scanner()
            optimized=optimize(scanner,maxsteps=120,assert_convergence=True,
                convergence_energy=1e-8,convergence_grms=3e-5,convergence_gmax=5e-5,
                convergence_drms=1e-4,convergence_dmax=2e-4)
            final=atoms.copy();final.positions=optimized.atom_coords(unit='Angstrom');write(work/'opt.xyz',final)
            del mf,grad,scanner;_free_gpu_pool()
            mol,mf=build(final.positions);mf.kernel()
            if not mf.converged:raise RuntimeError('refined SCF failed')
            gobj=mf.nuc_grad_method();gobj.grid_response=True;g=cpu(gobj.kernel());diag=_orbital_diag(mf,1)
            # Test fine CPHF integration; coarse response grids were inaccurate for this metal.
            mf.cphf_grids=mf.grids.copy()
            hobj=mf.Hessian();hobj.grid_response=True;hobj.max_cycle=100
            ha=cpu(hobj.kernel());A=ha.transpose(0,2,1,3).reshape(3*n,3*n)
            np.save(work/'hessian_matrix.npy',A)
            data=thermo.harmonic_analysis(mol,ha,imaginary_freq=False);freq=np.asarray(data['freq_wavenumber'])
            e_opt=float(mf.e_tot);del mf,hobj,gobj;_free_gpu_pool()
            checks=[]
            for ia,axis in ((0,1),(next(i for i,s in enumerate(atoms.get_chemical_symbols()) if s=='Fe'),2)):
                for step in (.0005,.001):
                    gs=[]
                    for direction in (-1,1):
                        coords=final.positions.copy();coords[ia,axis]+=direction*step*.52917721067
                        _,fm=build(coords);fm.kernel()
                        if not fm.converged:raise RuntimeError('independent gradient failed')
                        gg=fm.nuc_grad_method();gg.grid_response=True;gs.append(cpu(gg.kernel()));del fm,gg;_free_gpu_pool()
                    numerical=(gs[1]-gs[0])/(2*step);err=float(np.linalg.norm(numerical-ha[ia,:,axis,:])/np.linalg.norm(numerical))
                    checks.append(dict(atom=ia,axis=axis,step_Bohr=step,relative_error=err))
            values={};write_xtb_hessian(work/'hessian',(A+A.T)/2)
            for cutoff in (25,50,100):
                log=work/f'thermo_{cutoff}.log'
                with log.open('w') as f:
                    proc=subprocess.run(['/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb','thermo',
                        '--sthr',str(cutoff),'--temp','298.15','opt.xyz','hessian'],cwd=work,stdout=f,stderr=subprocess.STDOUT,timeout=120)
                if proc.returncode:raise RuntimeError('refined qRRHO failed')
                values[str(cutoff)]=parse_xtb_thermal_eV(log.read_text())
            passed=(min(freq)>=-5 and max(abs(g).ravel())<5e-5 and all(c['relative_error']<.01 for c in checks)
                and abs(diag.get('spin_contam',99))<.1)
            record=dict(sign=sign,geometry=str(work/'opt.xyz'),electronic_opt_energy_Ha=e_opt,
                frequencies_cm=freq.tolist(),minimum_frequency_cm=float(freq.min()),gradient_max=float(abs(g).max()),
                smd_spin_diagnostics=diag,derivative_checks=checks,nuclear_by_cutoff_eV=values,
                cphf_grid_level=3,minimum_and_derivatives_verified=passed,software=software_versions(),production_modified=False)
            if passed:
                energy=dft_smd(work/'opt.xyz',1,2,backend='gpu',do_opt=False,do_freq=False)
                old=load_record(DFT/'ferrocene/ox',1,2)
                record.update(active_protocol_energy=energy,delta_active_smd_eV=energy['e_smd_eV']-old['e_smd_eV'])
            (work/'result.json').write_text(json.dumps(record,indent=2));status['branches'].append(record);save()
            print('FC_REFINED_BRANCH',sign,record['minimum_frequency_cm'],record['gradient_max'],passed,flush=True)
        good=[r for r in status['branches'] if r['minimum_and_derivatives_verified']]
        status.update(status='refined_reference_verified' if good else 'refinement_failed',finished=time.time())
        if good:status['selected']=min(good,key=lambda r:r['electronic_opt_energy_Ha'])
        save()
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
