#!/usr/bin/env python
"""Checkpointed solution DFT Hessian from central differences of full SMD+D4 gradients.

Single existing singlet case. No production geometry or thermochemical record is overwritten.
Analytic Hessian failed two-size independent gradient tests; this is the numerical reference.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import time

import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo

from redox.core.common import DFT, ROOT, HARTREE_EV
from redox.qm.dft import _build_smd_mf, _free_gpu_pool

OUT=ROOT/'results/validation/physical_thermal_spin_2026-10-10/numerical_solution_hessian_aq_red2_tight'
XYZ=DFT/'aq_benzylamino/red2/opt.xyz'
STEP=5e-4


def write_hessian(path,matrix):
    from redox.qm.vibrations import write_xtb_hessian
    write_xtb_hessian(path,matrix)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    atoms=read(XYZ);symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
    code_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    identity=dict(geometry_sha256=hashlib.sha256(XYZ.read_bytes()).hexdigest(),charge=-2,mult=1,
                  xc='r2scan',basis='def2-svpd',disp='d4',solvent='SMD(acetonitrile)',
                  density_fit=True,conv_tol=1e-12,conv_tol_grad=1e-8,step_Bohr=STEP,code_sha256=code_hash)
    status=dict(status='running',identity=identity,host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().split()[21],
                cpu_affinity=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),
                started=time.time(),columns_complete=[],expected_columns=3*len(atoms))
    original_source=ROOT/'results/validation/physical_thermal_spin_2026-10-10/numerical_hessian_v1_source.py'
    original_hash=hashlib.sha256(original_source.read_bytes()).hexdigest() if original_source.exists() else None
    def save():
        p=OUT/'status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    lock=(OUT/'writer.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    save()
    baseline_dm=None
    def gradient(pos,use_warm=True):
        nonlocal baseline_dm
        mol=gto.M(atom=list(zip(symbols,pos.tolist())),unit='Bohr',charge=-2,spin=0,
                  basis='def2-svpd',verbose=0,max_memory=8000)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200
        # A converged center density warm-starts tiny Cartesian displacements.
        # Every displaced SCF still converges independently in its own AO basis.
        e=float(mf.kernel(dm0=baseline_dm) if baseline_dm is not None and use_warm else mf.kernel())
        if not mf.converged:raise RuntimeError('unconverged SCF in finite-difference reference')
        raw=mf.nuc_grad_method().kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw)
        if not np.isfinite(g).all():raise RuntimeError('nonfinite gradient')
        if baseline_dm is None:
            dm=mf.make_rdm1();baseline_dm=np.asarray(dm.get() if hasattr(dm,'get') else dm)
        del mf;_free_gpu_pool()
        return e,g,mol
    try:
        base_e,base_g,mol=gradient(coords);status.update(base_e_Ha=base_e,base_gradient_max=float(abs(base_g).max()))
        np.save(OUT/'base_gradient.npy',base_g);save()
        # Check warm-start derivatives against the two independently cold-started
        # columns already computed for the analytic Hessian validation.
        pilot_checks=[]
        previous=ROOT/'results/validation/physical_thermal_spin_2026-10-10/rks_hessian_aq_red2_v2'
        for ia,axis in ((1,1),(16,2)):
            gs=[];cold_gs=[]
            for sign in (-1,1):
                pos=coords.copy();pos[ia,axis]+=sign*STEP
                _,g,_=gradient(pos);gs.append(g)
                _,cg,_=gradient(pos,use_warm=False);cold_gs.append(cg)
            warm=(gs[1]-gs[0])/(2*STEP)
            cold=(cold_gs[1]-cold_gs[0])/(2*STEP)
            np.save(OUT/f'cold_reference_atom{ia}_axis{axis}.npy',cold)
            error=float(np.linalg.norm(warm-cold)/np.linalg.norm(cold))
            pilot_checks.append(dict(atom=ia,axis=axis,relative_warm_cold_error=error))
        status['warm_start_checks']=pilot_checks;save()
        if any(row['relative_warm_cold_error']>1e-3 for row in pilot_checks):
            raise RuntimeError('warm start altered independently verified derivative')
        for column in range(3*len(atoms)):
            p=OUT/f'column_{column:03d}.json';npfile=p.with_suffix('.npy')
            if p.exists() and npfile.exists():
                old=json.loads(p.read_text())
                before=dict(old.get('identity',{}));producer=before.pop('code_sha256',None)
                current=dict(identity);current.pop('code_sha256')
                if before!=current or producer not in (original_hash,code_hash):
                    raise RuntimeError(f'incompatible cached column {p}')
            else:
                ia,axis=divmod(column,3);gradients=[];energies=[];started=time.time()
                for sign in (-1,1):
                    pos=coords.copy();pos[ia,axis]+=sign*STEP
                    energy,g,_=gradient(pos);energies.append(energy);gradients.append(g)
                numeric=(gradients[1]-gradients[0])/(2*STEP)
                np.save(npfile,numeric)
                p.write_text(json.dumps(dict(identity=identity,column=column,energies_Ha=energies,
                                             seconds=time.time()-started),indent=2))
                print(f'COLUMN_COMPLETE {column+1}/{3*len(atoms)} seconds={time.time()-started:.1f}',flush=True)
            status['columns_complete'].append(column);save()
        J=np.stack([np.load(OUT/f'column_{i:03d}.npy').ravel() for i in range(3*len(atoms))],axis=1)
        status['hessian_max_asymmetry_Ha_Bohr2']=float(abs(J-J.T).max())
        H=(J+J.T)/2;np.save(OUT/'hessian_matrix.npy',H)
        tensor=H.reshape(len(atoms),3,len(atoms),3).transpose(0,2,1,3)
        analysis=thermo.harmonic_analysis(mol,tensor,exclude_trans=True,exclude_rot=True,imaginary_freq=False)
        freq=np.asarray(analysis['freq_wavenumber']);np.save(OUT/'normal_modes.npy',analysis['norm_mode'])
        status.update(projected_frequencies_cm=freq.tolist(),minimum_frequency_cm=float(freq.min()),
                      n_imag_lt_minus5=int((freq<-5).sum()),stage='thermal_postprocessing')
        write_hessian(OUT/'hessian',H)
        (OUT/'mol.xyz').write_bytes(XYZ.read_bytes())
        thermal={}
        for cutoff in (25,50,100):
            log=OUT/f'thermo_{cutoff}.log'
            argv=[str(Path(os.environ.get('XTB_BIN','/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb'))),
                  'thermo','--sthr',str(cutoff),'--temp','298.15','mol.xyz','hessian']
            with log.open('w') as f:proc=subprocess.run(argv,cwd=OUT,stdout=f,stderr=subprocess.STDOUT,timeout=120)
            if proc.returncode:raise RuntimeError(f'xTB thermo failed: {log}')
            from redox.qm.vibrations import parse_xtb_thermal_eV
            thermal[str(cutoff)]=parse_xtb_thermal_eV(log.read_text())
        checks=[]
        for ia,axis in ((1,1),(16,2)):
            gs=[]
            for sign in (-1,1):
                pos=coords.copy();pos[ia,axis]+=sign*1e-3
                _,g,_=gradient(pos);gs.append(g)
            coarse=(gs[1]-gs[0])/(2e-3)
            fine=np.load(OUT/f'column_{3*ia+axis:03d}.npy')
            checks.append(dict(atom=ia,axis=axis,relative_step_error=float(np.linalg.norm(fine-coarse)/np.linalg.norm(fine))))
        status.update(status='verified_numerical' if all(c['relative_step_error']<.01 for c in checks)
                      and status['hessian_max_asymmetry_Ha_Bohr2']<1e-3 else 'validation_failed',
                      step_checks=checks,g_thermal_qrrho_by_cutoff_eV=thermal,
                      stationary_minimum=bool(freq.min()>=-5),finished=time.time(),production_modified=False)
        save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise
    finally:lock.close()


if __name__=='__main__':main()
