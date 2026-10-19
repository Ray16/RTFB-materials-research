#!/usr/bin/env python
"""Check the metal basis/ECP determinant omitted by the earlier OROP level crossing.

Diagnostic only: existing production SMD geometries, assumed D3(0), and PySCF default
PCM cavity. The source states LANL2DZ for metals; exact TeraChem cavity/geometry is
not reproduced here. No production/reference records are changed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time
import types

import numpy as np
from pyscf import gto, dft
from redox.core.common import DFT, HARTREE_EV, read_result


def cached_cpcm(cm):
    """Native symmetric C-PCM energy/field with one cavity factorization per geometry."""
    from scipy.linalg import lu_factor,lu_solve
    cm.build();K=cm._intermediates['K'];f=(cm.eps-1)/cm.eps
    factor=lu_factor(K,check_finite=False)
    if not np.allclose(K,K.T,atol=1e-12):raise RuntimeError('C-PCM matrix is not symmetric')
    def vind(self,density):
        nao=density.shape[-1];dm=np.asarray(density).reshape(-1,nao,nao)
        if len(dm)==2:dm=(dm[0]+dm[1]).reshape(1,nao,nao)
        V=self.v_grids_n-self._get_v(dm)
        q=-f*lu_solve(factor,V.T,check_finite=False).T
        energy=.5*np.dot(q[0],V[0]);matrix=self._get_vmat(q)[0]
        self._intermediates.update(q=q[0],q_sym=q[0],v_grids=V[0],dm=dm)
        return energy,matrix
    cm._get_vind=types.MethodType(vind,cm)
    return cm


def validate_cached_cpcm():
    from pyscf import scf
    from pyscf.solvent.pcm import PCM
    mol=gto.M(atom='O 0 0 0; H 0 0 .96; H .92 0 -.24',basis='sto-3g',verbose=0)
    dm=scf.RHF(mol).get_init_guess()
    native=PCM(mol);native.method='C-PCM';native.eps=37.5;native.lebedev_order=11
    fast=PCM(mol);fast.method='C-PCM';fast.eps=37.5;fast.lebedev_order=11;cached_cpcm(fast)
    residuals=[]
    for scale in (1,.95):
        en,vn=native._get_vind(dm*scale);ef,vf=fast._get_vind(dm*scale)
        residuals.append(dict(energy_Ha=float(abs(en-ef)),field_Ha=float(abs(vn-vf).max())))
    if any(max(r.values())>1e-11 for r in residuals):raise RuntimeError('CPCM LU validation failed')
    return residuals


def solve_state(state):
    xyz=DFT/'ferrocene'/state/'opt.xyz'
    q,m=(0,1) if state=='neu' else (1,2)
    lines=xyz.read_text().splitlines();atoms='\n'.join(lines[2:2+int(lines[0])])
    def build(charge,spin):
        mol=gto.M(atom=atoms,basis={'default':'6-31g*','Fe':'lanl2dz'},
                  ecp={'Fe':'lanl2dz'},charge=charge,spin=spin,verbose=0,max_memory=8000)
        mf=(dft.UKS if spin else dft.RKS)(mol).density_fit().PCM()
        mf.with_solvent.method='C-PCM';mf.with_solvent.eps=37.5
        cached_cpcm(mf.with_solvent)
        mf.xc='b3lyp';mf.disp='d3zero';mf.conv_tol=1e-9;mf.max_cycle=200
        return mf
    guesses=[]
    for name in ('minao','atom','huckel') if m>1 else ('minao',):
        print('SCF_START',state,name,flush=True)
        mf=build(q,m-1);mf.init_guess=name;e=float(mf.kernel())
        guesses.append(dict(guess=name,e_Ha=e,converged=bool(mf.converged)))
        print('SCF_COMPLETE',state,guesses[-1],flush=True)
    if m==2:
        # Occupation seeds are real independent density starts, all at this same ECP level.
        parent=build(q-1,0);parent.kernel()
        if parent.converged:
            C=parent.mo_coeff;nocc=int((parent.mo_occ>0).sum())
            for k in (0,1):
                oa=parent.mo_occ/2;ob=oa.copy();ob[nocc-1-k]=0
                mf=build(q,1);dm=mf.make_rdm1((C,C),(oa,ob));e=float(mf.kernel(dm0=dm))
                guesses.append(dict(guess=f'rm{k}',e_Ha=e,converged=bool(mf.converged),s2=float(mf.spin_square()[0])))
                print('SCF_COMPLETE',state,guesses[-1],flush=True)
    valid=[g for g in guesses if g['converged'] and np.isfinite(g['e_Ha'])]
    if not valid:raise RuntimeError(f'{state}: no converged solution')
    best=min(valid,key=lambda x:x['e_Ha'])
    return dict(state=state,charge=q,mult=m,geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                guesses=guesses,best=best,thermal_eV=read_result('ferrocene',state)['g_thermal_eV'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    if a.out.exists():ap.error('choose a new output; prior diagnostic preserved')
    a.out.mkdir(parents=True)
    started=time.time();validation=validate_cached_cpcm();print('CPCM_REFERENCE_PASS',validation,flush=True)
    rows=[]
    for state in ('neu','ox'):
        row=solve_state(state);rows.append(row)
        (a.out/f'{state}.json').write_text(json.dumps(row,indent=2)+'\n')
    electronic=(rows[1]['best']['e_Ha']-rows[0]['best']['e_Ha'])*HARTREE_EV
    report=dict(status='completed',host=os.uname().nodename,pid=os.getpid(),
                affinity=sorted(os.sched_getaffinity(0)),elapsed_seconds=time.time()-started,
                source='https://doi.org/10.1021/acs.jctc.1c01040',
                method=dict(xc='b3lyp',disp='d3zero (assumed)',basis={'C,H':'6-31g*','Fe':'lanl2dz'},
                            ecp={'Fe':'lanl2dz'},density_fit=True,solvent='C-PCM eps37.5',
                            cavity='PySCF default modified Bondi; not exact TeraChem reproduction'),
                absolute_electronic_Fc_V=electronic,
                absolute_with_legacy_thermal_Fc_V=electronic+rows[1]['thermal_eV']-rows[0]['thermal_eV'],
                published_OROP_Fc_V=4.66224854035885,states=rows,production_modified=False,
                native_cpcm_solver_checks=validation)
    (a.out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
