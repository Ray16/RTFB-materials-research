#!/usr/bin/env python
"""Phase-matched occupation seeds on preserved refined Fc+ reference geometries."""
import hashlib,importlib.util,json,os,socket,time
from pathlib import Path
import numpy as np
import cupy
from redox.core.common import ROOT,DFT,HARTREE_EV
from redox.core.protocol import ACTIVE_SP,load_record,software_versions
from redox.qm.dft import _kernel_robust,_orbital_diag,_free_gpu_pool
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'


def main():
    spec=importlib.util.spec_from_file_location('S',ROOT/'scripts/validation/reference/scf_level_check.py')
    S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
    source=BASE/'fc_conservative_minimum_refinement';summary=json.loads((source/'status.json').read_text())
    if summary['status']!='refined_reference_verified':raise RuntimeError('refined minima not verified')
    out=BASE/'refined_fc_smd_occupation_checks';out.mkdir(exist_ok=True)
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
        cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),protocol=ACTIVE_SP,
        software=software_versions(),production_modified=False,cases=[])
    p=out/'status.json'
    def save():
        tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    save()
    try:
        cases=[('original',DFT/'ferrocene/ox/opt.xyz',None)]+[
            ('minus' if r['sign']<0 else 'plus',Path(r['geometry']),r)
            for r in summary['branches'] if r['minimum_and_derivatives_verified']]
        old=load_record(DFT/'ferrocene/ox',1,2)
        for name,geometry,branch in cases:
            ref=S._mf(S._mol(geometry,0,0,ACTIVE_SP['basis']),S.LEVELS['prod'],True,'gpu')
            _,conv=_kernel_robust(ref)
            if not conv:raise RuntimeError('phase-matched neutral occupation reference failed')
            C=ref.mo_coeff;occ=cupy.asnumpy(ref.mo_occ);nocc=int((occ>0).sum())
            target=S._mf(S._mol(geometry,1,1,ACTIVE_SP['basis']),S.LEVELS['prod'],True,'gpu')
            densities={}
            for k in (0,1):
                oa=occ/2;ob=oa.copy();ob[nocc-1-k]=0
                densities[f'smd_rm{k}']=cupy.asnumpy(target.make_rdm1((C,C),(cupy.asarray(oa),cupy.asarray(ob))))
            del ref,target,C;_free_gpu_pool()
            guesses=[]
            for label,dm in densities.items():
                mf=S._mf(S._mol(geometry,1,1,ACTIVE_SP['basis']),S.LEVELS['prod'],True,'gpu')
                e,conv=_kernel_robust(mf,dm0=dm);diag=_orbital_diag(mf,1)
                guesses.append(dict(guess=label,e_Ha=e,e_eV=e*HARTREE_EV,converged=conv,smd_diagnostics=diag))
                del mf;_free_gpu_pool()
            valid=[r for r in guesses if r['converged'] and abs(r['smd_diagnostics'].get('spin_contam',99))<.1]
            if not valid:raise RuntimeError('no valid occupation-seeded SMD minimum')
            best=min(valid,key=lambda r:r['e_Ha'])
            record=dict(name=name,geometry=str(geometry),geometry_sha256=hashlib.sha256(geometry.read_bytes()).hexdigest(),
                guesses=guesses,best=best,delta_smd_from_original_eV=best['e_eV']-old['e_smd_eV'],
                standard_three_guess_energy_eV=branch['active_protocol_energy']['e_smd_eV'] if branch else old['e_smd_eV'])
            record['three_guess_above_seeded_eV']=record['standard_three_guess_energy_eV']-best['e_eV']
            (out/f'{name}.json').write_text(json.dumps(record,indent=2));status['cases'].append(record);save()
            print('SMD_OCCUPATION_CHECK',name,best['e_eV'],record['three_guess_above_seeded_eV'],flush=True)
        status.update(status='smd_occupation_checks_verified',finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
