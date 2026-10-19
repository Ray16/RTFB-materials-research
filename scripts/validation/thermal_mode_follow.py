#!/usr/bin/env python
"""Bounded GFN2-ALPB stationary-point follow-up for two problematic AQ states.

Moves in both directions along the projected imaginary mode, then optimizes and
recomputes the Hessian at tighter electronic accuracy. No production replacement.
"""
import hashlib,json,os,re,socket,subprocess,time
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,DFT,read_manifest,read_result
from redox.core.common import smiles_to_xyz_map
from redox.properties.integrity import connectivity_change
from redox.qm.vibrations import parse_xtb_thermal_eV

BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
OUT=BASE/'mode_follow_2026-10-11'
XTB='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    summary=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
                 affinity=sorted(os.sched_getaffinity(0)),production_modified=False,states=[])
    def save():
        p=OUT/'status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(summary,indent=2));t.replace(p)
    save()
    try:
        for state in ('neu','red2'):
            sd=next((BASE/'stationary/aq_benzylamino'/state).glob('alpb_relaxed__*'))
            parent=json.loads((sd/'result.json').read_text())
            if not parent['connectivity_intact']:raise RuntimeError('source is chemically changed')
            atoms=read(parent['geom']);q,m=parent['charge'],parent['mult'];n=len(atoms)
            numbers=[float(v) for line in Path(parent['hessian']).read_text().splitlines()
                     if not line.startswith('$') for v in line.split()]
            H=np.array(numbers).reshape(3*n,3*n).T
            mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),
                      basis='sto-3g',charge=q,spin=m-1,verbose=0)
            data=thermo.harmonic_analysis(mol,((H+H.T)/2).reshape(n,3,n,3).transpose(0,2,1,3),
                                         exclude_trans=True,exclude_rot=True,imaginary_freq=False)
            index=int(np.argmin(data['freq_wavenumber']));mode=np.asarray(data['norm_mode'][index])
            mode/=np.max(np.linalg.norm(mode,axis=1))
            cases=[]
            for sign in (-1,1):
                work=OUT/state/('minus' if sign<0 else 'plus');result=work/'result.json'
                if result.exists():cases.append(json.loads(result.read_text()));continue
                if work.exists():raise RuntimeError(f'partial branch must be inspected: {work}')
                work.mkdir(parents=True)
                seed=atoms.copy();seed.positions+=sign*.10*mode;write(work/'seed.xyz',seed)
                (work/'thermo.inp').write_text('$thermo\n temp=298.15\n sthr=50\n$end\n')
                argv=[XTB,'seed.xyz','--gfn','2','--chrg',str(q),'--uhf',str(m-1),'--alpb','acetonitrile',
                      '--ohess','verytight','--acc','0.1','--input','thermo.inp']
                with (work/'xtb.log').open('w') as file:
                    proc=subprocess.run(argv,cwd=work,stdout=file,stderr=subprocess.STDOUT,timeout=1200)
                text=(work/'xtb.log').read_text()
                if proc.returncode or not (work/'hessian').exists() or not (work/'xtbopt.xyz').exists():
                    raise RuntimeError(f'branch incomplete: {work}')
                freq=[]
                for line in (work/'vibspectrum').read_text().splitlines():
                    vals=line.split()
                    if vals and vals[0].isdigit():
                        for token in vals[1:]:
                            try:freq.append(float(token));break
                            except ValueError:continue
                if len(freq)!=3*n:raise RuntimeError('spectrum dimension mismatch')
                ni=re.findall(r'#\s*imaginary freq\.\s+(\d+)',text)
                if not ni:raise RuntimeError('missing imaginary-mode count')
                intact=connectivity_change(DFT/'aq_benzylamino'/state/'opt.xyz',work/'xtbopt.xyz')==0
                smi=next(r['smiles'] for r in read_manifest() if r['id']=='aq_benzylamino')
                intact=intact and smiles_to_xyz_map(smi,work/'xtbopt.xyz') is not None
                converged='GEOMETRY OPTIMIZATION CONVERGED' in text.upper()
                energy=re.findall(r'TOTAL ENERGY\s+([-+0-9.Ee]+)',text)
                rec=dict(id='aq_benzylamino',state=state,displacement_sign=sign,
                         source_hessian_sha256=hashlib.sha256(Path(parent['hessian']).read_bytes()).hexdigest(),
                         projected_source_frequency_cm=float(data['freq_wavenumber'][index]),
                         seed_max_displacement_A=.10,argv=argv,model='GFN2-ALPB(MeCN)',
                         total_energy_Ha=float(energy[-1]) if energy else None,
                         n_imag=int(ni[-1]),minimum_frequency_cm=min(freq),connectivity_intact=intact,
                         optimization_converged=converged,g_thermal_eV=parse_xtb_thermal_eV(text),
                         stationary_hessian_valid=converged and intact and int(ni[-1])==0 and min(freq)>=-5,
                         high_level_thermal_accuracy_validated=False,production_modified=False)
                result.write_text(json.dumps(rec,indent=2)+'\n');cases.append(rec)
                print('MODE_BRANCH_COMPLETE',state,sign,rec['minimum_frequency_cm'],rec['stationary_hessian_valid'],flush=True)
            summary['states'].append(dict(state=state,branches=cases));save()
        summary.update(status='model_stationary_checks_completed',finished=time.time());save()
    except Exception as exc:
        summary.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
