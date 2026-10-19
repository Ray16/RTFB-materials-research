#!/usr/bin/env python
"""Bounded preserved solution-model Hessians and mode following for selected states."""
import argparse,hashlib,json,os,re,socket,subprocess,time
from pathlib import Path
import numpy as np
from ase.io import read,write
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import ROOT,DFT,read_manifest,read_result,smiles_to_xyz_map
from redox.properties.integrity import connectivity_change,kabsch_rmsd_heavy
from redox.qm.thermal import SOLUTION_THERMAL,thermal_identity,thermal_hash,spectrum_values,stationary_verdict
from redox.qm.vibrations import parse_xtb_thermal_eV
XTB='/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--only',required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--threads',type=int,default=2);a=ap.parse_args()
    if len(os.sched_getaffinity(0))!=a.threads or os.getloadavg()[0]+a.threads>os.cpu_count():raise RuntimeError('CPU capacity/affinity insufficient')
    ids=set(a.only.split(','))|{'ferrocene'};rows=[r for r in read_manifest() if r['id'] in ids]
    if {r['id'] for r in rows}!=ids:raise ValueError('unknown selected molecule')
    a.out.mkdir(parents=True,exist_ok=True)
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],
        affinity=sorted(os.sched_getaffinity(0)),started=time.time(),protocol=SOLUTION_THERMAL,
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),states=[],production_modified=False)
    def save():
        p=a.out/'status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    save()
    def branch(row,identity,path,atoms,displacement=None):
        p=path/'result.json'
        if p.exists():
            rec=json.loads(p.read_text())
            if rec['identity']!=identity:raise RuntimeError('existing thermal identity mismatch')
            return rec
        if path.exists():raise RuntimeError(f'partial thermal output requires inspection: {path}')
        path.mkdir(parents=True);write(path/'seed.xyz',atoms)
        (path/'thermal.inp').write_text('$thermo\n temp=298.15\n sthr=50\n$end\n')
        argv=[XTB,'seed.xyz','--gfn','2','--chrg',str(identity['charge']),'--uhf',str(identity['mult']-1),
            '--alpb','acetonitrile','--ohess','verytight','--acc','0.1','--input','thermal.inp']
        print('START',row['id'],row['state'],path.name,flush=True)
        with (path/'xtb.log').open('w') as f:
            proc=subprocess.run(argv,cwd=path,stdout=f,stderr=subprocess.STDOUT,timeout=900)
        text=(path/'xtb.log').read_text()
        if proc.returncode or not (path/'xtbopt.xyz').exists() or not (path/'hessian').exists():raise RuntimeError(f'xTB incomplete {path}')
        frequencies=spectrum_values((path/'vibspectrum').read_text());assert len(frequencies)==3*len(atoms)
        n=re.findall(r'#\s*imaginary freq\.\s+(\d+)',text);assert n
        src=DFT/row['id']/row['state']/'opt.xyz';geom=path/'xtbopt.xyz'
        intact=connectivity_change(src,geom)==0 and smiles_to_xyz_map(row['smiles'],geom) is not None
        energies=re.findall(r'TOTAL ENERGY\s+([-+0-9.Ee]+)',text)
        rec=dict(identity=identity,id=row['id'],state=row['state'],returncode=proc.returncode,
            optimization_converged='GEOMETRY OPTIMIZATION CONVERGED' in text.upper(),connectivity_intact=intact,
            frequencies_cm=frequencies,n_imag=int(n[-1]),g_thermal_eV=parse_xtb_thermal_eV(text),
            total_energy_Ha=float(energies[-1]),geom=str(geom),hessian=str(path/'hessian'),stdout=str(path/'xtb.log'),
            branch_displacement=displacement,heavy_rmsd_from_DFT_A=kabsch_rmsd_heavy(src,geom),
            transfer_accuracy_validated=False,production_modified=False)
        rec['stationary_valid'],rec['issues']=stationary_verdict(rec,identity)
        p.write_text(json.dumps(rec,indent=2));return rec
    def imaginary_seed(rec):
        atoms=read(rec['geom']);n=len(atoms)
        values=[float(v) for l in Path(rec['hessian']).read_text().splitlines() if not l.startswith('$') for v in l.split()]
        H=np.array(values).reshape(3*n,3*n).T
        mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='sto-3g',
            charge=rec['identity']['charge'],spin=rec['identity']['mult']-1,verbose=0)
        data=thermo.harmonic_analysis(mol,((H+H.T)/2).reshape(n,3,n,3).transpose(0,2,1,3),imaginary_freq=False)
        k=int(np.argmin(data['freq_wavenumber']));mode=np.array(data['norm_mode'][k]);mode/=np.max(np.linalg.norm(mode,axis=1))
        return atoms,mode,float(data['freq_wavenumber'][k])
    try:
        for row in rows:
            src=DFT/row['id']/row['state']/'opt.xyz';raw=read_result(row['id'],row['state'],raw=True)
            identity=thermal_identity(src,raw['charge'],raw['mult']);sd=a.out/row['id']/row['state']/thermal_hash(identity)
            rec=branch(row,identity,sd/'initial',read(src));cases=[rec]
            for level in range(2):
                invalid=[r for r in cases if r['connectivity_intact'] and r['optimization_converged'] and not r['stationary_valid']]
                if not invalid:break
                parent=min(invalid,key=lambda r:r['total_energy_Ha']);atoms,mode,freq=imaginary_seed(parent)
                for sign in (-1,1):
                    seed=atoms.copy();seed.positions+=sign*.1*mode
                    cases.append(branch(row,identity,sd/f'follow{level+1}_{sign:+}',seed,
                        dict(parent_hessian_sha256=hashlib.sha256(Path(parent['hessian']).read_bytes()).hexdigest(),
                            projected_frequency_cm=freq,max_displacement_A=.1,sign=sign)))
                stable=[r for r in cases if r['stationary_valid']]
                if stable and min(r['total_energy_Ha'] for r in stable)<=min(r['total_energy_Ha'] for r in invalid)+1e-5:break
            stable=[r for r in cases if r['stationary_valid']]
            result=dict(id=row['id'],state=row['state'],identity=identity,branches=cases,status='stable_model_minimum' if stable else 'no_stable_model_minimum')
            if stable:
                chosen=min(stable,key=lambda r:r['total_energy_Ha']);result['selected']=chosen
                result['stable_branch_thermal_spread_eV']=max(r['g_thermal_eV'] for r in stable)-min(r['g_thermal_eV'] for r in stable)
                work=Path(chosen['geom']).parent;thermal={}
                for cutoff in (25,50,100):
                    log=work/f'precise_thermo_{cutoff}.log'
                    with log.open('w') as f:
                        p=subprocess.run([XTB,'thermo','--sthr',str(cutoff),'--temp','298.15','xtbopt.xyz','hessian'],cwd=work,stdout=f,stderr=subprocess.STDOUT,timeout=120)
                    if p.returncode:raise RuntimeError('thermo postprocessing failed')
                    thermal[str(cutoff)]=parse_xtb_thermal_eV(log.read_text())
                result['g_thermal_by_cutoff_eV']=thermal
            (sd/'selected.json').write_text(json.dumps(result,indent=2));status['states'].append(dict(id=row['id'],state=row['state'],status=result['status'],selected=str(sd/'selected.json')));save()
        status.update(status='model_validation_complete',finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
