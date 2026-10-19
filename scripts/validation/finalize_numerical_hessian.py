#!/usr/bin/env python
"""Bounded finalization of the numerical Hessian, including two-size gradient checks.

Waits for the exact owned producer; fixes only postprocessing format in a separate
directory. Invalid/missing gradient data never count as a completed reference.
"""
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time

import numpy as np
from ase.io import read
from pyscf import gto
from pyscf.hessian import thermo
from redox.core.common import DFT,ROOT
from redox.qm.dft import _build_smd_mf,_free_gpu_pool
from redox.qm.vibrations import write_xtb_hessian,parse_xtb_thermal_eV

BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
SOURCE=BASE/'numerical_solution_hessian_aq_red2_tight'
OUT=BASE/'numerical_solution_hessian_verified'
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'


def validate():
    # This function runs inside the GPU reservation after the original worker exits.
    inputs=json.loads((SOURCE/'status.json').read_text())
    geometry=DFT/'aq_benzylamino/red2/opt.xyz'
    identity=inputs['identity']
    if hashlib.sha256(geometry.read_bytes()).hexdigest()!=identity['geometry_sha256']:
        raise RuntimeError('geometry changed while the Hessian was computed')
    atoms=read(geometry);symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
    n=3*len(atoms)
    if len(inputs.get('columns_complete',[]))!=n:raise RuntimeError('numerical columns incomplete')
    columns=[]
    for i in range(n):
        info=json.loads((SOURCE/f'column_{i:03d}.json').read_text())
        column=np.load(SOURCE/f'column_{i:03d}.npy')
        if info['identity']!=identity or column.shape!=(len(atoms),3) or not np.isfinite(column).all():
            raise RuntimeError(f'invalid column {i}')
        columns.append(column.ravel())
    J=np.stack(columns,axis=1);asym=float(abs(J-J.T).max());H=(J+J.T)/2
    if asym>1e-3:raise RuntimeError(f'finite Hessian is asymmetric: {asym}')
    if any(r['relative_warm_cold_error']>1e-3 for r in inputs['warm_start_checks']):
        raise RuntimeError('independent initial-density check failed')
    checks=[]
    def grad(pos):
        mol=gto.M(atom=list(zip(symbols,pos.tolist())),unit='Bohr',charge=-2,spin=0,
                  basis='def2-svpd',verbose=0,max_memory=8000)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
        mf.grids.level=identity.get('grid_level',3)
        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;mf.kernel()
        if not mf.converged:raise RuntimeError('independent gradient SCF failed')
        gradient=mf.nuc_grad_method();gradient.grid_response=identity.get('grid_response',False)
        raw=gradient.kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw)
        del mf;_free_gpu_pool();return g,mol
    for ia,axis in ((1,1),(16,2)):
        gradients=[]
        for sign in (-1,1):
            x=coords.copy();x[ia,axis]+=sign*1e-3;g,mol=grad(x);gradients.append(g)
        coarse=(gradients[1]-gradients[0])/(2e-3)
        fine=np.load(SOURCE/f'column_{3*ia+axis:03d}.npy')
        error=float(np.linalg.norm(coarse-fine)/np.linalg.norm(fine))
        checks.append(dict(atom=ia,axis=axis,relative_step_error=error))
        if error>.01:raise RuntimeError(f'step-size sensitivity too high {checks[-1]}')
    # Harmonic projection requires center coordinates, not the last displaced geometry.
    mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=-2,spin=0,basis='def2-svpd',verbose=0)
    tensor=H.reshape(len(atoms),3,len(atoms),3).transpose(0,2,1,3)
    result=thermo.harmonic_analysis(mol,tensor,exclude_trans=True,exclude_rot=True,imaginary_freq=False)
    freq=np.asarray(result['freq_wavenumber'])
    OUT.mkdir(parents=True,exist_ok=True);np.save(OUT/'hessian_matrix.npy',H);np.save(OUT/'normal_modes.npy',result['norm_mode'])
    write_xtb_hessian(OUT/'hessian',H);(OUT/'mol.xyz').write_bytes(geometry.read_bytes())
    values={}
    for cutoff in (25,50,100):
        log=OUT/f'thermo_{cutoff}.log'
        argv=['/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb','thermo','--sthr',str(cutoff),
              '--temp','298.15','mol.xyz','hessian']
        with log.open('w') as file:p=subprocess.run(argv,cwd=OUT,stdout=file,stderr=subprocess.STDOUT,timeout=120)
        if p.returncode:raise RuntimeError(f'thermo failed: {log}')
        values[str(cutoff)]=parse_xtb_thermal_eV(log.read_text())
    rec=dict(status='numerical_reference_verified',identity=identity,
             maximum_asymmetry_Ha_Bohr2=asym,step_checks=checks,warm_start_checks=inputs['warm_start_checks'],
             projected_frequencies_cm=freq.tolist(),minimum_frequency_cm=float(freq.min()),
             stationary_minimum=bool(freq.min()>=-5),n_imag_lt_minus5=int((freq<-5).sum()),
             g_thermal_qrrho_by_cutoff_eV=values,
             source_run_status=inputs['status'],source_dir=str(SOURCE),
             production_modified=False,finished=time.time())
    (OUT/'result.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True)


def main():
    global SOURCE,OUT
    if '--grid-response' in sys.argv:
        SOURCE=BASE/'numerical_solution_hessian_gridresponse'
        OUT=BASE/'numerical_solution_hessian_gridresponse_verified'
    if '--validate' in sys.argv:validate();return
    OUT.mkdir(parents=True,exist_ok=True)
    status=dict(status='waiting_for_numerical_producer',host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().split()[21],started=time.time(),max_wait_seconds=21600)
    p=OUT/'coordinator_status.json'
    def save():
        tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    child=None
    def stop(signum,frame):raise KeyboardInterrupt(f'termination signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    try:
        save();deadline=time.time()+21600
        while time.time()<deadline:
            raw=json.loads((SOURCE/'status.json').read_text())
            shards=[SOURCE/f'shard_{i}_status.json' for i in (0,1)]
            if (SOURCE/'shard_plan.json').exists() and not all(p.exists() for p in shards):
                time.sleep(15);continue
            if all(p.exists() for p in shards):
                progress=[json.loads(p.read_text()) for p in shards]
                if any(d['status']=='failed' for d in progress):raise RuntimeError('numerical shard failed')
                if all(d['status']=='columns_verified' for d in progress):
                    raw.update(status='numerical_columns_verified',columns_complete=sorted(
                        set().union(*(set(d['columns']) for d in progress))),shard_owners=[d['owner'] for d in progress],
                        identity=progress[0]['identity'],warm_start_checks=progress[0].get('warm_start_checks',raw.get('warm_start_checks',[])))
                    tmp=SOURCE/'status.shards.tmp';tmp.write_text(json.dumps(raw,indent=2));tmp.replace(SOURCE/'status.json')
                    break
                time.sleep(15);continue
            pid=raw['pid'];proc=Path(f'/proc/{pid}/stat')
            alive=proc.exists() and proc.read_text().split()[21]==raw['process_start']
            if not alive:
                if len(raw.get('columns_complete',[]))!=raw.get('expected_columns'):
                    raise RuntimeError('producer ended without all numerical outputs')
                break
            time.sleep(15)
        else:raise RuntimeError('producer deadline expired')
        while time.time()<deadline:
            slots=subprocess.check_output([GATE,'list','--hosts','lambda1'],text=True)
            if 'lambda1 0' in slots.splitlines():break
            time.sleep(15)
        else:raise RuntimeError('GPU reservation deadline expired')
        argv=[GATE,'run','0','--label','redox-numeric-hessian-postcheck-20261011','--','taskset','-c','0,1',
              'timeout','--signal=TERM','--kill-after=30','1h',sys.executable,'-u',str(Path(__file__).resolve()),'--validate']
        if '--grid-response' in sys.argv:argv.append('--grid-response')
        status.update(status='validating_postprocessing',argv=argv);save()
        with (OUT/'validation.log').open('w') as file:
            child=subprocess.Popen(argv,cwd=ROOT,stdout=file,stderr=subprocess.STDOUT,start_new_session=True)
            status['gate_pid']=child.pid;save();rc=child.wait()
        if rc:raise RuntimeError(f'postcheck failed {rc}; inspect validation.log')
        result=json.loads((OUT/'result.json').read_text())
        status.update(status=result['status'],stationary_minimum=result['stationary_minimum'],finished=time.time());save()
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()


if __name__=='__main__':main()
