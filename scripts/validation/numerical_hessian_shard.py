#!/usr/bin/env python
"""Two disjoint resumable column shards of the already-validated numerical Hessian."""
import argparse,ast,hashlib,json,os,socket,textwrap,time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from redox.core.common import ROOT,DFT
from redox.qm.dft import _build_smd_mf,_free_gpu_pool

BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
OUT=BASE/'numerical_solution_hessian_aq_red2_tight'


def main():
    global OUT
    ap=argparse.ArgumentParser();ap.add_argument('--shard',type=int,choices=[0,1],required=True)
    ap.add_argument('--grid-response',action='store_true');a=ap.parse_args()
    prior=json.loads((OUT/'status_before_sharding.json').read_text());identity=prior['identity']
    source=BASE/'numerical_hessian_tight_runtime_source.py';text=source.read_text()
    if hashlib.sha256(source.read_bytes()).hexdigest()!=identity['code_sha256']:
        raise RuntimeError('preserved scientific producer does not match original column identity')
    original=ast.parse(text);mainfn=next(n for n in original.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    gradientfn=next(n for n in mainfn.body if isinstance(n,ast.FunctionDef) and n.name=='gradient')
    gradient_source=textwrap.dedent(ast.get_source_segment(text,gradientfn))
    if a.grid_response:
        OUT=BASE/'numerical_solution_hessian_gridresponse';OUT.mkdir(parents=True,exist_ok=True)
        old="raw=mf.nuc_grad_method().kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw)"
        new="grad=mf.nuc_grad_method();grad.grid_response=True;raw=grad.kernel();g=np.asarray(raw.get() if hasattr(raw,'get') else raw)"
        if old not in gradient_source:raise RuntimeError('preserved gradient statement changed')
        gradient_source=gradient_source.replace(old,new)
        identity=dict(identity,grid_response=True,grid_level=3,
                      gradient_algorithm_sha256=hashlib.sha256(gradient_source.encode()).hexdigest())
    factory_source='def factory(symbols):\n    baseline_dm=None\n'+textwrap.indent(gradient_source,'    ')+'\n    return gradient\n'
    namespace=dict(np=np,gto=gto,_build_smd_mf=_build_smd_mf,_free_gpu_pool=_free_gpu_pool)
    exec(compile(factory_source,'preserved_gradient_producer','exec'),namespace)
    geometry=DFT/'aq_benzylamino/red2/opt.xyz'
    if hashlib.sha256(geometry.read_bytes()).hexdigest()!=identity['geometry_sha256']:
        raise RuntimeError('source geometry changed')
    atoms=read(geometry);coords=atoms.positions/.52917721067
    gradient=namespace['factory'](atoms.get_chemical_symbols());base_e,_,_=gradient(coords)
    owner=dict(host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().split()[21],
               cpus=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'))
    status=dict(status='running',owner=owner,shard=a.shard,started=time.time(),identity=identity,
                actual_worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                unchanged_gradient_ast_sha256=hashlib.sha256(ast.dump(gradientfn,include_attributes=False).encode()).hexdigest(),
                columns=[],base_e_Ha=base_e)
    report=OUT/f'shard_{a.shard}_status.json'
    def save():
        tmp=report.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(report)
    save()
    try:
        if a.grid_response:
            pilot=[]
            for column in (70,73):
                ia,axis=divmod(column,3);gs=[]
                for sign in (-1,1):
                    x=coords.copy();x[ia,axis]+=sign*identity['step_Bohr'];_,g,_=gradient(x);gs.append(g)
                numeric=((gs[1]-gs[0])/(2*identity['step_Bohr'])).ravel()
                cold=np.load(BASE/'numerical_grid_checks'/f'column{column}_level3.npy')
                error=float(np.linalg.norm(numeric-cold)/np.linalg.norm(cold))
                if error>1e-3:raise RuntimeError('full-response warm/cold derivative check failed')
                pilot.append(dict(atom=ia,axis=axis,relative_warm_cold_error=error))
            status['warm_start_checks']=pilot;save()
        for column in range(a.shard,3*len(atoms),2):
            path=OUT/f'column_{column:03d}.json';data=path.with_suffix('.npy')
            if path.exists() and data.exists():
                old=json.loads(path.read_text());array=np.load(data)
                if old['identity']!=identity or not np.isfinite(array).all():raise RuntimeError('invalid existing column')
                status['columns'].append(column);save();continue
            claim=OUT/f'column_{column:03d}.claim';claim.mkdir()
            (claim/'owner.json').write_text(json.dumps(owner))
            try:
                ia,axis=divmod(column,3);gs=[];energies=[];started=time.time()
                for sign in (-1,1):
                    pos=coords.copy();pos[ia,axis]+=sign*identity['step_Bohr']
                    energy,g,_=gradient(pos);gs.append(g);energies.append(energy)
                result=(gs[1]-gs[0])/(2*identity['step_Bohr'])
                temporary=data.with_suffix(f'.{os.getpid()}.tmp')
                with temporary.open('wb') as f:np.save(f,result)
                temporary.replace(data)
                rec=dict(identity=identity,column=column,energies_Ha=energies,seconds=time.time()-started,
                         actual_worker_sha256=status['actual_worker_sha256'],source_experiment_reused=True)
                tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(rec,indent=2));tmp.replace(path)
                print('COLUMN_COMPLETE',a.shard,column+1,'seconds',time.time()-started,flush=True)
                status['columns'].append(column);save()
            finally:
                (claim/'owner.json').unlink(missing_ok=True);claim.rmdir()
        status.update(status='columns_verified',finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
