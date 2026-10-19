#!/usr/bin/env python
"""Independent grid-response/energy-derivative checks at the worst asymmetric Hessian pair."""
import hashlib,json,os,socket,time
from pathlib import Path
import numpy as np
from ase.io import read
from pyscf import gto
from redox.core.common import ROOT,DFT
from redox.qm.dft import _build_smd_mf,_free_gpu_pool

OUT=ROOT/'results/validation/physical_thermal_spin_2026-10-10/numerical_grid_checks'
SOURCE=OUT.parent/'numerical_solution_hessian_aq_red2_tight'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    atoms=read(DFT/'aq_benzylamino/red2/opt.xyz');symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
    reference=np.stack([np.load(SOURCE/f'column_{i:03d}.npy').ravel() for i in range(3*len(atoms))],axis=1)
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),started=time.time(),
                affinity=sorted(os.sched_getaffinity(0)),gpu=os.getenv('CUDA_VISIBLE_DEVICES'),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),models=[])
    def save():
        p=OUT/'status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    def solve(pos,level):
        mol=gto.M(atom=list(zip(symbols,pos.tolist())),unit='Bohr',basis='def2-svpd',charge=-2,spin=0,verbose=0)
        mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4');mf.grids.level=level
        mf.conv_tol=1e-12;mf.conv_tol_grad=1e-8;mf.max_cycle=200;e=float(mf.kernel())
        if not mf.converged:raise RuntimeError('SCF not converged')
        grad=mf.nuc_grad_method();grad.grid_response=True;raw=grad.kernel()
        g=np.asarray(raw.get() if hasattr(raw,'get') else raw)
        del mf,grad;_free_gpu_pool();return e,g
    save()
    try:
        for level in (3,4):
            e0,g0=solve(coords,level);np.save(OUT/f'base_gradient_level{level}.npy',g0)
            item=dict(grid_level=level,grid_response=True,base_e_Ha=e0,base_gradient_max=float(abs(g0).max()),columns=[])
            status['models'].append(item);save()
            pair={}
            for column in (70,73):
                ia,axis=divmod(column,3);es=[];gs=[]
                for sign in (-1,1):
                    x=coords.copy();x[ia,axis]+=sign*5e-4;e,g=solve(x,level);es.append(e);gs.append(g)
                vector=((gs[1]-gs[0])/1e-3).ravel();pair[column]=vector
                np.save(OUT/f'column{column}_level{level}.npy',vector)
                item['columns'].append(dict(column=column,relative_to_original_column=float(
                    np.linalg.norm(vector-reference[:,column])/np.linalg.norm(vector)),
                    energy_derivative_error_Ha_Bohr=float((es[1]-es[0])/1e-3-g0[ia,axis])))
                save();print('GRID_COLUMN_COMPLETE',level,column,item['columns'][-1],flush=True)
            item['worst_pair_asymmetry_Ha_Bohr2']=float(abs(pair[70][73]-pair[73][70]));save()
        status.update(status='grid_sensitivity_completed',finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
