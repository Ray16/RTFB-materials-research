#!/usr/bin/env python
"""Use the actual solvent Hessian operator and complete first-atom H1.

The inherited checkpoint is the pure-solute make_h1 intermediate. The PCM
nuclear derivative must be added, and equilibrium solvent response enabled,
exactly as in the original full SMD Hessian driver. No SCF is repeated.
"""
import inspect
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def solvated_source(body, reorder):
    body=reorder(body)
    old="            result, e1 = local['solve_mo1'](mf,mf.mo_energy,mf.mo_coeff,mf.mo_occ,H1,max_cycle=100)"
    new="""            from gpu4pyscf.solvent.hessian import pcm
            mf.with_solvent.equilibrium_solvation=True
            hessian=mf.Hessian(); hessian.grid_response=True
            pcm_H1=pcm.analytical_grad_vmat(mf.with_solvent,mf.make_rdm1(),mf.mo_coeff,mf.mo_occ,atmlst=[0])
            H1=np.array(H1[:1])+cp.asnumpy(pcm_H1)
            np.save(raw/'first_atom_complete_H1.npy',H1)
            del pcm_H1
            dft._free_gpu_pool()
            result, e1 = local['solve_mo1'](mf,mf.mo_energy,mf.mo_coeff,mf.mo_occ,H1,
                fx=hessian.gen_vind(mf.mo_coeff,mf.mo_occ),max_cycle=100)"""
    if body.count(old)!=1:raise ValueError('Fixed-state solve template changed')
    body=body.replace(old,new,1)
    old='rhf.gen_vind(mf,mf.mo_coeff,mf.mo_occ)'
    if body.count(old)!=1:raise ValueError('Independent response template changed')
    return body.replace(old,'hessian.gen_vind(mf.mo_coeff,mf.mo_occ)',1)


def main():
    source=ROOT/'scripts/validation/probe_fixed_orbital_response.py'
    ordering=ROOT/'scripts/validation/probe_reordered_response.py'
    parent=api.load_module('complete_solvated_response_probe',source)
    reorder=api.load_module('complete_response_original_build_order',ordering)
    if sys.argv[1]=='run':
        namespace=dict(vars(parent));local={}
        exec(compile(solvated_source(inspect.getsource(parent.run),reorder.reordered_source),
            str(Path(__file__).resolve()),'exec'),namespace,local)
        parent.run=local['run']
    parent.main()
    if sys.argv[1]=='prepare':
        campaign=Path(sys.argv[sys.argv.index('--campaign')+1]).resolve()
        plan=api.read_plan(campaign)
        plan['execution_driver']=str(Path(__file__).resolve())
        plan['producer_sha256'].update({str(p):api.digest(p) for p in (Path(__file__).resolve(),ordering)})
        plan.update(construction_order='PCM before DF selection',equilibrium_solvation=True,
            PCM_nuclear_H1_required=True,only_atom_zero_solved=True)
        atomic_write_json(campaign/'plan.json',plan)


if __name__=='__main__':
    main()
