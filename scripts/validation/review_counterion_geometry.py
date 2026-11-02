#!/usr/bin/env python3
"""Review a chloride-contact geometry without assigning association G or yield."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import h5py
from ase.io import read as atoms_read
from pyscf.lib import param
from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load

SOURCE_BYTES=Path(__file__).read_bytes()


def review(namespace):
    from redox.validation.execution_review import review_closure,freeze_source,SOURCE_BYTES as closure_source
    from redox.validation.scf_frame import checkpoint_frame,checkpoint_spin,freeze_reviewer
    from redox.validation.protonation import hydrogen_topology_status
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(p):used.append(p);return json.loads(p.read_text())
    plan=read(case/'plan.json')
    if (plan['role']!='primary_counterion_geometry_control' or (plan['charge'],plan['mult'])!=(1,1)
            or Chem.GetFormalCharge(Chem.MolFromSmiles(plan['smiles']))!=1 or plan['parent_formal_core_charge']!=2):
        raise ValueError('Exact oxidized-core chloride microstate required, distinct from a reduced radical')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','chloride_geometry_owner')
    stage=review_closure(case,read,owners.alive,expected_status='counterion_geometry_control_finished_independent_review_pending',
        host_GiB=32,GPU_MiB=31744,worker=sources/'worker.py')
    used += [freeze_source(case,Path(__file__),SOURCE_BYTES),freeze_reviewer(case),
        freeze_source(case,ROOT/'src/redox/validation/execution_review.py',closure_source)]
    result=read(case/'optimizer_result.json');final=read(case/'final_evaluation.json');xyz=case/'optimizer_return.xyz';used.append(xyz)
    if (not result['optimized'] or not result['geometry_optimization_converged'] or not result['optimization_SCF_converged']
            or not final['SCF_converged'] or hydrogen_topology_status(plan['smiles'],xyz)!='ok'):
        raise ValueError('Actual optimizer convergence and exact final ionic/covalent graph required')
    atoms=atoms_read(xyz);gradient=np.asarray(final['gradient_Ha_Bohr'])
    if gradient.shape!=(len(atoms),3) or not np.isfinite(gradient).all() or abs(gradient).max()>.00045:
        raise ValueError('Fresh optimizer-return stationarity check failed')
    paths=sorted((case/'SCF_evaluations').glob('*/*/SCF.h5'))
    if not paths:raise ValueError('Every accepted geometry needs actual final-SCF evidence')
    checkpoint=paths[-1];settings=read(checkpoint.with_name('settings.json'));solved=read(checkpoint.with_name('SCF_result.json'))
    if (settings['xc']!='r2scan' or settings['disp']!='d4' or settings['conv_tol']!=1e-12
            or settings['conv_tol_grad']!=1e-8 or settings['max_cycle']!=200 or settings['grid_level']!=3
            or settings['solvent_class']!='SMD' or not solved['converged'] or abs(solved['energy_Ha']-final['E_Ha'])>1e-10):
        raise ValueError('Actual final primary-model SCF settings/convergence differ')
    frame=checkpoint_frame(checkpoint,charge=1,mult=1,basis='def2-svp',coordinates_Bohr=atoms.positions/param.BOHR,symbols=atoms.get_chemical_symbols())
    with h5py.File(checkpoint,'r') as stream:actual_energy=float(stream['scf/e_tot'][()])
    if (not all(math.isfinite(value) for value in (actual_energy,solved['energy_Ha'],final['E_Ha']))
            or abs(actual_energy-solved['energy_Ha'])>1e-10 or abs(actual_energy-final['E_Ha'])>1e-10):
        raise ValueError('Actual final SCF checkpoint energy differs from geometry records')
    spin=checkpoint_spin(checkpoint);stored=np.load(checkpoint.with_name('gradient.npy'),allow_pickle=False)
    if not np.allclose(stored,gradient,rtol=0,atol=1e-12):raise ValueError('Actual final gradient differs')
    used.extend([checkpoint,checkpoint.with_name('gradient.npy')])
    verify(namespace)
    host=read(case/'host_guard/status.json')
    return dict(status='primary_counterion_geometry_independently_checked',namespace=str(namespace),key=plan['source_id'],
        smiles=plan['smiles'],parent_smiles=plan['parent_smiles'],total_charge=1,multiplicity=1,
        formal_core_charge=2,geometry_sha256=sha256_file(xyz),nao=plan['nao'],method=plan['method'],
        actual_final_SCF_energy_Ha=actual_energy,
        maximum_gradient_Ha_Bohr=float(abs(gradient).max()),actual_frame=frame,actual_occupied_spin=spin,
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes']),
        native_G_qualified=False,association_free_energy_computed=False,experimental_speciation_inferred=False,material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier counterion geometry reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','key','maximum_gradient_Ha_Bohr','peak_GPU_MiB')}))
