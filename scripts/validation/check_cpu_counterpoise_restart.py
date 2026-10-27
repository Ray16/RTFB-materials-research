#!/usr/bin/env python3
"""Check CPU SCF checkpoint restart against the reviewed water-dimer control."""
import argparse
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from pyscf import dft,df,lib

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def check(out):
    if out.exists():raise FileExistsError('Preserve existing restart checks')
    report_path=ROOT/'results/reviews/2026-10-26/counterpoise_control/independent_review.json'
    report=json.loads(report_path.read_text());api.verify_sources(report['source_sha256'])
    if report['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Unqualified baseline')
    core=api.load_module('cpu_restart_baseline_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    if core.SETTINGS!=report['settings']:raise ValueError('Protocol differs from reviewed baseline')
    base=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'
    checkpoint=base/'A_own/SCF.h5';reference=json.loads((base/'A_own/energy.json').read_text())
    summary=json.loads((base/'summary.json').read_text());symbols=summary['source_symbols'];coordinates=summary['source_coordinates_A'];a=summary['A_indices']
    for path in (checkpoint,base/'A_own/energy.json',base/'grid_coordinates_AU.npy',base/'grid_weights.npy'):
        if report['source_sha256'].get(str(path))!=api.digest(path):raise ValueError('Baseline source changed')
    mols=core.molecules(symbols,coordinates,a);mol=mols['A_own'];full=mols['AB'];saved=lib.chkfile.load_mol(str(checkpoint))
    if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Checkpoint molecular identity differs')
    coords=np.load(base/'grid_coordinates_AU.npy');weights=np.load(base/'grid_weights.npy')
    out.mkdir(parents=True);mf=dft.RKS(mol).density_fit(only_dfj=False);mf.with_df.auxmol=df.addons.make_auxmol(full,None)
    mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.
    mf.conv_tol=core.SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=core.SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=core.SETTINGS['SCF_max_cycle']
    mf.grids.level=3;mf.nlcgrids.level=3
    for grids in (mf.grids,mf.nlcgrids):grids.coords=coords.copy();grids.weights=weights.copy()
    fixed=api.load_module('cpu_restart_fixed_auxiliary',ROOT/'scripts/validation/run_cpu_counterpoise_reference.py')
    mf.chkfile=str(out/'SCF.h5');start=time.time();initial=mf.from_chk(str(checkpoint),project=False)
    with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):
        energy=float(mf.kernel(dm0=initial))
    auxiliary=df.addons.make_auxmol(full,None)
    if (mf.with_df.auxmol.nao_nr()!=auxiliary.nao_nr()
            or not np.array_equal(mf.with_df.auxmol.atom_charges(),auxiliary.atom_charges())
            or not np.allclose(mf.with_df.auxmol.atom_coords(),auxiliary.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Fixed full-complex auxiliary basis was not retained')
    atomic_write_json(out/'diagnostic.json',dict(energy_Ha=energy,reference_energy_Ha=reference['E_Ha'],
        absolute_error_Ha=abs(energy-reference['E_Ha']),converged=bool(mf.converged),qualified=False))
    if not mf.converged or not math.isfinite(energy) or abs(energy-reference['E_Ha'])>1e-8:
        raise ValueError('Checkpoint restart did not reproduce reviewed energy')
    result=dict(status='CPU_counterpoise_checkpoint_restart_reproduced',role='A_own',energy_Ha=energy,
        reference_energy_Ha=reference['E_Ha'],absolute_error_Ha=abs(energy-reference['E_Ha']),seconds=time.time()-start,
        reference_checkpoint_sha256=api.digest(checkpoint),new_checkpoint_sha256=api.digest(out/'SCF.h5'),
        scientific_scope='Water-dimer CPU control, one term; not the 38-atom pilot or full counterpoise')
    atomic_write_json(out/'summary.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(2);os.sched_setaffinity(0,cpu_allocation(2))
    if any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two CPU threads and capped native pools required')
    print(json.dumps(check(a.out.resolve())))
