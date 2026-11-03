"""Reconstruct a failed primary gradient from its actual SCF and wavefunction."""
import argparse
import json
import os
from pathlib import Path
import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.execution_review import freeze_source
from redox.validation.failed_execution import physical_failure
from redox.validation.primary_checkpoints import checkpoint_series
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.route_primary_native import validate_plan


def diagnose(namespace):
    spec=verify(namespace);case=Path(spec['case']);used=[namespace]
    def read(p):used.append(Path(p));return json.loads(Path(p).read_text())
    plan=read(case/'plan.json');validate_plan(plan)
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','failed_stationarity_owners')
    physical_failure(case,read,owners.alive,status='route_primary_native_component_failed_preserved',
        error="RuntimeError('Source is not stationary under the conservative solution gradient; inspect before replacing geometry')",host_GiB=40)
    geometry=case/'state/opt.xyz';used.append(geometry)
    series=checkpoint_series(case/'SCF_evaluations',geometry,
        dict(charge=plan['charge'],mult=plan['mult'],basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,
            conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),read,used,expected_count=1)
    gradients=list((case/'SCF_evaluations').glob('*/*/gradient.npy'))
    if len(gradients)!=1 or list((case/'state/native_thermal').glob('*/hessian_cartesian.npy')):
        raise ValueError('One actual pre-Hessian gradient failure required')
    used.append(gradients[0]);gradient=np.load(gradients[0],allow_pickle=False)
    if gradient.shape!=(39,3) or not np.isfinite(gradient).all() or np.max(abs(gradient))<=.00045:
        raise ValueError('The actual finite gradient must exceed the unchanged stationarity gate')
    for fn in (diagnose,physical_failure,checkpoint_series,validate_plan):used.append(freeze_source(case,Path(fn.__code__.co_filename)))
    verify(namespace)
    return dict(status='route_primary_stationarity_failure_independently_reconstructed',namespace=str(namespace),
        candidate_alias=plan['candidate_alias'],source_id=plan['source_id'],smiles=plan['smiles'],
        charge=plan['charge'],mult=plan['mult'],geometry=str(geometry),geometry_sha256=sha256_file(geometry),
        maximum_gradient_Ha_Bohr=float(abs(gradient).max()),stationarity_limit_Ha_Bohr=.00045,
        actual_SCF_frames_checked=series['actual_SCF_frames_checked'],Hessian_computed=False,
        native_G_qualified=False,redox_potential_qualified=False,material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve prior failed-gradient diagnoses')
    result=diagnose(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','maximum_gradient_Ha_Bohr','actual_SCF_frames_checked')}))
