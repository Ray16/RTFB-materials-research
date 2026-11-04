"""Read-only as-used diagnosis of the saved failed neutral-minus Hessian."""
import json
import os
from pathlib import Path

from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation

with allocation_guard():
    host_capacity(1, minimum_ram_gib=4)
    os.sched_setaffinity(0, cpu_allocation(1))

import numpy as np
from ase.io import read as atoms_read
from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.workflows.molecular_native_namespace import verify, load
from redox.workflows.route_primary_native import validate_plan
from redox.validation.failed_execution import physical_failure
from redox.validation.primary_checkpoints import checkpoint_series

namespace=Path('/nfs/lambda_stor_01/homes/rzhu/redox_candidate_continuations/2026-11-02_phenothiazine_refined_native/phenothiazine_neu_mode_minus/namespace.json')
spec=verify(namespace); case=Path(spec['case']); used=[namespace,Path(__file__)]
def read(path):
    used.append(Path(path))
    return json.loads(Path(path).read_text())

plan=read(case/'plan.json'); validate_plan(plan)
owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','saved_failed_hessian_owners')
physical_failure(case,read,owners.alive,status='route_primary_native_component_failed_preserved',
    error="RuntimeError('Native Hessian symmetry/finite-value failure')",host_GiB=40)
paths=list((case/'state/native_thermal').glob('*/hessian_cartesian.npy'))
assert len(paths)==1
path=paths[0]; used.append(path); matrix=np.load(path,allow_pickle=False)
assert matrix.shape==(117,117) and matrix.dtype==np.float64 and np.isfinite(matrix).all()
asym=float(abs(matrix-matrix.T).max()); assert asym>1e-6
diag=read(path.with_name('hessian_diagnostics.json'))
assert diag['maximum_asymmetry']==asym and diag['all_finite']
xyz=case/'state/opt.xyz'; used.append(xyz); atoms=atoms_read(xyz)
series=checkpoint_series(case/'SCF_evaluations',xyz,
    dict(charge=0,mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,
         conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),read,used,expected_count=1)
i,j=np.unravel_index(abs(matrix-matrix.T).argmax(),matrix.shape)
result=dict(status='failed_refined_native_Hessian_symmetry_independently_reconstructed',
    namespace=str(namespace),key=plan['source_id'],maximum_asymmetry_Ha_Bohr2=asym,
    unchanged_symmetry_limit_Ha_Bohr2=1e-6,all_finite=True,shape=list(matrix.shape),
    maximum_pair=dict(first_atom=int(i//3),first_axis=int(i%3),first_element=atoms[int(i//3)].symbol,
                      second_atom=int(j//3),second_axis=int(j%3),second_element=atoms[int(j//3)].symbol),
    actual_center_SCF_frames_checked=series['actual_SCF_frames_checked'],
    sampled_derivatives_completed=False,stationary_minimum_qualified=False,native_G_qualified=False,
    redox_potential_qualified=False,physical_thresholds_relaxed=False,scientific_retry_performed=False,
    source_sha256={str(p):sha256_file(p) for p in used})
output=Path(__file__).with_name('neutral_minus_Hessian_symmetry_failure_02.json')
if output.exists(): raise FileExistsError('Preserve earlier diagnoses')
output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','maximum_asymmetry_Ha_Bohr2','native_G_qualified')}))
