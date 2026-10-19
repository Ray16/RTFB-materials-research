"""Attribute a preserved full-gradient failure to its separately evaluated host CDS term."""
import ast
import ctypes
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
from ase.io import read
from pyscf import gto
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

base = Path(__file__).resolve().parent
out = base / 'cpu_cds_diagnosis'
with allocation_guard():
    resources = host_capacity(1)
    cpus = cpu_allocation(1)
    os.sched_setaffinity(0, cpus)
assert os.environ['OMP_NUM_THREADS'] == os.environ['OPENBLAS_NUM_THREADS'] == '1'
control = json.loads((out / 'report.json').read_text())
assert control['branches'][0]['maximum_difference_from_preserved_GPU_wrapper_CDS'] == 0
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest() == h
           for p,h in control['producer_sha256'].items())
package = Path(importlib.util.find_spec('gpu4pyscf').origin).parent
source = package / 'solvent/smd.py'
tree = ast.parse(source.read_text())
values = {}
for node in tree.body:
    if not isinstance(node,ast.Assign) or not isinstance(node.targets[0],ast.Name):
        continue
    if node.targets[0].id == 'hartree2kcal':
        values['hartree2kcal'] = ast.literal_eval(node.value)
    elif node.targets[0].id == 'solvent_db':
        values['solvent_db'] = {'acetonitrile': next(ast.literal_eval(v)
            for k,v in zip(node.value.keys,node.value.values) if ast.literal_eval(k)=='acetonitrile')}
function = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='get_cds_legacy')
namespace = dict(np=np,ctypes=ctypes,libsolvent=ctypes.CDLL(str(package/'lib/libsolvent.so')),**values)
exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),namespace)
atoms = read(base/'frozen_minus/opt.xyz')
coords = atoms.positions/.52917721067
mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(),coords.tolist())),unit='Bohr',
            charge=-1,spin=1,basis='sto-3g',verbose=0)
column,step = 120,.0005
gradients = []
for sign in (-1,1):
    x=coords.copy();x.flat[column]+=sign*step;mol.set_geom_(x,unit='B')
    gradients.append(namespace['get_cds_legacy'](SimpleNamespace(mol=mol,solvent='acetonitrile',
                                                               solvent_descriptors=None))[1].ravel().copy())
cds_fd = (gradients[1]-gradients[0])/(2*step)
cds_hessian = np.load(out/'minus_fine_CDS.npy')
full_hessian_path = next((base/'minus/work/native_thermal').glob('*/hessian_cartesian.npy'))
full_hessian = np.load(full_hessian_path)
full_gradients = [np.load(base/'minus'/f'gradient_{column}_{step}_{sign}.npy') for sign in (-1,1)]
full_fd = (full_gradients[1]-full_gradients[0])/(2*step)
full_error = full_fd-full_hessian[column]
cds_error = cds_fd-cds_hessian[column]
remaining = full_error-cds_error
relative = float(np.linalg.norm(full_error)/np.linalg.norm(full_fd))
expected = json.loads((base/'minus/status.json').read_text())['worst_contact_checks'][0]['relative_error']
assert abs(relative-expected)<1e-12
report=dict(status='preserved_full_residual_attributed',column=column,step_Bohr=step,
    original_relative_error=relative,full_error_norm_Ha_Bohr2=float(np.linalg.norm(full_error)),
    CDS_error_norm_Ha_Bohr2=float(np.linalg.norm(cds_error)),
    residual_after_subtracting_CDS_norm_Ha_Bohr2=float(np.linalg.norm(remaining)),
    residual_after_subtracting_CDS_relative_to_full_FD=float(np.linalg.norm(remaining)/np.linalg.norm(full_fd)),
    CPU_affinity=cpus,resources=resources,GPU_used=False,quantum_SCF_repeated=False,
    thresholds_changed=False,production_admission=False,completed=time.time(),
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
      (Path(__file__),out/'report.json',full_hessian_path,base/'minus/status.json',
       base/'minus'/f'gradient_{column}_{step}_-1.npy',base/'minus'/f'gradient_{column}_{step}_1.npy')})
np.savez(out/'minus_residual_vectors.npz',full_error=full_error,CDS_error=cds_error,remainder=remaining)
atomic_write_json(out/'residual_attribution.json',report)
print(json.dumps(report,indent=2))
