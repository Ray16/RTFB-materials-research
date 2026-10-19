"""Host-only CDS derivative sweep; no SCF, GPU context or production admission."""
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
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

base = Path(__file__).resolve().parent
out = base / 'cpu_cds_diagnosis'
out.mkdir(exist_ok=False)
with allocation_guard():
    resources = host_capacity(1)
    cpus = cpu_allocation(1)
    os.sched_setaffinity(0, cpus)
assert os.environ['OMP_NUM_THREADS'] == os.environ['OPENBLAS_NUM_THREADS'] == '1'
# Extract only the installed host Fortran interface. Importing gpu4pyscf is
# unnecessary, and no CUDA runtime API is invoked in this diagnostic.
package = Path(importlib.util.find_spec('gpu4pyscf').origin).parent
source = package / 'solvent/smd.py'
library = package / 'lib/libsolvent.so'
tree = ast.parse(source.read_text())
assignments = {}
for node in tree.body:
    if not isinstance(node, ast.Assign) or not isinstance(node.targets[0], ast.Name):
        continue
    name = node.targets[0].id
    if name == 'hartree2kcal':
        assignments[name] = ast.literal_eval(node.value)
    elif name == 'solvent_db':
        assignments[name] = {'acetonitrile': next(ast.literal_eval(value)
                             for key,value in zip(node.value.keys,node.value.values)
                             if ast.literal_eval(key) == 'acetonitrile')}
function = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                and node.name == 'get_cds_legacy')
namespace = dict(np=np, ctypes=ctypes, libsolvent=ctypes.CDLL(str(library)), **assignments)
exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), namespace)
cds = namespace['get_cds_legacy']
results = []
for branch in ('minus', 'plus'):
    xyz = base / f'frozen_{branch}/opt.xyz'
    atoms = read(xyz)
    coordinates = atoms.positions / .52917721067
    mol = gto.M(atom=list(zip(atoms.get_chemical_symbols(), coordinates.tolist())),
                unit='Bohr', charge=-1, spin=1, basis='sto-3g', verbose=0)
    # Only atom charges and coordinates enter the host CDS interface; the basis
    # constructs a Mole object but no integrals, electronic density or SCF.
    coords = mol.atom_coords(unit='B').copy()

    def evaluate(x):
        mol.set_geom_(x, unit='B')
        return cds(SimpleNamespace(mol=mol, solvent='acetonitrile', solvent_descriptors=None))

    def hessian(step):
        matrix = np.empty((coords.size, coords.size))
        for column in range(coords.size):
            x = coords.copy(); x.flat[column] += step
            gp = evaluate(x)[1].copy().ravel()
            x = coords.copy(); x.flat[column] -= step
            gm = evaluate(x)[1].copy().ravel()
            matrix[column] = (gp - gm) / (2 * step)
        return matrix

    coarse = hessian(1e-4)
    fine = hessian(1e-6)
    np.save(out/f'{branch}_coarse_CDS.npy',coarse)
    np.save(out/f'{branch}_fine_CDS.npy',fine)
    i, j = map(int, np.unravel_index(np.argmax(abs(coarse-coarse.T)), coarse.shape))
    details = dict(branch=branch, geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                   worst_coarse_indices=[i,j], coarse_maximum_asymmetry=float(abs(coarse-coarse.T).max()),
                   fine_maximum_asymmetry=float(abs(fine-fine.T).max()), columns=[])
    previous = base / branch / 'coarse_CDS_hessian.npy'
    if previous.exists():
        delta = float(abs(coarse-np.load(previous)).max())
        details['maximum_difference_from_preserved_GPU_wrapper_CDS'] = delta
        atomic_write_json(out/f'{branch}_control.json',dict(maximum_difference=delta,
            reproduced=delta<1e-9,coarse_asymmetry=details['coarse_maximum_asymmetry']))
        assert delta < 1e-9, f'Host interface did not reproduce preserved CDS: {delta}'
    for column in sorted({i,j}):
        for step in (1e-6, 1e-5, 1e-4, .0005, .001):
            x = coords.copy(); x.flat[column] += step
            ep,gp = evaluate(x); gp = gp.copy().ravel()
            x = coords.copy(); x.flat[column] -= step
            em,gm = evaluate(x); gm = gm.copy().ravel()
            numerical = (gp-gm)/(2*step)
            details['columns'].append(dict(column=column, step_Bohr=step,
                 norm_Ha_Bohr2=float(np.linalg.norm(numerical)),
                 difference_from_fine_column_Ha_Bohr2=float(np.linalg.norm(numerical-fine[column])),
                 derivative_energy_Ha_Bohr=float((ep-em)/(2*step))))
    np.save(out/f'{branch}_coarse_CDS.npy',coarse)
    np.save(out/f'{branch}_fine_CDS.npy',fine)
    results.append(details)
report = dict(status='host_CDS_diagnostic_completed', completed=time.time(),
              resources=resources, cpu_affinity=cpus, branches=results,
              GPU_used=False, quantum_SCF_repeated=False, production_admission=False,
              producer_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (Path(__file__),source,library)},
              interpretation='Step-dependent host CDS curvature is diagnostic only; existing full quantum failure gates remain unchanged.')
atomic_write_json(out/'report.json',report)
print(json.dumps(results,indent=2))
