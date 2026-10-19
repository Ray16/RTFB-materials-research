"""Measure the installed CDS gradient across a specific solvent-sphere tangency."""
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
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

base=Path(__file__).resolve().parent
out=base/'contact_sweep_refined';out.mkdir(exist_ok=False)
with allocation_guard():
    resources=host_capacity(1);cpus=cpu_allocation(1);os.sched_setaffinity(0,cpus)
assert os.environ['OMP_NUM_THREADS']==os.environ['OPENBLAS_NUM_THREADS']=='1'
package=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
source=package/'solvent/smd.py';library=package/'lib/libsolvent.so'
tree=ast.parse(source.read_text());constants={}
for node in tree.body:
    if not isinstance(node,ast.Assign) or not isinstance(node.targets[0],ast.Name):continue
    name=node.targets[0].id
    if name=='hartree2kcal':constants[name]=ast.literal_eval(node.value)
    if name=='solvent_db':constants[name]={'acetonitrile':next(ast.literal_eval(v) for k,v in zip(node.value.keys,node.value.values) if ast.literal_eval(k)=='acetonitrile')}
function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='get_cds_legacy')
namespace=dict(np=np,ctypes=ctypes,libsolvent=ctypes.CDLL(str(library)),**constants)
exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),namespace)
cds=namespace['get_cds_legacy']
previous=ROOT/'results/validation/hexyl_finer_derivatives_2026-10-19'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
records=[]
# The tagged upstream source uses 1.2-A H radii plus a 0.4-A solvent radius.
# Its legacy Angstrom/Bohr constant differs slightly from the pipeline constant.
legacy_to_A=.52917724924;pair_threshold_legacy_A=2*(float(np.float32(1.20))+.4)
for branch in ('minus','plus'):
    path=previous/f'frozen_{branch}/opt.xyz';atoms=read(path);i,j=40,71
    assert atoms[i].symbol==atoms[j].symbol=='H'
    coords=atoms.positions/.52917721067
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),coords.tolist())),unit='Bohr',charge=-1,spin=1,basis='sto-3g',verbose=0)
    origin=mol.atom_coords(unit='B').copy();vector=origin[j]-origin[i];direction=vector/np.linalg.norm(vector)
    target=origin.copy();target[j]=origin[i]+direction*(pair_threshold_legacy_A/legacy_to_A)
    details=[]
    for step_A in (1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10):
        sides=[]
        for sign in (-1,1):
            x=target.copy();x[j]+=sign*direction*(step_A/legacy_to_A);mol.set_geom_(x,unit='B')
            energy,gradient=cds(SimpleNamespace(mol=mol,solvent='acetonitrile',solvent_descriptors=None))
            gradient=gradient.copy();assert np.isfinite(gradient).all() and np.isfinite(energy)
            np.save(out/f'{branch}_{step_A:g}_{sign}_gradient.npy',gradient)
            sides.append(dict(sign=sign,energy_Ha=float(energy),directional_gradient_Ha_Bohr=float(gradient[j]@direction),gradient=gradient))
        jump=sides[1]['gradient']-sides[0]['gradient']
        details.append(dict(step_legacy_A=step_A,step_Bohr=step_A/legacy_to_A,
            energy_difference_Ha=sides[1]['energy_Ha']-sides[0]['energy_Ha'],
            full_gradient_jump_norm_Ha_Bohr=float(np.linalg.norm(jump)),
            directional_gradient_jump_Ha_Bohr=sides[1]['directional_gradient_Ha_Bohr']-sides[0]['directional_gradient_Ha_Bohr'],
            cross_contact_directional_curvature_Ha_Bohr2=float((jump[j]@direction)/(2*step_A/legacy_to_A)),
            sides=[{k:v for k,v in s.items() if k!='gradient'} for s in sides]))
    records.append(dict(branch=branch,geometry_sha256=sha(path),atom_indices_zero_based=[i,j],
        original_distance_pipeline_A=float(np.linalg.norm(vector)*.52917721067),
        original_distance_legacy_A=float(np.linalg.norm(vector)*legacy_to_A),samples=details))
sources=[Path(__file__),source,library,previous/'frozen_minus/opt.xyz',previous/'frozen_plus/opt.xyz']
result=dict(status='installed_host_CDS_precision_resolved_contact_sweep_completed',records=records,
    radius_hypothesis_source='https://github.com/pyscf/gpu4pyscf/blob/v1.8.1/gpu4pyscf/lib/solvent/mnsol.F',
    pair_threshold_legacy_A=pair_threshold_legacy_A,legacy_A_per_Bohr=legacy_to_A,upstream_unsuffixed_real_literal_precision_hypothesis=True,source_binary_equivalence_not_proven=True,
    source_sha256={str(p):sha(p) for p in sources},GPU_used=False,DFT_SCF_repeated=False,
    canonical_values_changed=False,native_admission=False,resources=resources,CPU_affinity=cpus,checked=time.time())
atomic_write_json(out/'report.json',result)
print(json.dumps([{k:r[k] for k in ('branch','original_distance_pipeline_A')}|dict(last_samples=r['samples'][-2:]) for r in records]))
