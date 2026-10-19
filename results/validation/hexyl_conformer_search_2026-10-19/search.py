"""Bounded CPU-only conformer preparation and host-CDS smoothness diagnostics.

This selects different starting basins; it never admits a native Hessian, changes
the Hamiltonian or replaces a canonical calculation.
"""
import ast
import ctypes
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess
from types import SimpleNamespace
import time

import numpy as np
from ase.io import read
from ase.geometry import get_dihedrals
from pyscf import gto
from rdkit import Chem
from rdkit.Chem import Lipinski

from redox.core.common import ROOT,read_manifest
from redox.core.io import atomic_write_json
from redox.qm.uma import geometry_checks
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

base=Path(__file__).resolve().parent
assert not (base/'plan.json').exists(), 'Preserve prior preparation and attempts'
with allocation_guard():
    resources=host_capacity(2)
    cpus=cpu_allocation(2);os.sched_setaffinity(0,cpus)
assert all(os.environ[k]=='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
gid='exp_ndi_hexyl_pentyloxy'
row=next(r for r in read_manifest() if r['id']==gid and r['state']=='red1')
seeds=sorted((ROOT/'library'/gid/'conformers').glob('conf_*.xyz'));assert len(seeds)==10
previous=ROOT/'results/validation/hexyl_finer_derivatives_2026-10-19'
references={'canonical':ROOT/'calcs/dft'/gid/'red1/opt.xyz',
            'failed_minus':previous/'frozen_minus/opt.xyz',
            'failed_plus':previous/'frozen_plus/opt.xyz'}
package=Path(importlib.util.find_spec('gpu4pyscf').origin).parent
source=package/'solvent/smd.py';library=package/'lib/libsolvent.so'
xtb=ROOT.parent/'miniforge3/envs/redox/bin/xtb'
plan=dict(status='ten_ALPB_conformer_preparations_predeclared',id=gid,state='red1',charge=-1,mult=2,
    source_sha256={str(p):sha(p) for p in [Path(__file__),source,library,xtb,*seeds,*references.values()]},
    method='GFN2-xTB/ALPB(MeCN), --acc 1.0, standard --opt; then installed host SMD CDS only',
    distinct_basin_heuristic=dict(min_heavy_RMSD_A=.5,min_changed_acyclic_torsions=2,torsion_change_deg=60),
    diagnostic_smoothness_heuristic=dict(coarse_step_Bohr=1e-4,comparison_step_Bohr=.0005,
        max_CDS_asymmetry_Ha_Bohr2=1e-6,max_CDS_step_difference_Ha_Bohr2=1e-5),
    selection='At most two lowest ALPB-energy conformers satisfying graph, smoothness and distinct-basin heuristics, also distinct from each other',
    native_admission=False,native_thresholds_changed=False,canonical_outputs_modified=False,
    GPU_used=False,resources=resources,CPU_affinity=cpus,created=time.time())
atomic_write_json(base/'plan.json',plan)
tree=ast.parse(source.read_text());constants={}
for node in tree.body:
    if not isinstance(node,ast.Assign) or not isinstance(node.targets[0],ast.Name):continue
    key=node.targets[0].id
    if key=='hartree2kcal':constants[key]=ast.literal_eval(node.value)
    if key=='solvent_db':constants[key]={'acetonitrile':next(ast.literal_eval(v) for k,v in zip(node.value.keys,node.value.values) if ast.literal_eval(k)=='acetonitrile')}
function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='get_cds_legacy')
namespace=dict(np=np,ctypes=ctypes,libsolvent=ctypes.CDLL(str(library)),**constants)
exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),namespace)
cds=namespace['get_cds_legacy']


def matrices(path,steps):
    atoms=read(path);symbols=atoms.get_chemical_symbols();coords=atoms.positions/.52917721067
    mol=gto.M(atom=list(zip(symbols,coords.tolist())),unit='Bohr',charge=-1,spin=1,basis='sto-3g',verbose=0)
    origin=mol.atom_coords(unit='B').copy();results={}
    for step in steps:
        matrix=np.empty((origin.size,origin.size))
        for col in range(origin.size):
            x=origin.copy();x.flat[col]+=step;mol.set_geom_(x,unit='B')
            plus=cds(SimpleNamespace(mol=mol,solvent='acetonitrile',solvent_descriptors=None))[1].copy().ravel()
            if step==1e-4:x.flat[col]-=2*step  # Match the installed coarse operation sequence.
            else:x=origin.copy();x.flat[col]-=step
            mol.set_geom_(x,unit='B')
            minus=cds(SimpleNamespace(mol=mol,solvent='acetonitrile',solvent_descriptors=None))[1].copy().ravel()
            matrix[col]=(plus-minus)/(2*step)
        results[step]=matrix
    return results


controls=[]
for branch in ('minus','plus'):
    actual=matrices(references['failed_'+branch],[1e-6])[1e-6]
    expected=np.load(previous/'cpu_cds_diagnosis'/f'{branch}_fine_CDS.npy')
    delta=float(np.max(np.abs(actual-expected)));assert delta<1e-10,delta
    controls.append(dict(branch=branch,max_difference_from_preserved_CDS_Ha_Bohr2=delta))
atomic_write_json(base/'host_interface_controls.json',dict(status='preserved_failed_geometry_CDS_controls_reproduced',controls=controls))
template=read(references['canonical']);numbers=template.numbers;heavy=np.flatnonzero(numbers>1)
declared=Chem.MolFromSmiles(row['smiles'])
from redox.build.geometry import smiles_to_xyz_map
mapping=smiles_to_xyz_map(row['smiles'],references['canonical']);assert mapping is not None
# These stored seeds and DFT descendants preserve heavy-first atom order.
assert np.array_equal(heavy,np.arange(len(heavy)))
torsions=[]
for i,j in declared.GetSubstructMatches(Lipinski.RotatableBondSmarts):
    left=sorted(a.GetIdx() for a in declared.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx()!=j)
    right=sorted(a.GetIdx() for a in declared.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx()!=i)
    if left and right:torsions.append(tuple(mapping[k] for k in (left[0],i,j,right[0])))


def compare(path,reference):
    a,b=read(path),read(reference);assert np.array_equal(a.numbers,numbers) and np.array_equal(b.numbers,numbers)
    p,q=a.positions[heavy],b.positions[heavy];p=p-p.mean(axis=0);q=q-q.mean(axis=0)
    u,_,vt=np.linalg.svd(p.T@q);correction=np.eye(3);correction[-1,-1]=np.linalg.det(u@vt)
    rmsd=float(np.sqrt(np.mean(np.sum((p@u@correction@vt-q)**2,axis=1))))
    def angles(atoms):
        x=atoms.positions;return np.array([get_dihedrals((x[b]-x[a])[None,:],(x[c]-x[b])[None,:],(x[d]-x[c])[None,:])[0] for a,b,c,d in torsions])
    delta=np.abs((angles(a)-angles(b)+180)%360-180)
    return dict(heavy_RMSD_A=rmsd,changed_torsions_60deg=int(np.sum(delta>=60)),torsion_differences_deg=delta.tolist())


records=[]
for seed in seeds:
    out=base/seed.stem;out.mkdir()
    (out/'input.xyz').write_bytes(seed.read_bytes())
    argv=[str(xtb),'input.xyz','--gfn','2','--chrg','-1','--uhf','1','--opt','--alpb','acetonitrile','--acc','1.0']
    entry=dict(seed=str(seed),seed_sha256=sha(seed),argv=argv,started=time.time())
    with (out/'xtb.log').open('w') as log:
        process=subprocess.run(argv,cwd=out,stdout=log,stderr=subprocess.STDOUT)
    entry['returncode']=process.returncode
    optimized=out/'xtbopt.xyz';text=(out/'xtb.log').read_text()
    energies=re.findall(r'TOTAL ENERGY\s+([-+0-9.EeDd]+)',text)
    if process.returncode or not optimized.exists() or not energies or 'GEOMETRY OPTIMIZATION CONVERGED' not in text:
        entry.update(status='ALPB_failed_preserved',finished=time.time());records.append(entry)
        atomic_write_json(base/'progress.json',dict(records=records,completed=len(records),expected=10));continue
    energy=float(energies[-1].replace('D','E').replace('d','e'));assert math.isfinite(energy)
    assert not geometry_checks(row['smiles'],optimized,-1,2)
    assert hydrogen_topology_status(row['smiles'],optimized)=='ok'
    entry.update(status='ALPB_geometry_and_CDS_screen_completed',ALPB_energy_Ha=energy,geometry=str(optimized),geometry_sha256=sha(optimized))
    entry['comparisons']={name:compare(optimized,path) for name,path in references.items()}
    arrays=matrices(optimized,[1e-4,.0005]);coarse,fd=arrays[1e-4],arrays[.0005]
    np.save(out/'coarse_CDS_Ha_Bohr2.npy',coarse);np.save(out/'FD_CDS_Ha_Bohr2.npy',fd)
    asym=float(np.max(np.abs(coarse-coarse.T)));difference=float(np.max(np.abs(coarse-fd)))
    entry.update(CDS_asymmetry_Ha_Bohr2=asym,CDS_step_difference_Ha_Bohr2=difference,
        distinct_from_preserved_basins=all(c['heavy_RMSD_A']>=.5 and c['changed_torsions_60deg']>=2 for c in entry['comparisons'].values()),
        CDS_smoothness_heuristic_passed=asym<=1e-6 and difference<=1e-5,finished=time.time())
    records.append(entry);atomic_write_json(base/'progress.json',dict(records=records,completed=len(records),expected=10))
    print(seed.stem,energy,entry['distinct_from_preserved_basins'],entry['CDS_smoothness_heuristic_passed'],flush=True)
selected=[]
for entry in sorted((r for r in records if r.get('distinct_from_preserved_basins') and r.get('CDS_smoothness_heuristic_passed')),key=lambda r:r['ALPB_energy_Ha']):
    if any((c:=compare(entry['geometry'],r['geometry']))['heavy_RMSD_A']<.5 or c['changed_torsions_60deg']<2 for r in selected):continue
    selected.append(entry)
    if len(selected)==2:break
for p,h in plan['source_sha256'].items():assert sha(p)==h,p
atomic_write_json(base/'review.json',dict(status='distinct_seed_preparation_completed_not_native_admission',records=records,selected_for_future_DFT=[r['geometry'] for r in selected],all_selected_are_heuristic_seeds_only=True,canonical_outputs_modified=False,GPU_used=False,DFT_SCF_repeated=False,native_thresholds_changed=False,checked=time.time()))
print('Selected new seeds:',len(selected))
