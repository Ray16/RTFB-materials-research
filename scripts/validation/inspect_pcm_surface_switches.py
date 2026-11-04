#!/usr/bin/env python3
"""Inspect SWIG point-retention changes at saved native FD geometries on CPU.

This mirrors the installed GPU surface geometry formula in NumPy. It does not
evaluate electrostatic energies/forces and cannot establish their error cause.
"""
import argparse
import ast
import importlib.metadata
import json
import os
from pathlib import Path
import shutil

import numpy as np
from scipy.spatial.distance import cdist

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.archived_native_namespace import verify_spec


def retention_cutoff(path):
    tree=ast.parse(path.read_text());function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='gen_surface')
    comparisons=[n.value for n in ast.walk(function) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='idx' for t in n.targets)]
    if (len(comparisons)!=1 or not isinstance(comparisons[0],ast.Compare) or len(comparisons[0].ops)!=1
            or not isinstance(comparisons[0].ops[0],ast.Gt) or ast.unparse(comparisons[0].left)!='w * swf'
            or len(comparisons[0].comparators)!=1 or not isinstance(comparisons[0].comparators[0],ast.Constant)):
        raise ValueError('Unrecognized installed point-retention formula')
    return float(comparisons[0].comparators[0].value)


def weights(coords,numbers,rad,ng):
    from pyscf.dft import gen_grid
    from pyscf.solvent.pcm import switch_h
    unit=gen_grid.MakeAngularGrid(ng);radii=rad[numbers];width=radii*np.sqrt(14./ng)
    alpha=.5+radii/width-np.sqrt((radii/width)**2-1./28);inner=radii-alpha*width
    result=[]
    for atom in range(len(coords)):
        points=radii[atom]*unit[:,:3]+coords[atom]
        dij=(cdist(points,coords)-inner)/width;dij[:,atom]=1.;dij[dij<1e-8]=0.
        result.append(unit[:,3]*4*np.pi*np.prod(switch_h(dij),axis=1))
    return np.asarray(result)


def inspect(namespace,output):
    if output.exists():raise FileExistsError('Preserve prior surface diagnostics')
    spec=verify_spec(namespace);reader=NamespaceSources(spec);case=Path(spec['attempt'])
    plan=json.loads((case/'plan.json').read_text());raw=reader.path(spec['original_raw_root'])
    if plan['experimental_model']['electrostatic_surface_points']!=590:raise ValueError('Expected the declared 590-point electrostatic surface')
    paths=sorted((raw/'SCF_evaluations').glob('*/evaluation.json'))
    if not paths:raise ValueError('Saved native SCF geometries required')
    from pyscf.solvent import smd
    from pyscf.data.elements import charge
    installed=Path(importlib.metadata.distribution('gpu4pyscf-cuda12x').locate_file('gpu4pyscf/solvent/pcm.py'))
    cutoff=retention_cutoff(installed)
    if cutoff!=1e-12:raise ValueError('The diagnosed GPU cutoff differs from the expected installed build')
    symbols=[r.split()[0] for r in (raw/'geometry.xyz').read_text().splitlines()[2:]]
    numbers=np.array([charge(s) for s in symbols]);rad=smd.smd_radii(smd.solvent_db['acetonitrile'][2])
    output.mkdir(parents=True);source=output/'producer.py';shutil.copyfile(__file__,source)
    all_weights=[];masks=[]
    for path in paths:
        record=json.loads(path.read_text());coords=np.asarray(record['coordinates_Bohr'])
        if coords.shape!=(len(symbols),3) or not record['SCF_converged']:raise ValueError('Saved coordinate/state coverage differs')
        w=weights(coords,numbers,rad,590);all_weights.append(w);masks.append(w>cutoff)
    array=output/'surface_weights.npz';np.savez(array,weights=np.asarray(all_weights),retained=np.asarray(masks))
    rows=[]
    for index,path in enumerate(paths):
        changed=np.argwhere(masks[index]!=masks[0])
        rows.append(dict(point=path.parent.name,retained_points=int(masks[index].sum()),
            changed_from_base=[dict(atom=int(a),element=symbols[a],sphere_point=int(p),
                base_weight=float(all_weights[0][a,p]),weight=float(all_weights[index][a,p])) for a,p in changed]))
    sources={str(p):sha256_file(p) for p in [namespace,case/'plan.json',installed,Path(smd.__file__),source,array,*paths]}
    result=dict(status='CPU_mirror_of_installed_SMD_electrostatic_surface_retention',key=spec['key'],cutoff=cutoff,
        points_per_sphere=590,records=rows,source_sha256=sources,CPU_affinity=sorted(os.sched_getaffinity(0)),
        GPU_surface_equality_verified=False,electrostatic_forces_evaluated=False,error_cause_established=False,native_G_qualified=False)
    atomic_write_json(output/'summary.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=inspect(a.spec.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],records=result['records'])))
