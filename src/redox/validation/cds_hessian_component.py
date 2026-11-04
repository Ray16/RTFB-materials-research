"""CPU-only step convergence of the unchanged legacy SMD CDS component.

Loads only the installed pure-CPU function and its Fortran library, without
importing CuPy or initializing CUDA. Component reassembly is diagnostic only;
it cannot qualify a full Hessian, minimum or nuclear free energy.
"""
import argparse
import ast
import ctypes
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.runtime_paths import gpu4pyscf_source_root
from redox.qm.cds_finite import gradient_hessian, energy_gradient
from redox.workflows.molecular_native_namespace import verify


def cpu_legacy():
    root=gpu4pyscf_source_root();source=root/'solvent/smd.py';library=root/'lib/libsolvent.so'
    content=source.read_text();tree=ast.parse(content)
    assignments={n.targets[0].id:n.value for n in tree.body if isinstance(n,ast.Assign)
                 and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)}
    # Other solvents can contain arithmetic expressions; bind only the exact
    # literal MeCN descriptor vector required by this declared component study.
    database=assignments['solvent_db']
    descriptors={'acetonitrile':ast.literal_eval(next(value for key,value in zip(database.keys,database.values)
        if isinstance(key,ast.Constant) and key.value=='acetonitrile'))}
    conversion=ast.literal_eval(assignments['hartree2kcal'])
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='get_cds_legacy')
    namespace=dict(np=np,ctypes=ctypes,solvent_db=descriptors,hartree2kcal=conversion,
                   libsolvent=ctypes.CDLL(str(library)))
    exec(compile(ast.Module(body=[function],type_ignores=[]),str(source),'exec'),namespace)
    if 'cupy' in sys.modules:raise ValueError('This component study must remain CPU-only')
    return namespace['get_cds_legacy'],{str(p):sha256_file(p) for p in (source,library)},descriptors['acetonitrile']


def analyze(namespace_path, output):
    if output.exists() or output.with_suffix('.npz').exists():raise FileExistsError('Preserve prior component studies')
    from pyscf import gto
    spec=verify(namespace_path);case=Path(spec['case']);plan_path=case/'plan.json'
    plan=json.loads(plan_path.read_text());stage=case/'status.json'
    status=json.loads(stage.read_text())
    if (not plan.get('refined_route_native_component') or plan['charge']!=0 or plan['mult']!=1
            or status['status']!='route_primary_native_component_failed_preserved'
            or status['error']!="RuntimeError('Native Hessian symmetry/finite-value failure')"):
        raise ValueError('Only the exact preserved failed refined neutral Hessian is scoped here')
    settings=list((case/'SCF_evaluations').glob('*/*/settings.json'))
    matrices=list((case/'state/native_thermal').glob('*/hessian_cartesian.npy'))
    if len(settings)!=1 or len(matrices)!=1:raise ValueError('One actual center and complete saved matrix required')
    point=json.loads(settings[0].read_text());coordinates=np.asarray(point['coordinates_Bohr'],dtype=float)
    from ase.io import read
    geometry=case/'state/opt.xyz';atoms=read(geometry)
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),coordinates.tolist())),unit='Bohr',
              basis='def2-svp',charge=0,spin=0,verbose=0)
    if mol.natm!=39 or mol.nao!=397:raise ValueError('Exact checked phenothiazine component dimension required')
    legacy,sources,descriptors=cpu_legacy()
    state=SimpleNamespace(mol=mol,solvent='acetonitrile',solvent_descriptors=None)
    def evaluate(coords):
        mol.set_geom_(coords,unit='Bohr')
        value,gradient=legacy(state)
        if not np.isfinite(value) or not np.isfinite(gradient).all():raise ValueError('Finite CDS energy/gradient required')
        return value,np.asarray(gradient).copy()
    baseline=np.load(matrices[0],allow_pickle=False)
    if baseline.shape!=(117,117) or not np.isfinite(baseline).all():raise ValueError('Finite saved full matrix required')
    energy,gradient=evaluate(coordinates)
    parts={};metrics=[]
    for step in (1e-4,1e-5,1e-6,1e-7):
        tensor=gradient_hessian(coordinates,lambda p:evaluate(p)[1],step)
        part=tensor.transpose(0,2,1,3).reshape(117,117);parts[str(step)]=part
        residual=baseline-parts['0.0001'];candidate=residual+part
        metrics.append(dict(step_Bohr=step,CDS_maximum_asymmetry_Ha_Bohr2=float(abs(part-part.T).max()),
            reassembled_maximum_asymmetry_Ha_Bohr2=float(abs(candidate-candidate.T).max()),
            maximum_change_from_default_Ha_Bohr2=float(abs(part-parts['0.0001']).max())))
    residual=baseline-parts['0.0001']
    energy_checks=[]
    for step in (1e-4,1e-5):
        numeric=energy_gradient(coordinates,lambda p:evaluate(p)[0],step)
        energy_checks.append(dict(step_Bohr=step,tested_components=117,
            maximum_energy_gradient_error_Ha_Bohr=float(abs(numeric-gradient).max())))
    evaluate(coordinates)
    arrays=output.with_suffix('.npz');np.savez_compressed(arrays,saved_full=baseline,residual_without_default_CDS=residual,
        **{'CDS_'+key:value for key,value in parts.items()})
    paths=[namespace_path,plan_path,stage,geometry,settings[0],matrices[0],Path(__file__),
           Path(gradient_hessian.__code__.co_filename),arrays]
    sources.update({str(p):sha256_file(p) for p in paths})
    result=dict(status='legacy_CDS_Hessian_step_component_study_not_full_QC',key=plan['source_id'],
        namespace=str(namespace_path),geometry_sha256=sha256_file(geometry),coordinates_source=str(settings[0]),
        solvent='acetonitrile',solvent_descriptors=descriptors,CDS_center_energy_Ha=energy,
        saved_full_maximum_asymmetry_Ha_Bohr2=float(abs(baseline-baseline.T).max()),
        residual_without_default_CDS_maximum_asymmetry_Ha_Bohr2=float(abs(residual-residual.T).max()),
        component_steps=metrics,energy_gradient_checks=energy_checks,
        fine_pair_maximum_change_Ha_Bohr2=float(abs(parts['1e-06']-parts['1e-07']).max()),
        numerical_components_only=True,CDS_energy_and_gradient_functions_changed=False,
        electrostatic_model_changed=False,full_DFT_Hessian_recomputed=False,CUDA_initialized=False,
        native_G_qualified=False,redox_potential_qualified=False,source_sha256=sources)
    verify(namespace_path)
    if 'cupy' in sys.modules:raise ValueError('Unexpected CUDA dependency in CPU component analysis')
    atomic_write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('spec','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=analyze(a.spec.resolve(),a.out.resolve())
    print(json.dumps({k:r[k] for k in ('key','saved_full_maximum_asymmetry_Ha_Bohr2',
        'residual_without_default_CDS_maximum_asymmetry_Ha_Bohr2','component_steps','energy_gradient_checks')}))
