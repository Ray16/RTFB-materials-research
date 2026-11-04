#!/usr/bin/env python3
"""Bind a failed optimizer's converged return to a separate fresh final solve."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

from ase.io import read

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);plan=json.loads((case/'plan.json').read_text())
    if not plan.get('optimizer_returned_input'):raise ValueError('Explicit optimizer-return finalization source required')
    fixed=load(ROOT/'scripts/validation/review_charged_polymer_pilot.py','fixed_ionic_point_review').review(namespace)
    geometry=load(ROOT/'scripts/validation/review_charged_polymer_optimization.py','final_ionic_geometry_checks')
    source=plan['source_optimization'];old=Path(source['source_case']);source_spec=verify(Path(source['namespace']))
    returned=json.loads((sources/'source_optimizer_return.json').read_text());point=json.loads((case/'evaluation.json').read_text())
    original_plan=json.loads((old/'plan.json').read_text())
    if (not returned['converged'] or returned['geometry_evaluations']!=len(list((old/'geometry_evaluations').glob('*.json')))
            or original_plan['prepared_case']!=plan['prepared_case']):raise ValueError('Exact converged source geometry history required')
    for path,digest in source['source_sha256'].items():
        if sha256_file(path)!=digest:raise ValueError('Optimizer-return source changed')
    final=dict(point,optimizer_converged=True)
    # The fixed-point review already binds the actual SCF checkpoint energy and
    # input geometry. Reuse those checked values for the joint state binding.
    positions=read(sources/'input.xyz').positions
    physical=geometry.physical_checks(original_plan,returned,final,positions,fixed['checkpoint_coordinates_A'],fixed['E_Ha'])
    retained=geometry.graph_checks(plan['prepared_case']['smiles'],Path(source_spec['sources'])/'input.xyz',sources/'input.xyz')
    return dict(status='charged_polymer_solution_stationary_geometry_independently_checked',namespace=str(namespace),
        candidate_alias=plan['exact_candidate_alias'],polymer_smiles=plan['prepared_case']['smiles'],
        graph_key=plan['prepared_case']['key'],charge=plan['charge'],mult=plan['mult'],atoms=plan['atoms'],E_Ha=fixed['E_Ha'],
        geometry_sha256=sha256_file(sources/'input.xyz'),quantum_settings=plan['quantum_settings'],
        solvent_model=fixed['solvent_model'],calculation_ion_convention=fixed['calculation_ion_convention'],
        geometry_evaluations=returned['geometry_evaluations'],**physical,**retained,peak_GPU_MiB=fixed['peak_GPU_MiB'],
        source_optimizer_namespace=source['namespace'],Hessian_minimum_check_complete=False,conformer_sampling_qualified=False,
        interaction_strength_computed=False,native_G_qualified=False,material_admission=False,
        source_sha256={**fixed['source_sha256'],**source['source_sha256']},producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve existing final-state reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','graph_key','E_Ha','maximum_atom_gradient_Ha_Bohr')}))
