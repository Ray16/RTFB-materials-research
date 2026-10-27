#!/usr/bin/env python3
"""Assemble five independently stationary CPU terms into gas counterpoise.

Every energy remains bound to its exact geometry, common grid, method,
checkpoint review and orbital-stationarity result. No solvent correction or
association free energy is inferred.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

ROLES=('AB','A_own','B_own','A_ghost','B_ghost')


def assemble(rows,core):
    if len(rows)!=5 or {r['role'] for r in rows}!=set(ROLES):raise ValueError('Exactly five unique counterpoise terms required')
    settings=rows[0]['settings'];grid=rows[0]['common_grid_hashes']
    if settings!=core.SETTINGS or set(grid)!={'grid_coordinates_AU.npy','grid_weights.npy'}:
        raise ValueError('Counterpoise method or common-grid inventory differs')
    for row in rows:
        if (row['settings']!=settings or row['common_grid_hashes']!=grid or not row['stationarity_passed']
                or not isinstance(row['E_Ha'],(float,int)) or isinstance(row['E_Ha'],bool) or not math.isfinite(row['E_Ha'])):
            raise ValueError('Term method, grid, energy or stationarity differs')
    return core.sums({row['role']:row['E_Ha'] for row in rows})


def collect(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,remaining_root,reviews_root,out):
    if out.exists():raise FileExistsError('Preserve previous collected interaction')
    worker=api.load_module('remaining_collection_dependency_gate',ROOT/'scripts/pipeline/counterpoise_cpu_remaining.py')
    plan,sources,_,_=worker.inputs(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity)
    records=[]
    for role,raw,gradient_path,expected_gradient in (
        ('A_own',own,own_stationarity,'CPU_pilot_energy_and_orbital_stationarity_checked'),
        ('A_ghost',ghost,ghost_stationarity,'CPU_ghost_energy_and_orbital_stationarity_checked')):
        state=json.loads((raw/role/'energy.json').read_text());gradient=json.loads(gradient_path.read_text())
        if gradient['status']!=expected_gradient or not gradient['passed']:raise ValueError('A-term stationarity not qualified')
        records.append(dict(role=role,E_Ha=state['E_Ha'],settings=state['settings'],common_grid_hashes=state['common_grid_hashes'],stationarity_passed=True))
    for role in ('AB','B_own','B_ghost'):
        raw=remaining_root/role;state_path=raw/role/'energy.json';launch_path=raw/'launch.json';status_path=raw/'status.json'
        state=json.loads(state_path.read_text());status=json.loads(status_path.read_text());launch=json.loads(launch_path.read_text())
        review_path=reviews_root/role/'review.json';gradient_path=reviews_root/role/'stationarity/summary.json'
        review=json.loads(review_path.read_text());gradient=json.loads(gradient_path.read_text())
        if (state!=status or state['status']!='CPU_remaining_term_converged_independent_review_pending' or state['role']!=role
                or state['source_sha256']!=launch['source_sha256'] or state['qualified']
                or review['status']!='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending' or review['role']!=role
                or gradient['status']!='CPU_remaining_term_energy_and_orbital_stationarity_checked' or gradient['role']!=role or not gradient['passed']
                or abs(review['energy_Ha']-state['E_Ha'])>1e-9 or abs(gradient['reference_energy_Ha']-state['E_Ha'])>1e-9
                or gradient['full_counterpoise_computed'] or gradient['solvent_counterpoise_computed']):
            raise ValueError('Remaining term lacks matching independent review or stationarity: '+role)
        api.verify_sources(state['source_sha256']);api.verify_sources(review['source_sha256']);api.verify_sources(gradient['source_sha256'])
        bound=[launch_path,status_path,state_path,raw/role/'SCF.h5',raw/'grid_coordinates_AU.npy',raw/'grid_weights.npy']
        if any(gradient['source_sha256'].get(str(p))!=api.digest(p) or review['source_sha256'].get(str(p))!=api.digest(p) for p in bound):
            raise ValueError('Review and gradient do not bind the exact term checkpoint: '+role)
        sources.update({str(p):api.digest(p) for p in [*bound,review_path,gradient_path]})
        records.append(dict(role=role,E_Ha=state['E_Ha'],settings=state['settings'],common_grid_hashes=state['common_grid_hashes'],stationarity_passed=True))
    core=api.load_module('remaining_cpu_counterpoise_sums',ROOT/'scripts/validation/local_counterpoise.py');sums=assemble(records,core)
    sources.update({str(p):api.digest(p) for p in (Path(__file__).resolve(),ROOT/'scripts/validation/local_counterpoise.py',ROOT/'scripts/pipeline/counterpoise_cpu_remaining.py')})
    api.verify_sources(sources)
    result=dict(status='five_CPU_gas_counterpoise_terms_stationary',checked=time.time(),source_construction_key=plan['source_construction_key'],
        terms=records,sums=sums,source_sha256=sources,full_counterpoise_computed=True,solvent_counterpoise_computed=False,
        thermal_free_energy_computed=False,primary_candidate_admission=False,
        scope='Frozen local capped-fragment gas electronic counterpoise only; not a MeCN correction, free energy, coupling or polymer ensemble')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('case','geometry-review','own','own-stationarity','ghost','ghost-stationarity','remaining-root','reviews-root','out'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();values=[getattr(a,k).resolve() for k in ('case','geometry_review','own','own_stationarity','ghost','ghost_stationarity','remaining_root','reviews_root','out')]
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=collect(*values);print(json.dumps(dict(status=result['status'],sums=result['sums'])))
