#!/usr/bin/env python3
"""Independent projected-partition and sampled-FD replay for a route reagent."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def review(plan_path,output):
    if output.exists():raise FileExistsError('Preserve existing reagent reviews')
    from redox.core.common import ROOT,read_result
    from redox.core.protocol import load_record,record_path
    from redox.qm.native_thermal import load_native_thermal,cache_path
    from redox.validation.protonation import hydrogen_topology_status
    plan=json.loads(plan_path.read_text());case=plan_path.parent;produced=json.loads((case/'result.json').read_text())
    if plan['role']!='reaction_component' or produced['role']!='reaction_component':raise ValueError('Explicit reagent role required')
    if any(sha(k)!=v for k,v in plan['source_sha256'].items()):raise ValueError('Prepared source hashes changed')
    gid,state=plan['source_id'],plan['state'];directory=ROOT/'calcs/dft'/gid/state
    raw=read_result(gid,state,raw=True);selected=read_result(gid,state)
    if ((raw['charge'],raw['mult'])!=(plan['charge'],plan['mult'])
            or selected.get('geometry_spin_qc','ok')!='ok' or selected.get('thermal_qc')!='ok'
            or selected.get('spin_ensemble_status') not in ('not_selected',None)
            or hydrogen_topology_status(plan['smiles'],directory/'opt.xyz')!='ok'):
        raise ValueError('Reagent identity, geometry, spin or thermal gate failed')
    energy=load_record(directory,plan['charge'],plan['mult']);thermal=load_native_thermal(directory,plan['charge'],plan['mult'])
    if energy is None or thermal is None:raise ValueError('Checked native and active-energy records required')
    surface,_=cache_path(directory,plan['charge'],plan['mult']);record=json.loads(surface.read_text())
    audit=load(ROOT/'scripts/validation/audit/native_partition.py','independent_reagent_partition')
    partition=audit.read_partition(surface,directory/'opt.xyz')
    if partition is None:raise ValueError('Independent projected partition failed')
    matrix_path=surface.with_suffix('')/'hessian_cartesian.npy';matrix=np.load(matrix_path,allow_pickle=False)
    native_review=load(case.parent/'sources/native_review.py','independent_reagent_FD')
    checks,paths=native_review.finite_difference_review(case,directory/'opt.xyz',record,matrix)
    values={str(c):energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][str(c)]
        +partition['g_thermal_total_eV']-partition['nuclear_by_cutoff_eV']['50'] for c in (25,50,100)}
    if any(not np.isfinite(v) or abs(v-produced['G_by_cutoff_eV'][c])>1e-8 for c,v in values.items()):
        raise ValueError('Reagent G does not independently reconstruct')
    files={plan_path,Path(__file__).resolve(),case/'result.json',directory/'opt.xyz',directory/'result.json',
        directory/'thermal_policy.json',record_path(directory,plan['charge'],plan['mult']),surface,matrix_path,
        ROOT/'scripts/validation/audit/native_partition.py',*paths}
    sources={str(p):sha(p) for p in files};sources.update(partition['evidence_sha256'])
    result=dict(status='reaction_component_native_numerically_reconstructed_resource_review_pending',role='reaction_component',
        source_id=gid,state=state,charge=plan['charge'],mult=plan['mult'],inchikey=plan['inchikey'],
        G_by_cutoff_eV=values,active_E_smd_eV=energy['e_smd_eV'],sampled_derivative_checks=checks,
        minimum_frequency_cm=record['min_frequency_cm'] if 'min_frequency_cm' in record else None,
        final_geometry_sha256=sha(directory/'opt.xyz'),source_sha256=sources,
        global_conformer_convergence_proven=False,material_admission=False)
    output.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=review(a.plan.resolve(),a.out.resolve());print(json.dumps(dict(status=r['status'],source_id=r['source_id'])))
