#!/usr/bin/env python3
"""Reconstruct a private molecular native state and its saved FD columns."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def finite_difference_review(case,geometry,record,matrix):
    """Recompute the sampled errors from saved displaced SCF/gradient returns."""
    from ase.io import read
    coordinates=read(geometry).positions/.52917721067
    entries=[]
    for path in sorted((case/'SCF_evaluations').glob('*/*/settings.json')):
        settings=json.loads(path.read_text());gradient=path.with_name('gradient.npy');solved=path.with_name('SCF_result.json')
        if gradient.exists() and solved.exists() and json.loads(solved.read_text())['converged']:
            if (settings['charge'],settings['mult'])!=(record['identity']['charge'],record['identity']['mult']):
                raise ValueError('Checkpoint charge/spin differs from the native state')
            entries.append((path,settings,np.load(gradient,allow_pickle=False)))
    def gradient_at(target):
        matches=[entry for entry in entries if np.max(np.abs(np.asarray(entry[1]['coordinates_Bohr'])-target))<3e-9]
        if not matches:raise ValueError('Missing displaced gradient checkpoint')
        if any(not np.isfinite(r[2]).all() or np.max(np.abs(r[2]-matches[0][2]))>1e-8 for r in matches):
            raise ValueError('Repeated checkpoint gradients disagree')
        return matches[0][2],matches
    base,base_rows=gradient_at(coordinates)
    if abs(float(np.abs(base).max())-record['gradient_max'])>1e-10:
        raise ValueError('Native gradient maximum does not reconstruct')
    checks=[];used=[r[0] for r in base_rows]
    for row in record['derivative_checks']:
        atom,axis,step=row['atom_index'],row['axis'],row['step_Bohr']
        positive=coordinates.copy();positive[atom,axis]+=step
        negative=coordinates.copy();negative[atom,axis]-=step
        plus,rows_plus=gradient_at(positive);minus,rows_minus=gradient_at(negative)
        numerical=((plus-minus)/(2*step)).reshape(-1)
        error=float(np.linalg.norm(numerical-matrix[3*atom+axis])/np.linalg.norm(numerical))
        if not math.isfinite(error) or error>.01 or abs(error-row['relative_error'])>1e-10:
            raise ValueError('Sampled finite-gradient Hessian column does not reconstruct')
        checks.append(dict(atom_index=atom,axis=axis,step_Bohr=step,relative_error=error))
        used.extend(r[0] for r in rows_plus+rows_minus)
    expected={(0,1,h) for h in (.0005,.001)}
    symbols=read(geometry).get_chemical_symbols()
    other=next((i for i,s in enumerate(symbols) if s=='O'),next((i for i,s in enumerate(symbols) if s=='Fe'),len(symbols)-1))
    expected|={(other,2,h) for h in (.0005,.001)}
    if {(c['atom_index'],c['axis'],c['step_Bohr']) for c in checks}!=expected or len(checks)!=4:
        raise ValueError('Two complete sampled columns at two steps are required')
    paths=set()
    for path in used:
        paths.update([path,path.with_name('gradient.npy'),path.with_name('SCF_result.json')])
    return checks,paths


def review(plan_path,output):
    if output.exists():raise FileExistsError('Preserve existing numerical reviews')
    from redox.core.common import ROOT,read_result
    from redox.core.protocol import load_record,record_path
    from redox.qm.native_thermal import load_native_thermal,cache_path
    from redox.validation.protonation import hydrogen_topology_status
    from ase.io import read
    plan=json.loads(plan_path.read_text());case=plan_path.parent
    gid,state=plan['source_id'],plan['state'];directory=ROOT/'calcs/dft'/gid/state
    raw=read_result(gid,state,raw=True);selected=read_result(gid,state)
    if ((raw['charge'],raw['mult'])!=(plan['charge'],plan['mult'])
            or selected.get('geometry_spin_qc','ok')!='ok' or selected.get('thermal_qc')!='ok'
            or selected.get('spin_ensemble_status') not in ('not_selected',None)
            or hydrogen_topology_status(plan['product_smiles'],directory/'opt.xyz')!='ok'):
        raise ValueError('Exact native geometry, thermal selection or spin gate failed')
    energy=load_record(directory,plan['charge'],plan['mult']);thermal=load_native_thermal(directory,plan['charge'],plan['mult'])
    if energy is None or thermal is None:raise ValueError('Validated native surface and active energy required')
    surface,_=cache_path(directory,plan['charge'],plan['mult']);record=json.loads(surface.read_text())
    definition=importlib.util.spec_from_file_location('molecular_independent_nuclear_partition',ROOT/'scripts/validation/audit/native_partition.py')
    audit=importlib.util.module_from_spec(definition);definition.loader.exec_module(audit)
    partition=audit.read_partition(surface,directory/'opt.xyz')
    if partition is None:raise ValueError('Independent projected nuclear partition failed')
    matrix_path=surface.with_suffix('')/'hessian_cartesian.npy';matrix=np.load(matrix_path,allow_pickle=False)
    checks,checkpoint_paths=finite_difference_review(case,directory/'opt.xyz',record,matrix)
    values={str(c):energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][str(c)]
        +partition['g_thermal_total_eV']-partition['nuclear_by_cutoff_eV']['50'] for c in (25,50,100)}
    produced=json.loads((case/'result.json').read_text())
    if any(not math.isfinite(v) or abs(v-produced['G_by_cutoff_eV'][c])>1e-8 for c,v in values.items()):
        raise ValueError('Native free-energy sums do not reconstruct')
    files={plan_path,Path(__file__).resolve(),case/'result.json',directory/'opt.xyz',directory/'result.json',
        directory/'thermal_policy.json',record_path(directory,plan['charge'],plan['mult']),surface,matrix_path,
        ROOT/'scripts/validation/audit/native_partition.py',*checkpoint_paths}
    sources={str(p):sha(p) for p in files};sources.update(partition['evidence_sha256'])
    if any(sha(p)!=h for p,h in sources.items()):raise ValueError('Numerical review source changed')
    result=dict(status='molecular_native_numerically_reconstructed_resource_review_pending',
        candidate_alias=plan['candidate_alias'],source_id=gid,state=state,charge=plan['charge'],mult=plan['mult'],
        product_inchikey=plan['product_inchikey'],G_by_cutoff_eV=values,
        active_E_smd_eV=energy['e_smd_eV'],sampled_derivative_checks=checks,
        final_geometry_sha256=sha(directory/'opt.xyz'),source_sha256=sources,
        repaired_geometry_used=sha(directory/'opt.xyz')!=plan['source_geometry_sha256'],
        global_conformer_convergence_proven=False,material_admission=False,
        experimental_accuracy_inferred=False)
    output.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();result=review(a.plan.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],state=result['state'])))
