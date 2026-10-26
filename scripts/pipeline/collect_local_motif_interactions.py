#!/usr/bin/env python
"""Join independently checked frozen states without borrowing missing energies."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

EV_KJ_MOL=1.602176634e-19*6.02214076e23/1000


def interaction(states):
    if set(states)!={'A','B','AB'}:raise ValueError('Three exact local energy parents required')
    if any(row['settings']!=states['A']['settings'] for row in states.values()):raise ValueError('Local electronic protocols differ')
    if any(row['task']['charge']!=0 or row['task']['mult']!=1 for row in states.values()):
        raise ValueError('This first local comparison is scoped to neutral closed-shell controls')
    result={}
    for medium in ('gas','acetonitrile'):
        values={name:row['energies'][medium] for name,row in states.items()}
        if not all(type(v) in (int,float) and math.isfinite(v) for v in values.values()):raise ValueError('Finite matched electronic energies required')
        result[medium]=math.fsum([values['AB'],-values['A'],-values['B']])*EV_KJ_MOL
    return result


def collect(construction,reviews,out):
    manifest_path=construction/'manifest.json';manifest=json.loads(manifest_path.read_text());api.verify_sources(manifest['source_sha256'])
    catalog={};proof={**manifest['source_sha256'],str(manifest_path):api.digest(manifest_path),str(Path(__file__).resolve()):api.digest(__file__)}
    for path in reviews:
        row=json.loads(path.read_text())
        if (row['status']!='frozen_local_electronic_state_and_resource_pilot_independently_checked'
                or row['primary_candidate_admission'] or row['thermochemistry_computed']):
            raise ValueError('Independent frozen-state review required')
        api.verify_sources(row['source_sha256']);plan=api.read_plan(Path(row['campaign']))
        key=plan['geometry_sha256']
        if key in catalog:raise ValueError('Ambiguous duplicate local geometry/state review')
        catalog[key]=row;proof[str(path)]=api.digest(path);proof.update(row['source_sha256'])
    rows=[]
    for case in manifest['entries']:
        found={};missing=[]
        for molecule in case['molecules']:
            path=Path(molecule['geometry'])
            if api.digest(path)!=molecule['geometry_sha256']:raise ValueError('Declared local geometry changed')
            proof[str(path)]=api.digest(path)
            row=catalog.get(molecule['geometry_sha256'])
            if row is None:missing.append(molecule['id']);continue
            if row['task']['atoms']!=molecule['atoms']:raise ValueError('Local state atom count differs')
            found[molecule['id']]=row
        values=None if missing else interaction(found)
        rows.append(dict(key=case['key'],parent_key=case['parent_key'],cap_length_A=case['cap_length_A'],
            parent_contact=case['contact'],missing_states=missing,
            interaction_energy_kJ_mol_by_medium=values,
            exact_state_energies_eV={name:row['energies'] for name,row in found.items()},
            free_energy_computed=False,standard_state_correction_applied=False,basis_superposition_correction_computed=False))
    result=dict(status='checked_frozen_local_interaction_sums_remaining_inputs_pending',checked=time.time(),
        declared_triplets=len(rows),completed_triplets=sum(r['interaction_energy_kJ_mol_by_medium'] is not None for r in rows),
        checked_electronic_geometries=len(catalog),states=rows,source_sha256=proof,conversion_kJ_mol_per_eV=EV_KJ_MOL,
        definition='E_AB-E_A-E_B with unchanged coordinates, electronic protocol, state and medium',
        primary_candidate_admission=False,free_energies_computed=False,electronic_coupling_computed=False,
        population_or_loading_inferred=False,full_polymer_nonadditivity_supplied=False,
        scope='Uncorrected frozen methyl-capped local motif interaction, including caps and continuum response. No pure-pi binding, thermal association free energy, full-backbone decomposition or synthesis eligibility inferred.')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--construction',type=Path,required=True);p.add_argument('--review',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=collect(a.construction.resolve(),[r.resolve() for r in a.review],a.out.resolve());print(json.dumps(dict(completed=row['completed_triplets'],
        values=[(r['key'],r['interaction_energy_kJ_mol_by_medium']) for r in row['states'] if r['interaction_energy_kJ_mol_by_medium'] is not None])))


if __name__=='__main__':main()
