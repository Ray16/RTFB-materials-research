#!/usr/bin/env python3
"""Assemble the exact lead's fully reviewed primary G implementation control."""
import argparse
from copy import deepcopy
import json
import math
import os
from pathlib import Path

import numpy as np
from ase.io import read as read_atoms

from redox.core.common import ROOT,HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import load,verify
from redox.workflows.ionic_stationary_cycles import screening_signature
from redox.validation.scf_frame import checkpoint_spin
from redox.validation.native_component_sensitivity import reconstruct


def collect(package,reviews,atomic_energy_review,output):
    manifest=json.loads((package/'manifest.json').read_text());components={};hashes={str(package/'manifest.json'):sha256_file(package/'manifest.json')}
    equations=load(ROOT/'scripts/validation/audit/native_partition.py','primary_cohort_atomic_nuclear_equations')
    context_base=dict(quantity='composite_solution_G',unit='eV',electronic_protocol=manifest['active_protocol'],
        thermal_definition=dict(native_protocol=manifest['native_protocol'],rotor_cutoffs_cm=[25,50,100],
            spin_partition='unsplit ground spin multiplet'),solvent='acetonitrile',solvent_model='SMD',
        temperature_K=298.15,standard_state='source xTB 1-atm ideal-translation convention; no 1-M conversion applied',
        reference=dict(energy='absolute component within primary composite protocol',
            ions='separated species; omitted spectators cancel only within declared cycles'))
    for state in manifest['states']:
        key=state['id'];geometry=package/state['files']['geometry']
        if sha256_file(geometry)!=manifest['files_sha256'][state['files']['geometry']]:raise ValueError('Frozen cohort geometry changed')
        atoms=read_atoms(geometry);symbols=atoms.get_chemical_symbols()
        if len(atoms)==1:
            if state['mult']!=1:raise ValueError('This atomic reference is an explicit singlet')
            path=atomic_energy_review;review=json.loads(path.read_text())
            if review['key']!=key or not review['actual_checkpoint_frames_verified'] or any(sha256_file(p)!=h for p,h in review['source_sha256'].items()):
                raise ValueError('Exact checked atomic electronic parent required')
            ns=Path(review['namespace']);spec=verify(ns);case=Path(spec['case']);energy_path=case/'corrected/result.json'
            energy=json.loads(energy_path.read_text());partition_path=package/state['files']['partition'];reference=json.loads(partition_path.read_text())
            if (sha256_file(partition_path)!=manifest['files_sha256'][state['files']['partition']]
                    or reference['temperature_K']!=298.15 or reference['internal_modes']!=0
                    or '1 atm' not in reference['definition']):raise ValueError('Exact monatomic translation convention required')
            value=float(equations.nuclear_equations(np.array([]),atoms.positions,np.array([reference['mass_amu']]),1,50))
            if abs(value-reference['G_thermal_eV'])>1e-8:raise ValueError('Independent monatomic translation differs')
            spin_rows=[]
            for settings_path in (case/'corrected/SCF_evaluations').glob('*/*/settings.json'):
                settings=json.loads(settings_path.read_text())
                if (settings['charge'],settings['mult'])==(state['charge'],state['mult']):spin_rows.append(checkpoint_spin(settings_path.with_name('SCF.h5')))
            if not spin_rows:raise ValueError('Actual atomic SCF wavefunction backing required')
            implementation=case/'corrected/implementation/implementation.json'
            signature=screening_signature(json.loads(implementation.read_text()),Path(spec['sources']),implementation.with_name('screening_kernel.cu'))
            atomic=dict(status='atomic_translation_and_saved_wavefunction_independently_checked',key=key,
                charge=state['charge'],mult=state['mult'],nuclear_eV=value,corrected_electronic_eV=energy['e_smd_eV'],
                G_by_cutoff_eV={c:energy['e_smd_eV']+value for c in ('25','50','100')},
                actual_SCF_wavefunction_spins_checked=len(spin_rows),no_internal_modes=True,new_Hessian_computed=False)
            atomic_path=output/'atomic_reference.json';atomic_write_json(atomic_path,atomic)
            values=atomic['G_by_cutoff_eV'];bound={str(path):sha256_file(path),str(atomic_path):sha256_file(atomic_path),str(partition_path):sha256_file(partition_path),str(energy_path):sha256_file(energy_path)}
            kind='analytic_monatomic_translation'
        else:
            path=reviews/(key+'.json');review=json.loads(path.read_text())
            if (review['status']!='primary_molecular_native_control_independently_reconstructed'
                    or review['key']!=key or (review['charge'],review['mult'])!=(state['charge'],state['mult'])
                    or not review['native_component_requalified'] or review['material_admission']
                    or review['actual_SCF_checkpoint_frames_checked']!=9 or review['actual_SCF_wavefunction_spins_checked']!=9
                    or any(sha256_file(p)!=h for p,h in review['source_sha256'].items())):
                raise ValueError('Exact fully reviewed primary native component required')
            spec=verify(Path(review['namespace']));case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
            if sha256_file(case/'state/opt.xyz')!=sha256_file(geometry):raise ValueError('Native cohort source geometry differs')
            source_native=Path(spec['archive'])/'src/redox/qm/native_thermal.py'
            native_module=load(source_native,'primary_cohort_bound_native_protocol_'+key)
            if native_module.NATIVE_THERMAL!=manifest['native_protocol']:raise ValueError('Actual native protocol differs from source cohort')
            values=review['G_by_cutoff_eV'];signature=review['numerical_implementation'];bound={str(path):sha256_file(path)};kind='new_corrected_primary_native_Hessian'
        context=deepcopy(context_base);context['reference']['numerical_implementation']=signature
        components[key]=dict(charge=state['charge'],mult=state['mult'],symbols=symbols,
            G_by_cutoff_eV=values,context=context,sources=bound,qc='independently_checked',
            native_component_requalified=True,component_kind=kind)
        hashes.update(bound)
    inputs=dict(package=manifest,components=components)
    result=reconstruct(manifest,components)
    atomic_write_json(output/'scalar_inputs.json',inputs);atomic_write_json(output/'cycles.json',result)
    hashes[str(output/'scalar_inputs.json')]=sha256_file(output/'scalar_inputs.json')
    hashes[str(output/'cycles.json')]=sha256_file(output/'cycles.json')
    summary=dict(status=result['status'],states=len(components),
        new_molecular_Hessians=sum(r['component_kind']=='new_corrected_primary_native_Hessian' for r in components.values()),
        analytic_monatomic_components=sum(r['component_kind']=='analytic_monatomic_translation' for r in components.values()),
        cycles=len(result['cycles']),source_sha256=hashes,
        maximum_potential_change_V=max(abs(v) for r in result['cycles'] if r['unit']=='V' for v in r['corrected_minus_source_by_cutoff'].values()),
        maximum_reaction_change_kJ_mol=max(abs(v) for r in result['cycles'] if r['unit']=='kJ/mol' for v in r['corrected_minus_source_by_cutoff'].values()),
        occupied_polymer_G_cycles_qualified=False,conformer_ensemble_qualified=False,experimental_accuracy_validated=False,
        producer_sha256=sha256_file(__file__))
    atomic_write_json(output/'summary.json',summary);return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('package','reviews','atomic-energy-review','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier full-cohort replays')
    a.out.mkdir(parents=True);r=collect(a.package.resolve(),a.reviews.resolve(),a.atomic_energy_review.resolve(),a.out.resolve())
    print(json.dumps({k:r[k] for k in ('status','states','cycles','maximum_potential_change_V','maximum_reaction_change_kJ_mol')}))
