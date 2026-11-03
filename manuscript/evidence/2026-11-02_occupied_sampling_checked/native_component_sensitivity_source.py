"""Context-checked reconstruction of a complete primary-model G control cohort.

Callers independently review the actual SCFs, Hessians and partitions first.
This reader checks the declared state set, atom/charge balance and arithmetic;
it cannot turn an incomplete or alternate-model cohort into a primary result.
"""
from collections import Counter
from copy import deepcopy
import math

from redox.core.common import EV_KJ
from redox.properties.context import combine_properties

CUTOFFS=('25','50','100')


def reconstruct(package,components):
    states={r['id']:r for r in package['states']}
    if not states or not package['reactions'] or len(states)!=len(package['states']) or set(states)!=set(components):
        raise ValueError('Exactly one independently checked component for every package state required')
    if len({r['id'] for r in package['reactions']})!=len(package['reactions']):raise ValueError('Distinct declared cycle identities required')
    for key,row in components.items():
        if (row.get('qc')!='independently_checked' or not row.get('sources')
                or (row['charge'],row['mult'])!=(states[key]['charge'],states[key]['mult'])
                or row.get('native_component_requalified') is not True or not row.get('symbols')
                or row['context']['quantity']!='composite_solution_G' or row['context']['unit']!='eV'
                or row['context']['solvent']!=package['active_protocol']['solvent']
                or row['context']['temperature_K']!=package['native_protocol']['temperature_K']
                or row['context']['electronic_protocol']!=package['active_protocol']
                or row['context']['thermal_definition']['native_protocol']!=package['native_protocol']
                or row['context']['solvent_model']!='SMD'
                or set(row['G_by_cutoff_eV'])!=set(CUTOFFS)):
            raise ValueError('Exact reviewed primary-model state, method and complete G values required')
        for cutoff,value in row['G_by_cutoff_eV'].items():
            if type(value) not in (int,float) or not math.isfinite(value):raise ValueError('Finite native G required')
    rows=[]
    for reaction in package['reactions']:
        coefficients=reaction['coefficients'];elements=Counter();charge=0
        if (not coefficients or not set(coefficients)<=set(states)
                or any(type(c) is not int or c==0 for c in coefficients.values())):
            raise ValueError('Exact nonzero integer state coefficients required')
        for key,c in coefficients.items():
            for symbol in components[key]['symbols']:elements[symbol]+=c
            charge+=c*components[key]['charge']
        if charge or any(elements.values()) or sum(coefficients.values())!=0:
            raise ValueError('This reference-cancelled cohort requires atom, charge and component-count balance')
        if reaction['kind']=='one_electron_potential_vs_Fc':
            if reaction['electrons']!=1:raise ValueError('Exact one-electron Fc-reference convention required')
            factor,unit=-1.,'V'
        elif reaction['kind'] in ('grafting_free_energy','reaction_free_energy'):
            factor,unit=EV_KJ,'kJ/mol'
        else:raise ValueError('Unsupported primary cohort quantity')
        values,shifts,assemblies={},{},{}
        for cutoff in CUTOFFS:
            inputs={}
            for key in coefficients:
                row=components[key];context=deepcopy(row['context'])
                context['thermal_definition']['selected_rotor_cutoff_cm']=int(cutoff)
                inputs[key]=dict(value=row['G_by_cutoff_eV'][cutoff],context=context,
                    state=dict(identity=dict(id=key,smiles=states[key]['smiles']),charge=row['charge'],multiplicity=row['mult']),
                    qc='independently_checked',sources=row['sources'])
            assembled=combine_properties(inputs,coefficients)
            reference=math.fsum(c*states[k]['checked_G_by_cutoff_eV'][cutoff] for k,c in coefficients.items())*factor
            if abs(reference-reaction['expected_by_cutoff'][cutoff])>(1e-9 if unit=='V' else 1e-7):
                raise ValueError('Frozen reference cycle does not reproduce its declared value')
            values[cutoff]=assembled['value']*factor
            shifts[cutoff]=values[cutoff]-reference;assemblies[cutoff]=assembled
        rows.append(dict(id=reaction['id'],kind=reaction['kind'],unit=unit,
            corrected_by_cutoff=values,corrected_minus_source_by_cutoff=shifts,property_assemblies=assemblies))
    return dict(status='complete_primary_native_cohort_context_and_cycles_reconstructed',
        states=len(states),cycles=rows,corrected_native_components_recomputed=True,
        total_free_energy_cycles_numerically_requalified=True,
        conformer_ensemble_qualified=False,experimental_accuracy_validated=False,
        speciation_accuracy_validated=False,material_admission=False)
