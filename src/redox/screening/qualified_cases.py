"""Exact-graph attachment of qualified computational cases to candidate registries.

Replays declared arithmetic and checks case identity; it does not repeat raw
quantum QC or turn a finite molecular result into a material admission.
"""
from copy import deepcopy
import hashlib
import math
import json
from pathlib import Path

from rdkit import Chem

from redox.properties.context import CONTEXT_FIELDS
from redox.properties.reference_redox import neutral_oxidation_vs_fc
from redox.screening.feasibility import assess_candidate
from redox.core.hashing import sha256_file

SCHEMA='qualified_neutral_oxidation_case_v2'
MODEL_SCOPE='exact capped local graft; no full polymer or explicit counterion model'


def identity(smiles):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None:raise ValueError('Valid exact product graph required')
    canonical=Chem.MolToSmiles(mol,isomericSmiles=True)
    return canonical,Chem.MolToInchiKey(mol),hashlib.sha256(canonical.encode()).hexdigest()


def validate_case(case):
    if (case.get('schema')!=SCHEMA or case.get('status')!='qualified_exact_finite_redox_case'
            or case.get('material_admission') is not False
            or case.get('experimental_accuracy_validated') is not False
            or case.get('conformer_ensemble_qualified') is not False
            or case.get('model_scope')!=MODEL_SCOPE
            or case.get('all_axis_decision',{}).get('status')!='pending'
            or not case.get('source_bindings')):
        raise ValueError('Exact finite-case qualification and unresolved material scope required')
    smiles,key,graph=identity(case['product_smiles'])
    if (case['product_inchikey'],case['registry_graph_key'])!=(key,graph):
        raise ValueError('Product graph identities differ')
    cycles=case['redox_cycles']
    if set(cycles)!={'minus','plus'}:raise ValueError('Both declared neutral branches required')
    for cycle in cycles.values():
        if neutral_oxidation_vs_fc(cycle['components'])!=cycle:raise ValueError('Redox cycle differs from replay')
        for role in ('graft_neu','graft_ox'):
            if identity(cycle['components'][role]['state']['graph_identity'])[0]!=smiles:
                raise ValueError('A source-parent or different graft property cannot be transferred')
    prop=case['properties']['first_oxidation']
    context=prop['context']
    expected=cycles['minus']['potential_by_cutoff_V']['50']
    values=[v for c in cycles.values() for v in c['potential_by_cutoff_V'].values()]
    if (prop['value_QC']!='qualified_finite_model' or prop['unit']!='V'
            or type(prop['value']) not in (int,float) or not math.isfinite(prop['value'])
            or prop['value']!=expected or prop['scenario_bounds_V']!=[min(values),max(values)]
            or any(k not in context or context[k] is None for k in CONTEXT_FIELDS)
            or context['quantity']!='one_electron_oxidation_potential' or context['unit']!='V'
            or context['reference']['electrode']!='Fc/Fc+'):
        raise ValueError('Property value, context, reference or scenario bounds differ')
    source_context=cycles['minus']['components']['graft_neu']['context']
    for key in ('electronic_protocol','solvent','solvent_model','temperature_K','standard_state'):
        if context[key]!=source_context[key]:raise ValueError('Published physical context differs')
    thermal=deepcopy(source_context['thermal_definition']);thermal['selected_rotor_cutoff_cm']=50
    reference=dict(electrode='Fc/Fc+',definition=cycles['minus']['definition'],electrons=1,
                   delta_n_solutes=0,source_component_reference=source_context['reference'])
    if context['thermal_definition']!=thermal or context['reference']!=reference:
        raise ValueError('Thermal or reference definition differs')
    if (context.get('model_scope')!=MODEL_SCOPE or prop.get('component_states')!={
            branch:{role:row['state'] for role,row in cycle['components'].items()}
            for branch,cycle in cycles.items()}):
        raise ValueError('Exact model scope and component charge/spin identities required')
    if prop['component_numerical_policies']!={k:r['native_numerical_policy'] for k,r in cycles['minus']['components'].items()}:
        raise ValueError('Component numerical policies were lost or changed')
    if any(value is not None for key,value in case['properties'].items() if key!='first_oxidation'):
        raise ValueError('This case does not qualify uncomputed properties')
    if set(case['properties'])!={'first_oxidation','second_oxidation','grafting_free_energy',
            'solution_relaxation','occupancy_nonadditivity','motif_association',
            'electronic_coupling','accessible_capacity'}:
        raise ValueError('Unresolved property fields must remain explicit')
    return case


def assemble(redox, registry, policy, bindings, alias='phenothiazine'):
    if (redox['status']!='phenothiazine_finite_oxidation_potential_context_and_arithmetic_checked'
            or not redox['redox_potential_qualified'] or redox['experimental_accuracy_validated']
            or redox['conformer_ensemble_qualified'] or redox['material_admission']):
        raise ValueError('Qualified finite oxidation evidence required')
    matches=[r for r in registry['candidates'] if alias in r['aliases']]
    if len(matches)!=1 or matches[0]['route_status']!='proceed':
        raise ValueError('One exact currently route-admitted registry graph required')
    candidate=matches[0];cycles=deepcopy(redox['cycles']);primary=cycles['minus']
    context=deepcopy(primary['components']['graft_neu']['context'])
    context.update(quantity='one_electron_oxidation_potential',unit='V',
        model_scope=MODEL_SCOPE,
        reference=dict(electrode='Fc/Fc+',definition=primary['definition'],electrons=1,
                       delta_n_solutes=0,source_component_reference=deepcopy(context['reference'])))
    context['thermal_definition']['selected_rotor_cutoff_cm']=50
    values=[v for c in cycles.values() for v in c['potential_by_cutoff_V'].values()]
    prop=dict(value=primary['potential_by_cutoff_V']['50'],unit='V',value_QC='qualified_finite_model',
        context=context,scenario_bounds_V=[min(values),max(values)],
        bounds_scope='Retained branches and rotor settings; not an experimental confidence interval',
        primary_neutral_branch='minus',neutral_branches_population_weighted=False,
        component_states={branch:{role:deepcopy(row['state']) for role,row in cycle['components'].items()}
                          for branch,cycle in cycles.items()},
        component_numerical_policies={k:r['native_numerical_policy'] for k,r in primary['components'].items()})
    source=lambda key:[dict(path=bindings[key]['original_path'],sha256=bindings[key]['sha256'])]
    axes={a['id']:dict(product_key=candidate['inchikey'],scope=a['scope'],status='missing',
        reason='Exact application criterion or evidence unresolved') for a in policy['required_axes']}
    axes['attachment_route'].update(status='qualified',value='proceed',sources=source('registry'))
    axes['redox_window'].update(status='partial',model_value_V_vs_Fc=prop['value'],sources=source('redox'),
        reason='Finite molecular potential qualified; application cutoff, exact-graft accuracy and material transfer unresolved')
    axes['state_stability'].update(status='partial',sources=source('redox'),
        reason='Declared-state minima qualified; global spin/ensemble ordering and degradation remain unresolved')
    decision=assess_candidate(dict(id=alias,product_key=candidate['inchikey'],role='candidate',
                                   electrode_pools=['catholyte'],axes=axes),policy)
    result=dict(schema=SCHEMA,status='qualified_exact_finite_redox_case',candidate_alias=alias,
        product_smiles=candidate['product_smiles'],product_inchikey=candidate['inchikey'],
        registry_graph_key=candidate['graph_key'],registry_aliases=candidate['aliases'],
        model_scope=MODEL_SCOPE,
        properties=dict(first_oxidation=prop,second_oxidation=None,grafting_free_energy=None,
                        solution_relaxation=None,occupancy_nonadditivity=None,motif_association=None,
                        electronic_coupling=None,accessible_capacity=None),
        redox_cycles=cycles,field_QC=dict(first_oxidation='qualified finite molecular value',
            remaining_properties='unqualified or uncomputed; no family or parent transfer'),
        all_axis_decision=decision,material_admission=False,experimental_accuracy_validated=False,
        conformer_ensemble_qualified=False,source_bindings=deepcopy(bindings))
    return validate_case(result)


def attach(records, case, source_name, source_digest):
    validate_case(case)
    matches=[r for r in records if r['graph_key']==case['registry_graph_key']]
    if len(matches)!=1:raise ValueError('Exactly one matching registry graph required')
    target=matches[0]
    if (identity(target['product_smiles'])!=identity(case['product_smiles'])
            or target['inchikey']!=case['product_inchikey'] or case['candidate_alias'] not in target['aliases']):
        raise ValueError('Case and registry identity differ')
    if any(r.get('source')==source_name for r in target['evaluations']):raise ValueError('Case source already attached')
    target['evaluations'].append(dict(source=source_name,row=0,kind=SCHEMA,source_sha256=source_digest,
        model_scope=case['model_scope'],
        reported_values={'first_oxidation_V_vs_Fc':case['properties']['first_oxidation']['value']},
        properties=deepcopy(case['properties']),coverage=deepcopy(case['field_QC']),
        raw_records_revalidated=False,derived_arithmetic_replayed=True,
        material_admission=False,route_admission_inferred_from_property=False,
        all_axis_status=case['all_axis_decision']['status']))


def read_case(path):
    """Reconstruct from frozen local inputs without reusing private runtime paths."""
    path=Path(path);case=json.loads(path.read_text());inputs={}
    for key,binding in case['source_bindings'].items():
        relative=Path(binding['snapshot'])
        if relative.is_absolute() or '..' in relative.parts:raise ValueError('Unsafe case snapshot path')
        source=path.parent/relative
        if not source.resolve().is_relative_to(path.parent.resolve()) or sha256_file(source)!=binding['sha256']:
            raise ValueError('Candidate case snapshot changed')
        inputs[key]=json.loads(source.read_text())
    if set(inputs)!={'redox','registry','policy'}:raise ValueError('Complete frozen case inputs required')
    expected=assemble(inputs['redox'],inputs['registry'],inputs['policy'],case['source_bindings'],case['candidate_alias'])
    if {k:v for k,v in case.items() if k!='publisher_sha256'}!=expected:
        raise ValueError('Published case differs from its frozen reconstruction')
    return case
