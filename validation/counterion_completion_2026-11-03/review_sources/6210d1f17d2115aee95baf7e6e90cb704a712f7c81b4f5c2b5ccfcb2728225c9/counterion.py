"""Matched incoming-chloride association of an exact oxidized graft.

Inputs must already have independent physical QC. This module checks declared
chemical/thermodynamic compatibility and stoichiometry, not quantum convergence.
"""
from copy import deepcopy
import math

from rdkit import Chem

from redox.core.common import EV_KJ
from redox.properties.context import combine_properties
from redox.workflows.ionic_stationary_cycles import validate_cycle_balance

COEFFICIENTS = dict(paired=1, graft=-1, chloride=-1)


def association_cycle(components):
    if set(components) != set(COEFFICIENTS):
        raise ValueError('Exact paired, isolated graft and chloride components required')
    expected = dict(paired=(1, 1), graft=(2, 1), chloride=(-1, 1))
    cases = {}
    for key, record in components.items():
        state = record['state']
        if (state['charge'], state['multiplicity']) != expected[key] or record.get('thermodynamic_entities') != 1:
            raise ValueError('Exact singlet states and one statistical entity per component required')
        if set(record.get('G_by_cutoff_eV', {})) != {'25','50','100'}:
            raise ValueError('Every declared rotor setting is required')
        cases[key] = dict(smiles=state['smiles'], charge=state['charge'])
    canonical = lambda value: Chem.MolToSmiles(Chem.MolFromSmiles(value))
    if (canonical(cases['chloride']['smiles']) != '[Cl-]'
            or canonical(cases['paired']['smiles']) != canonical(cases['graft']['smiles']+'.[Cl-]')):
        raise ValueError('Association must retain the exact graft graph and add only disconnected chloride')
    validate_cycle_balance(dict(coefficients=COEFFICIENTS), cases)
    delta_n = sum(COEFFICIENTS[key]*record['thermodynamic_entities'] for key, record in components.items())
    source_by_cutoff, target_by_cutoff, checked = {}, {}, {}
    for cutoff in ('25','50','100'):
        records = {}
        for key, original in components.items():
            record = deepcopy(original)
            record['value'] = record['G_by_cutoff_eV'][cutoff]
            record['context']['thermal_definition']['evaluated_rotor_cutoff_cm'] = int(cutoff)
            records[key] = record
        combined = combine_properties(records, COEFFICIENTS)
        context = combined['context']
        if (context['quantity'] != 'composite_solution_G' or context['unit'] != 'eV'
                or context['standard_state'] != 'source xTB 1-atm ideal-translation convention; no 1-M conversion applied'):
            raise ValueError('The declared source translation convention is required before conversion')
        temperature = context['temperature_K']
        correction = delta_n*.00831446261815324*temperature*math.log(.08205736608096*temperature)
        source_by_cutoff[cutoff] = combined['value']*EV_KJ
        target_by_cutoff[cutoff] = source_by_cutoff[cutoff]+correction
        checked[cutoff] = combined
    return dict(status='finite_primary_microstate_association_cycle_computed', coefficients=dict(COEFFICIENTS),
        delta_n_solutes=delta_n, source_convention_deltaG_kJ_mol=source_by_cutoff,
        one_M_deltaG_kJ_mol=target_by_cutoff, one_atm_to_one_M_correction_kJ_mol=correction,
        standard_state_conversion='delta_n R T ln(R_atm T); source ideal translation to 1 mol/L',
        component_contexts=checked, pending_values_are_null=True,
        population_weights_inferred=False, conformer_ensemble_qualified=False,
        reaction_rate_or_yield_inferred=False, material_admission=False)
