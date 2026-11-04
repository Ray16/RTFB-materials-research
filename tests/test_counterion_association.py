"""Reaction/speciation tests: cluster solute count is independent of graph dots."""
from copy import deepcopy
import math
import pytest

from redox.properties.counterion import association_cycle


def components():
    graft = 'Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1'
    context = dict(quantity='composite_solution_G', unit='eV', electronic_protocol={'test': 'declared'},
        thermal_definition={'test': 'declared'}, solvent='acetonitrile', solvent_model='primary_SMD',
        temperature_K=298.15, standard_state='source xTB 1-atm ideal-translation convention; no 1-M conversion applied',
        reference={'energy':'absolute test components','spectator':'same omitted ion'})
    rows = {}
    for key, smiles, charge, value in [('paired',graft+'.[Cl-]',1,-30.),('graft',graft,2,-20.),('chloride','[Cl-]',-1,-10.)]:
        rows[key] = dict(state=dict(identity=key, smiles=smiles, charge=charge, multiplicity=1),
            G_by_cutoff_eV={c:value for c in ('25','50','100')}, context=deepcopy(context),
            qc='independently_checked', sources={'scope':'arithmetic test fixture; no physical admission'},
            thermodynamic_entities=1)
    return rows


def test_one_cluster_standard_state_correction_is_not_graph_fragment_count():
    result = association_cycle(components())
    assert result['delta_n_solutes'] == -1
    assert result['source_convention_deltaG_kJ_mol'] == {'25':0.,'50':0.,'100':0.}
    assert -7.93 < result['one_atm_to_one_M_correction_kJ_mol'] < -7.92
    assert result['one_M_deltaG_kJ_mol']['50'] == result['one_atm_to_one_M_correction_kJ_mol']
    assert not result['population_weights_inferred'] and not result['reaction_rate_or_yield_inferred']


@pytest.mark.parametrize('mutation', ['reduced_radical','lost_chloride','changed_core','mixed_model',
    'mixed_method','pending','missing_cutoff','nonfinite','cluster_count','missing_reference','already_converted'])
def test_incompatible_or_unqualified_components_cannot_be_combined(mutation):
    rows = components()
    if mutation == 'reduced_radical': rows['paired']['state']['multiplicity'] = 2
    if mutation == 'lost_chloride': rows['paired']['state']['smiles'] = rows['graft']['state']['smiles']
    if mutation == 'changed_core': rows['paired']['state']['smiles'] = rows['paired']['state']['smiles'].replace('Cc1','CCc1')
    if mutation == 'mixed_model': rows['paired']['context']['solvent_model'] = 'experimental_SWIG'
    if mutation == 'mixed_method': rows['paired']['context']['electronic_protocol'] = {'xc':'other'}
    if mutation == 'pending': rows['paired']['qc'] = 'pending'
    if mutation == 'missing_cutoff': rows['paired']['G_by_cutoff_eV'].pop('25')
    if mutation == 'nonfinite': rows['paired']['G_by_cutoff_eV']['50'] = math.nan
    if mutation == 'cluster_count': rows['paired']['thermodynamic_entities'] = 2
    if mutation == 'missing_reference': rows['paired']['context'].pop('reference')
    if mutation == 'already_converted':
        for row in rows.values(): row['context']['standard_state'] = '1 mol/L'
    with pytest.raises(ValueError): association_cycle(rows)


def test_input_context_and_values_are_not_changed_by_cycle_assembly():
    rows = components(); original = deepcopy(rows)
    association_cycle(rows)
    assert rows == original


def test_editing_a_result_does_not_change_the_next_reaction():
    first = association_cycle(components())
    first['coefficients']['paired'] = 2
    second = association_cycle(components())
    assert second['coefficients'] == {'paired':1,'graft':-1,'chloride':-1}
    assert second['delta_n_solutes'] == -1
