from copy import deepcopy
import json
from pathlib import Path
import pytest
from redox.validation.context_consequences import assess

NAMES=('base_assisted.json','pmdi_first_site_grafting.json','pmdi_capacity.json','spacing_conformers.json','tempo_solution_raw_result.json')


def inputs():
    root=Path(__file__).resolve().parents[1]/'manuscript/evidence/2026-10-31_xtb_component_checked'
    return [json.loads((root/name).read_text()) for name in NAMES]


def test_stated_context_changes_reconstruct_without_material_promotion():
    r=assess(*inputs())
    assert r['route_definition']['changes_at_declared_nonpositive_screen']==4
    assert r['route_definition']['complete']==5 and r['route_definition']['unchanged_pairwise_orders']==10
    assert r['polymer_context']['positive_penalties']==6 and r['polymer_context']['site_preference_changes_sign']
    assert r['capacity']['nominal_to_one_graft_bound_ratio']>2.5
    assert r['directional_response']['symmetric_parameter_qualified'] is False
    assert r['experimental_material_decisions_validated'] is False


@pytest.mark.parametrize('bad',['base_offset','cancellation','capacity','nan'])
def test_changed_chemical_or_numerical_interpretation_is_rejected(bad):
    x=inputs()
    if bad=='base_offset':x[0]['candidates'][0]['base_models']['dbu']['free_ions_kJ_mol']+=1
    if bad=='cancellation':x[1]['singles'][0]['by_cutoff']['50']['first_graft_kJ_mol']['contact_pair']+=1
    if bad=='capacity':x[2]['exact_graphs_at_handle_density_0_6'][1]['two_e_upper_bound_mAh_g']+=1
    if bad=='nan':x[4]['summary']['lambda_reduction_eV']=float('nan')
    with pytest.raises(ValueError):assess(*x)
