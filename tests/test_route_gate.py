"""Protect the prospective chemistry gate from descriptor-only admissions."""
import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('route_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)


def reviewed():
    return dict(exact_precursor_access_reviewed=True,attachment_selectivity_reviewed=True,
                conditions_defined=True,balanced_route_reviewed=True,redox_core_retention_reviewed=True,
                source_urls=['https://example.org/exact-route'],
                conditions=dict(solvent='acetonitrile',reactive_form='defined precursor',
                                base='defined base or explicitly no added base',counterions='defined'))


def test_good_computed_descriptors_do_not_admit_unknown_routes():
    row=dict(id='candidate',frag='[*:1]Oc1ccccc1',deltaG_reference_kJ_mol=-100,capacity=300)
    result=gate.screen(row,{}, {})
    assert result['route_status']=='hold'
    assert not result['proceed_to_full_redox']
    assert result['proceed_to_feasibility_model']


def test_reviewed_route_requires_sources_and_conditions():
    row=dict(id='candidate',frag='[*:1]Oc1ccccc1')
    record=reviewed()
    assert gate.screen(row,{'candidate':record},{})['proceed_to_full_redox']
    del record['source_urls']
    assert not gate.screen(row,{'candidate':record},{})['proceed_to_full_redox']


def test_synthesis_hold_and_actual_unfavorable_route_override_admission():
    row=dict(id='candidate',frag='[*:1]Oc1ccccc1');record=reviewed()
    held=gate.screen(row,{'candidate':record},{'candidates':{'candidate':{'status':'on_hold'}}})
    assert not held['proceed_to_full_redox']
    record['strongly_unfavorable_under_actual_conditions']=True
    excluded=gate.screen(row,{'candidate':record},{})
    assert excluded['route_status']=='exclude_route'
    assert not excluded['proceed_to_full_redox']


def test_unreviewed_attachment_class_is_not_submitted():
    result=gate.screen(dict(id='candidate',frag='[*:1]CC'),{}, {})
    assert result['route_status']=='exclude_unreviewed_attachment'
    assert not result['proceed_to_feasibility_model']


def test_base_required_retains_candidate_without_skipping_route_evidence():
    row=dict(id='candidate',frag='[*:1]Oc1ccccc1',deltaG_reference_kJ_mol=58)
    assessment={'candidates':{'candidate':{'status':'base_required'}}}
    result=gate.screen(row,{},assessment)
    assert result['proceed_to_feasibility_model']
    assert not result['proceed_to_full_redox']
    assert any('Base required; retained' in reason for reason in result['reasons'])
    # Positive HCl-reference energy cannot veto a separately reviewed base route.
    assert gate.screen(row,{'candidate':reviewed()},assessment)['proceed_to_full_redox']
