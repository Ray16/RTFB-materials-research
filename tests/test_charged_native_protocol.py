"""Charged native extension preserves the control defaults and sampling gates."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


@pytest.fixture
def worker():
    p=Path(__file__).resolve().parents[1]/'scripts/validation/run_screening_polymer_hessian.py'
    s=importlib.util.spec_from_file_location('checked_charged_native_protocol',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def test_neutral_and_positive_closed_shell_states_are_explicit(worker):
    assert worker.charge_state({})==(0,1)
    assert worker.charge_state(dict(charged_polymer_native_component=True,charge=2,mult=1))==(2,1)
    assert worker.charge_state(dict(charged_polymer_native_component=True,charge=0,mult=1,matched_empty_polymer_reference=True,prepared_case=dict(occupied_sites=[])))==(0,1)


@pytest.mark.parametrize('plan',[dict(charge=2),dict(charged_polymer_native_component=True,charge=-1,mult=1),
    dict(charged_polymer_native_component=True,charge=2,mult=3),dict(charge=True),
    dict(charged_polymer_native_component=True,charge=0,mult=1),
    dict(charged_polymer_native_component=True,charge=0,mult=1,matched_empty_polymer_reference=True,prepared_case=dict(occupied_sites=[2])),
    dict(charged_polymer_native_component=True,charge=2,mult=1,saved_full_Hessian={'matrix':'other'})])
def test_unsupported_or_ambiguous_states_are_rejected(worker,plan):
    with pytest.raises(ValueError):worker.charge_state(plan)


def test_only_explicit_checked_full_matrix_can_enter_charged_validation(worker):
    plan=dict(charged_polymer_native_component=True,charge=4,mult=1,
        charged_saved_Hessian_validation=True,Hessian_validation_continuation=True,
        saved_full_Hessian=dict(charge=4,mult=1,source_response_independently_checked=True))
    assert worker.charge_state(plan)==(4,1)
    for field in ('charged_saved_Hessian_validation','Hessian_validation_continuation'):
        with pytest.raises(ValueError):worker.charge_state(dict(plan,**{field:False}))
    for field in ('resume_components','resume_response','resume_partial_Hessian'):
        with pytest.raises(ValueError):worker.charge_state(dict(plan,**{field:{'source':'other'}}))


def test_component_sampling_adds_actual_graft_atoms_without_dropping_protocol_columns(worker):
    native=SimpleNamespace(fd_columns=lambda gradient,symbols:[(0,1),(2,2)])
    gradient=np.zeros((4,3));symbols=['C','N','C','H']
    assert worker.declared_columns(native,gradient,symbols,dict(extra_FD_columns=[[1,2],[2,2]]))==[(0,1),(1,2),(2,2)]
    assert worker.declared_columns(native,gradient,symbols,dict(failed_column=dict(atom=3,axis=0)))==[(0,1),(2,2),(3,0)]
    with pytest.raises(ValueError):worker.declared_columns(native,gradient,symbols,dict(extra_FD_columns=[[5,2]]))


def test_only_bounded_small_empty_inputs_get_a_smaller_output_envelope():
    from redox.workflows.charged_polymer_pilot import native_output_profile
    small=native_output_profile(91,881,177,True)
    assert small['maximum_durable_output_bytes']==2*2**30
    assert small['estimated_retained_bytes']*small['retention_margin_factor']<small['maximum_durable_output_bytes']
    assert native_output_profile(114,1100,213,False)['maximum_durable_output_bytes']==4*2**30
    assert native_output_profile(91,881,177,False)['maximum_durable_output_bytes']==4*2**30
    assert native_output_profile(119,1193,235,True)['maximum_durable_output_bytes']==4*2**30
