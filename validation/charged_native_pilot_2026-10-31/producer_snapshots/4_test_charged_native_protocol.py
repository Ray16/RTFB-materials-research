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


@pytest.mark.parametrize('plan',[dict(charge=2),dict(charged_polymer_native_component=True,charge=-1,mult=1),
    dict(charged_polymer_native_component=True,charge=2,mult=3),dict(charge=True),
    dict(charged_polymer_native_component=True,charge=2,mult=1,saved_full_Hessian={'matrix':'other'})])
def test_unsupported_or_ambiguous_states_are_rejected(worker,plan):
    with pytest.raises(ValueError):worker.charge_state(plan)


def test_component_sampling_adds_actual_graft_atoms_without_dropping_protocol_columns(worker):
    native=SimpleNamespace(fd_columns=lambda gradient,symbols:[(0,1),(2,2)])
    gradient=np.zeros((4,3));symbols=['C','N','C','H']
    assert worker.declared_columns(native,gradient,symbols,dict(extra_FD_columns=[[1,2],[2,2]]))==[(0,1),(1,2),(2,2)]
    assert worker.declared_columns(native,gradient,symbols,dict(failed_column=dict(atom=3,axis=0)))==[(0,1),(2,2),(3,0)]
    with pytest.raises(ValueError):worker.declared_columns(native,gradient,symbols,dict(extra_FD_columns=[[5,2]]))
