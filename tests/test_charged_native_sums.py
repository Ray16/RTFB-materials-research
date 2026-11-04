import copy
import importlib.util
from pathlib import Path

import pytest


@pytest.fixture
def record():
    p=Path(__file__).resolve().parents[1]/'scripts/validation/review_charged_polymer_native.py';s=importlib.util.spec_from_file_location('checked_native_component_sums',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
    nuclear=dict(nuclear_by_cutoff_eV={'25':1.,'50':2.,'100':3.});energy=-10.;parent=energy*27.211386245988
    component=dict(charge=2,mult=1,spin_partition_eV=0.,primary_model_admission=False,material_admission=False,
        E_geometry_level_Ha=energy,E_geometry_level_eV=parent,nuclear_by_cutoff_eV=copy.deepcopy(nuclear['nuclear_by_cutoff_eV']),G_by_cutoff_eV={k:parent+v for k,v in nuclear['nuclear_by_cutoff_eV'].items()})
    return m,component,nuclear,energy


def test_native_component_closes_at_each_cutoff(record):
    m,component,nuclear,energy=record;assert m.component_sums(component,nuclear,energy)==component['G_by_cutoff_eV']


def test_neutral_component_requires_explicit_matched_reference_role(record):
    m,component,nuclear,energy=record;component['charge']=0
    with pytest.raises(ValueError):m.component_sums(component,nuclear,energy)
    component['matched_empty_polymer_reference']=True
    assert m.component_sums(component,nuclear,energy)==component['G_by_cutoff_eV']


@pytest.mark.parametrize('change',['parent','nuclear','sum','spin','scope','nonfinite','nonfinite_sum','nonfinite_nuclear'])
def test_mixed_or_unsupported_component_cannot_be_admitted(record,change):
    m,component,nuclear,energy=record
    if change=='parent':component['E_geometry_level_eV']+=.01
    if change=='nuclear':component['nuclear_by_cutoff_eV']['25']+=.01
    if change=='sum':component['G_by_cutoff_eV']['100']+=.01
    if change=='spin':component['spin_partition_eV']=.1
    if change=='scope':component['primary_model_admission']=True
    if change=='nonfinite':energy=float('nan')
    if change=='nonfinite_sum':component['G_by_cutoff_eV']['50']=float('nan')
    if change=='nonfinite_nuclear':component['nuclear_by_cutoff_eV']['25']=float('nan')
    with pytest.raises(ValueError):m.component_sums(component,nuclear,energy)
