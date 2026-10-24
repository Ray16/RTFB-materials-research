"""A finite maximum cannot conceal NaN in another energy component."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('point_review',ROOT/'scripts/validation/review_polymer_resource_point.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def data():return dict(result=dict(e_smd_eV=-10.,e_gas_eV=-9.,converged_smd=True,converged_gas=True,charge=0,mult=1,solvent='acetonitrile'),gradient_max_Ha_Bohr=.1)


@pytest.mark.parametrize('key',['e_smd_eV','e_gas_eV'])
def test_each_nonfinite_energy_is_rejected_before_maximum(key):
    row=data();row['result'][key]=float('nan')
    with pytest.raises(ValueError):MODULE.checked_numbers(row)


def test_nonconverged_finite_energy_and_invalid_gradient_are_rejected():
    row=data();row['result']['converged_gas']=False
    with pytest.raises(ValueError):MODULE.checked_numbers(row)


def test_geometry_description_is_not_mistaken_for_the_bound_input_path(tmp_path):
    xyz=tmp_path/'input.xyz';xyz.write_text('1\npoint\nHe 0 0 0\n')
    row=data();row['geometry']='starting geometry only; resource qualification'
    row['result'].update(xyz=str(xyz),energy_geom_sha1=MODULE.geom_sha1(xyz))
    assert MODULE.point_geometry(row)==xyz
    xyz.write_text('changed')
    with pytest.raises(ValueError):MODULE.point_geometry(row)
    row=data();row['gradient_max_Ha_Bohr']=-.1
    with pytest.raises(ValueError):MODULE.checked_numbers(row)
