import copy
import importlib.util
from pathlib import Path

import pytest

spec=importlib.util.spec_from_file_location('numerical_repeat_control',Path(__file__).parents[1]/'scripts/pipeline/prepare_native_numerical_repeat.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def evidence():
    return dict(status='three_SCF_gradient_diagnostic_reconstructed_no_native_admission',full_native_G_qualified=False,
        original_failed_relative_error=.012,repeated_relative_error_against_source_Hessian=1e-5,
        points=[dict(energy_difference_Ha=1e-12,maximum_gradient_difference_Ha_Bohr=2e-6) for _ in range(3)])


def test_same_energy_nonreproducing_gradient_can_trigger_one_full_native_repeat(evidence):
    module.validate_control(evidence)


@pytest.mark.parametrize('change',['nan','incomplete','changed_energy','same_gradient','still_failed','already_native'])
def test_partial_changed_or_unsupported_controls_cannot_authorize_repeat(evidence,change):
    data=copy.deepcopy(evidence)
    if change=='nan':data['repeated_relative_error_against_source_Hessian']=float('nan')
    if change=='incomplete':data['points'].pop()
    if change=='changed_energy':data['points'][0]['energy_difference_Ha']=1e-5
    if change=='same_gradient':
        for row in data['points']:row['maximum_gradient_difference_Ha_Bohr']=1e-10
    if change=='still_failed':data['repeated_relative_error_against_source_Hessian']=.02
    if change=='already_native':data['full_native_G_qualified']=True
    with pytest.raises(ValueError):module.validate_control(data)
