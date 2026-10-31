"""Reject unsupported grid-order claims from altered frozen numerical evidence."""
import importlib.util
import json

import pytest

from redox.core.common import ROOT


def checker():
    spec = importlib.util.spec_from_file_location('screening_manuscript_check', ROOT/'manuscript/check_screening_controls.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.audit


@pytest.fixture
def evidence(tmp_path):
    for arch in ('V100', 'H100'):
        results = []
        differences = {}
        for variant, delta in [('installed', 1e-3), ('corrected', 1e-10)]:
            results.append(dict(variant=variant, values=[dict(sorted=True, grid_points=100,
                numeric_density_integral=1.+delta, XC_energy_Ha=-1.-delta),
                dict(sorted=False, grid_points=100, numeric_density_integral=1., XC_energy_Ha=-1.)]))
            differences[variant] = dict(density_integral_difference=(1.+delta)-1.,
                                       XC_energy_difference_Ha=(-1.-delta)+1.)
        (tmp_path/f'screening_{arch}_integration.json').write_text(json.dumps(dict(rows=results,native_G_qualified=False)))
        (tmp_path/f'screening_{arch}_review.json').write_text(json.dumps(dict(
            native_G_qualified=False,SCF_gradient_Hessian_qualification_complete=False,
            status='screening_masks_and_fixed_density_integration_independently_checked',
            upstream_absolute_limit=2e-9,integration_differences=differences)))
    return tmp_path


def test_frozen_control_arithmetic_reconstructs_both_devices(evidence):
    result=checker()(evidence)
    assert set(result['variants'])=={'V100','H100'} and not result['native_G_qualified']


@pytest.mark.parametrize('change',['altered_value','weakened_limit','unsupported_scope','missing_variant'])
def test_changed_or_promoted_control_is_rejected(evidence,change):
    suffix='integration' if change in ('altered_value','missing_variant') else 'review'
    path=evidence/f'screening_H100_{suffix}.json'
    data=json.loads(path.read_text())
    if change=='altered_value':data['rows'][1]['values'][0]['XC_energy_Ha']+=1e-5
    if change=='missing_variant':data['rows'].pop()
    if change=='weakened_limit':data['upstream_absolute_limit']=.1
    if change=='unsupported_scope':data['native_G_qualified']=True
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker()(evidence)
