"""A split CPU energy remains reusable without a second quantum calculation."""
from copy import deepcopy
import json

import pytest

from redox.core.protocol import ACTIVE_SP, geom_sha1, load_record
from redox.workflows.split_polymer import composite_result


def parents():
    geometry = dict(charge=0, mult=1, backend='gpu', optimized=True,
        geometry_optimization_converged=True, optimization_SCF_converged=True,
        geometry_only_result=True, opt_basis='def2-svp', opt_xc='r2scan')
    active = dict(protocol=deepcopy(ACTIVE_SP), geom_sha1='coordinates', charge=0, mult=1,
        source='computed:cpu', converged_smd=True, converged_gas=True,
        e_smd_eV=-10., e_gas_eV=-9.9, density_fit=True, conv_tol=1e-9)
    return geometry, active


def test_join_keeps_parents_unchanged_and_names_both_backends():
    geometry, active = parents()
    before = deepcopy((geometry, active))
    result = composite_result(geometry, active, 'coordinates', '/test/energy.json')
    assert (geometry, active) == before
    assert result['geometry_backend'] == 'gpu' and result['electronic_backend'] == 'cpu'
    assert result['optimized'] and not result['geometry_only_result']


@pytest.mark.parametrize('failure', ['geometry', 'protocol', 'state', 'unconverged', 'nonfinite'])
def test_mismatched_or_incomplete_parent_cannot_be_reused(failure):
    geometry, active = parents()
    if failure == 'geometry': active['geom_sha1'] = 'other coordinates'
    elif failure == 'protocol': active['protocol']['basis'] = 'def2-svp'
    elif failure == 'state': active['charge'] = -1
    elif failure == 'unconverged': geometry['optimized'] = False
    else: active['e_smd_eV'] = float('nan')
    with pytest.raises(ValueError): composite_result(geometry, active, 'coordinates', '/test/energy.json')


def test_native_sp_handoff_adopts_cpu_parent_without_recomputation(tmp_path, monkeypatch):
    from redox.qm import sp, dft
    (tmp_path/'opt.xyz').write_text('3\nsynthetic data-contract control\nO 0 0 0\nH .96 0 0\nH -.24 .93 0\n')
    geometry, active = parents()
    identity = geom_sha1(tmp_path/'opt.xyz')
    active['geom_sha1'] = identity
    result = composite_result(geometry, active, identity, '/test/original_cpu_record.json')
    (tmp_path/'result.json').write_text(json.dumps(result))
    def forbidden(*args, **kwargs):
        raise AssertionError('The checked CPU energy must not be recomputed')
    monkeypatch.setattr(dft, 'dft_smd', forbidden)
    status, record = sp.ensure('contract-control', 'rest', tmp_path, 0, 1, backend='gpu')
    assert status == 'adopted' and record['e_smd_eV'] == active['e_smd_eV']
    assert load_record(tmp_path, 0, 1)['e_gas_eV'] == active['e_gas_eV']
