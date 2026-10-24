from types import SimpleNamespace

import numpy as np
import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def setup():
    module = api.load_module('solute_workspace_storage_test', ROOT / 'scripts/pipeline/pcm_solute_workspace.py')
    matrices = {k: np.array([[1., -0.], [2., 3.]]) for k in ('S', 'D', 'K_LU')}
    solvent = SimpleNamespace(_intermediates=dict(matrices, q=np.ones(2)))
    xp = SimpleNamespace(float64=np.float64, asnumpy=lambda v: v.copy(), asarray=lambda v: v.copy())
    metrics = dict(suspensions=0, restorations=0, exact_matrix_round_trips=0, maximum_suspended_bytes=0)
    return module, solvent, xp, metrics, matrices


@pytest.mark.parametrize('access', ['getitem', 'get', 'normal'])
def test_round_trip_restores_after_success_or_forbidden_solvent_access(access):
    module, solvent, xp, metrics, original = setup()
    releases = []
    def exercise():
        with module.suspend_matrices(solvent, xp, metrics, lambda: releases.append(1)):
            assert set(solvent._intermediates) == {'q'}
            if access == 'getitem':
                solvent._intermediates['S']
            elif access == 'get':
                solvent._intermediates.get('S')
    if access == 'normal':
        exercise()
    else:
        with pytest.raises(RuntimeError, match='suspended solvent matrix'):
            exercise()
    for key in original:
        assert solvent._intermediates[key].tobytes() == original[key].tobytes()
    assert metrics['suspensions'] == metrics['restorations'] == 1
    assert metrics['exact_matrix_round_trips'] == 3 and len(releases) == 2


def test_round_trip_rejects_even_signed_zero_changes():
    module, solvent, xp, metrics, _ = setup()
    def altered(value):
        result = value.copy(); result[result == 0] = 0.
        return result
    xp.asarray = altered
    with pytest.raises(ValueError, match='round trip changed values'):
        with module.suspend_matrices(solvent, xp, metrics, lambda: None):
            pass


def test_unexpected_solvent_reset_fails_closed_and_preserves_original_matrices():
    module, solvent, xp, metrics, original = setup()
    with pytest.raises(RuntimeError, match='intermediates changed'):
        with module.suspend_matrices(solvent, xp, metrics, lambda: None):
            solvent._intermediates = {}
    for key in original:
        assert solvent._intermediates[key].tobytes() == original[key].tobytes()
