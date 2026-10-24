"""A requested fast mode must not silently change state or become another backend."""
from types import SimpleNamespace

import pytest

from redox.qm.uma_nitro import state_key, check_actual_nitro


def test_molecular_charge_and_multiplicity_are_explicit():
    assert state_key([8,1,1],0,1)[1:3] == (0,1)
    assert state_key([8,1,1],1,2)[1:3] == (1,2)
    assert state_key([8,1,1],-1,2)[1:3] == (-1,2)
    assert state_key([8,1,1],0,1) != state_key([8,1,1],0,3)


@pytest.mark.parametrize('charge,mult', [(0,2),(1,1),(0,0),(True,1),(0.,1)])
def test_inconsistent_or_implicit_state_is_rejected(charge,mult):
    with pytest.raises(ValueError):
        state_key([8,1,1],charge,mult)


def predictor(nitro, module):
    executor = type('Executor', (), {'__module__':module})()
    return SimpleNamespace(inference_settings=SimpleNamespace(nitro=nitro),
        model=SimpleNamespace(module=executor))


def test_actual_standalone_executor_is_required():
    assert check_actual_nitro(predictor(True,'fairchem.core.models.uma.uma_nitro'))['actual_mode']=='nitro'
    with pytest.raises(RuntimeError,match='fell back'):
        check_actual_nitro(predictor(False,'fairchem.core.models.uma.uma_nitro'))
    with pytest.raises(RuntimeError,match='not the standalone'):
        check_actual_nitro(predictor(True,'fairchem.core.models.uma.escn_md'))
