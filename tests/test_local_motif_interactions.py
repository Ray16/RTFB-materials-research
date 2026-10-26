from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('local_motif_interaction_test',ROOT/'scripts/pipeline/collect_local_motif_interactions.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def states():
    return {name:dict(settings={'basis':'same'},task={'charge':0,'mult':1},energies={'gas':value,'acetonitrile':value+.1})
            for name,value in [('A',-1000.),('B',-1001.),('AB',-2001.01)]}


def test_local_electronic_nonadditivity_keeps_medium_and_contains_no_standard_state_term():
    result=module.interaction(states())
    assert result['gas']==pytest.approx(-.01*module.EV_KJ_MOL)
    assert result['acetonitrile']==pytest.approx(-.11*module.EV_KJ_MOL)


@pytest.mark.parametrize('change',['missing','protocol','charge','spin','NaN'])
def test_incomplete_or_mixed_local_states_are_rejected(change):
    row=deepcopy(states())
    if change=='missing':row.pop('B')
    elif change=='protocol':row['B']['settings']['basis']='other'
    elif change=='charge':row['B']['task']['charge']=-1
    elif change=='spin':row['AB']['task']['mult']=3
    else:row['A']['energies']['gas']=float('nan')
    with pytest.raises(ValueError):module.interaction(row)
