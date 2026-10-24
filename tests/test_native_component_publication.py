from copy import deepcopy

import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


@pytest.fixture
def component():
    module = api.load_module('native_publication_test', ROOT / 'scripts/pipeline/publish_candidate_native_component.py')
    task = dict(key='pair__pose0', model_key='pair', occupied_sites=[2, 3])
    row = dict(key=task['key'], occupied_sites=[2, 3], state_dir='/raw/native', E_eV=-100., G_eV=-99.)
    review = dict(status='available_exact_polymer_native_components_independently_reconstructed',
        checked_components=1, states=[row], normal_admission_rechecked=True, QM_recomputed=False)
    energy = dict(task=deepcopy(task), E_smd_eV=-100.)
    thermal = dict(task=deepcopy(task), state_dir='/raw/native', G_eV=-99.)
    return module, review, energy, thermal


def test_checked_binding_returns_graph_key_not_pose_key(component):
    module, review, energy, thermal = component
    assert module.component_binding(review, energy, thermal) == 'pair'


@pytest.mark.parametrize('change', ['unchecked', 'wrong_site', 'wrong_state', 'nan_G', 'wrong_energy'])
def test_wrong_or_unchecked_native_component_cannot_reach_registry(component, change):
    module, review, energy, thermal = component
    if change == 'unchecked':
        review['normal_admission_rechecked'] = False
    elif change == 'wrong_site':
        review['states'][0]['occupied_sites'] = [2, 4]
    elif change == 'wrong_state':
        energy['task']['key'] = 'other_pose'
    elif change == 'nan_G':
        thermal['G_eV'] = float('nan')
    else:
        energy['E_smd_eV'] += .1
    with pytest.raises(ValueError):
        module.component_binding(review, energy, thermal)
