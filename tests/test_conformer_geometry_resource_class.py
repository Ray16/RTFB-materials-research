from copy import deepcopy

import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def case():
    return dict(smiles='C[C@H](Cl)CC', atoms=16, charge=0, mult=1, active_basis_functions=200)


def test_resource_class_can_cover_distinct_configurations_without_transferring_values():
    module = api.load_module('geometry_class_configuration_test', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    reference = case(); target = dict(reference, smiles='C[C@@H](Cl)CC')
    module.resource_compatible(target, reference, {'max_opt_steps': 400}, {'max_opt_steps': 400})
    assert target['smiles'] != reference['smiles']


@pytest.mark.parametrize('field,value', [('atoms', 17), ('charge', -1), ('mult', 3),
    ('active_basis_functions', 201), ('smiles', 'CCNCC')])
def test_geometry_resource_class_refuses_different_graph_state_or_basis(field, value):
    module = api.load_module('geometry_class_identity_test', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    reference = case(); target = deepcopy(reference); target[field] = value
    with pytest.raises(ValueError, match='resource class differs'):
        module.resource_compatible(target, reference, {}, {})


def test_geometry_resource_class_cannot_change_numerical_protocol():
    module = api.load_module('geometry_class_protocol_test', ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py')
    with pytest.raises(ValueError, match='resource class differs'):
        module.resource_compatible(case(), case(), {'SCF': 1e-8}, {'SCF': 1e-11})
