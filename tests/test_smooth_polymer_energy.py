"""Scientific scope gates for the alternate-model electronic stage."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('smooth_energy_test', ROOT / 'scripts/pipeline/smooth_polymer_energy.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture
def energy():
    return dict(charge=0, mult=1, solvent='acetonitrile', optimized=False,
        sp_xc='wb97m-v', sp_basis='def2-tzvpd', converged_smd=True, converged_gas=True,
        experimental_method_id='alternate-model', primary_candidate_admission=False,
        native_G_computed=False, e_smd_eV=-1000., e_gas_eV=-999.,
        smd_spin_contam=0., gas_spin_contam=0.)


def test_electronic_record_retains_alternate_scope(energy):
    module.check_energy(energy, dict(charge=0, mult=1), 'alternate-model')


@pytest.mark.parametrize('field,value', [('native_G_computed', True), ('primary_candidate_admission', True),
    ('G_eV', -1000.), ('g_thermal_eV', 0.), ('experimental_method_id', 'legacy-model'),
    ('converged_smd', False), ('charge', -1), ('smd_spin_contam', .11), ('e_smd_eV', float('nan'))])
def test_unqualified_electronic_state_cannot_be_promoted(energy, field, value):
    row = deepcopy(energy); row[field] = value
    with pytest.raises(ValueError):
        module.check_energy(row, dict(charge=0, mult=1), 'alternate-model')
