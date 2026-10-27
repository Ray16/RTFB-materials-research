import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('remaining_terms',ROOT/'scripts/validation/check_counterpoise_remaining_terms.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
BASE=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'


def test_rejects_wrong_role_checkpoint_before_energy_evaluation():
    data=json.loads((BASE/'summary.json').read_text())
    geometry=dict(symbols=data['source_symbols'],coordinates_A=data['source_coordinates_A'],A_indices=data['A_indices'])
    with pytest.raises(ValueError,match='Checkpoint nuclei'):
        module.evaluate(geometry,'AB',BASE/'B_own/SCF.h5',BASE/'grid_coordinates_AU.npy',BASE/'grid_weights.npy',0.)


def test_unplanned_role_is_rejected():
    with pytest.raises(ValueError,match='Unsupported'):
        module.evaluate({},'A_own',Path('unused'),Path('unused'),Path('unused'),0.)
