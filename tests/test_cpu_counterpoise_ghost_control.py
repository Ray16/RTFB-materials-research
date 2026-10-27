import importlib.util
import json
from pathlib import Path

import pytest
from pyscf.scf.chkfile import load_scf


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_ghost_control',ROOT/'scripts/validation/check_cpu_counterpoise_ghost_control.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
BASE=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'


def test_rejects_real_only_checkpoint_for_ghost_calculation():
    core=module.api.load_module('test_ghost_counterpoise_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    record=json.loads((BASE/'summary.json').read_text())
    mols=core.molecules(record['source_symbols'],record['source_coordinates_A'],record['A_indices'])
    ghost,_=load_scf(str(BASE/'A_ghost/SCF.h5'))
    own,_=load_scf(str(BASE/'A_own/SCF.h5'))
    module.verify_ghost_checkpoint(ghost,mols['A_ghost'],mols['AB'],3)
    with pytest.raises(ValueError,match='Ghost checkpoint'):
        module.verify_ghost_checkpoint(own,mols['A_ghost'],mols['AB'],3)
