import importlib.util
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('primary_relocation_test',ROOT/'scripts/release/validate_primary_relocation.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_replay_versions_must_match_declared_lock():
    expected={'numpy':'2.2.6','scipy':'1.17.1','ase':'3.29.0'}
    assert module.validate_versions(expected)==expected
    with pytest.raises(ValueError,match='dependencies differ'):
        module.validate_versions(dict(expected,ase='0.0'))
