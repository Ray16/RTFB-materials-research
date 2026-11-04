"""Runtime discovery must not hide scientific source or certificate changes."""
import hashlib
import json
from pathlib import Path

import pytest

from redox.qm.thermal_source import matches_qualified_source

ROOT = Path(__file__).resolve().parents[1]


def test_runtime_only_change_reconstructs_the_original_qualified_hash():
    proof = json.loads((ROOT/'config/projected_thermal_validation.json').read_text())
    source = (ROOT/'src/redox/qm/qrrho.py').read_bytes()
    assert matches_qualified_source(source, proof['qrrho_source_sha256'])
    assert matches_qualified_source(source, hashlib.sha256(source).hexdigest())


@pytest.mark.parametrize('old,new', [
    (b'HARTREE_EV = 27.211386245988', b'HARTREE_EV = 27.2'),
    (b'rotor_cutoffs_cm=[25, 50, 100]', b'rotor_cutoffs_cm=[20, 50, 100]'),
    (b'configured_library()', b'other_library()'),
    (b'import numpy as np', b'import numpy as np\nimport extra_dependency'),
])
def test_any_additional_source_change_is_rejected(old, new):
    proof = json.loads((ROOT/'config/projected_thermal_validation.json').read_text())
    source = (ROOT/'src/redox/qm/qrrho.py').read_bytes()
    assert old in source
    assert not matches_qualified_source(source.replace(old, new, 1), proof['qrrho_source_sha256'])


def test_original_certificate_and_actual_runtime_are_still_required(monkeypatch):
    from redox.qm import projected_thermal as module
    assert module.certificate()
    monkeypatch.setattr(module.Q, 'LIBRARY', ROOT/'absent-library')
    with pytest.raises((OSError, RuntimeError)):
        module.certificate()
