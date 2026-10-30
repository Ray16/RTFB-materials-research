import hashlib
import json
from pathlib import Path

import pytest

from redox.validation.archived_sources import ArchivedSources
from redox.validation.portable_energy_parent import relocate_plan


def fixture(tmp_path):
    archive = tmp_path/'archive'
    archive.mkdir()
    source = archive/'input.xyz'
    source.write_text('original coordinates')
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    plan = dict(tasks=[dict(key='exact', charge=0, mult=1)], experimental_method_id='method',
                quantum_settings=dict(convergence=1e-12), producer_sha256={'/old/input.xyz': digest},
                input_sha256={'/old/input.xyz': digest})
    path = tmp_path/'original_plan.json'
    path.write_text(json.dumps(plan))
    return plan, path, ArchivedSources('/old', archive)


def test_only_source_path_keys_change_and_original_bytes_remain(tmp_path):
    plan, path, reader = fixture(tmp_path)
    before = path.read_bytes()
    derived = relocate_plan(plan, reader, path)
    assert derived['tasks'] == plan['tasks']
    assert derived['quantum_settings'] == plan['quantum_settings']
    assert derived['input_sha256'] == {str(tmp_path/'archive/input.xyz'): next(iter(plan['input_sha256'].values()))}
    assert path.read_bytes() == before


def test_changed_numerical_plan_cannot_use_original_identity(tmp_path):
    plan, path, reader = fixture(tmp_path)
    plan['quantum_settings']['convergence'] = 1e-5
    with pytest.raises(ValueError, match='exact original'):
        relocate_plan(plan, reader, path)


def test_changed_backing_bytes_prevent_relocation(tmp_path):
    plan, path, reader = fixture(tmp_path)
    (tmp_path/'archive/input.xyz').write_text('changed coordinates')
    with pytest.raises(ValueError, match='differs'):
        relocate_plan(plan, reader, path)
