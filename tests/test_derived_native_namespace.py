import copy
import hashlib
import json
from pathlib import Path

import pytest

from redox.workflows.derived_native_namespace import namespace_command, portable_parent


def fixture(tmp_path):
    archive, attempt, raw = (tmp_path/name for name in ('archive', 'energy', 'raw'))
    for path in (archive, attempt, raw):
        path.mkdir()
    value = raw/'value.json'
    value.write_text('{"electronic": -10.0}')
    digest = hashlib.sha256(value.read_bytes()).hexdigest()
    spec = dict(archive=str(archive), prior_root='/old',
                original_campaign='/old/results/polymer_models/campaigns/energy/key', attempt=str(attempt),
                original_raw_parent='/old/calcs/polymer_models/energy', raw_parent=str(raw))
    report = dict(status='independently_checked', campaign=spec['original_campaign'],
                  task=dict(key='key', charge=0, mult=1), E_smd_eV=-10., primary_candidate_admission=False,
                  source_sha256={'/old/calcs/polymer_models/energy/value.json': digest})
    original = attempt/'independent_energy_review.json'
    original.write_text(json.dumps(report))
    return spec, report, original, value


def test_parent_relocation_changes_paths_only_and_preserves_original_report(tmp_path):
    spec, report, original, value = fixture(tmp_path)
    before = original.read_bytes()
    result = portable_parent(report, spec, original)
    assert result['campaign'] == spec['attempt']
    assert result['E_smd_eV'] == report['E_smd_eV']
    assert result['task'] == report['task']
    assert str(value) in result['source_sha256']
    assert str(original) in result['source_sha256']
    assert result['path_relocation']['numerical_or_chemical_fields_changed'] is False
    assert original.read_bytes() == before


def test_modified_parent_number_cannot_claim_the_original_review(tmp_path):
    spec, report, original, _ = fixture(tmp_path)
    report['E_smd_eV'] = -11.
    with pytest.raises(ValueError, match='exact independently checked report'):
        portable_parent(report, spec, original)


def test_changed_backing_energy_prevents_native_handoff(tmp_path):
    spec, report, original, value = fixture(tmp_path)
    value.write_text('changed')
    with pytest.raises(ValueError, match='differs'):
        portable_parent(report, spec, original)


def test_new_batches_use_private_parents_with_old_batches_read_only(tmp_path):
    archive = tmp_path/'archive'
    for parent in ('results/polymer_models/campaigns', 'calcs/polymer_models'):
        (archive/parent/'previous').mkdir(parents=True)
    spec = dict(prior_root='/old', archive=str(archive),
                original_campaign='/old/results/polymer_models/campaigns/new/key',
                original_raw_parent='/old/calcs/polymer_models/new',
                attempt=str(tmp_path/'new/campaign_parent/key'), raw_parent=str(tmp_path/'new/raw'))
    result = namespace_command(spec)
    assert result.count('--tmpfs') == 3
    assert result.count('--bind') == 2
    assert result.count('--ro-bind') == 4
    spec['attempt'] = str(archive/'results/key')
    with pytest.raises(ValueError, match='outside'):
        namespace_command(spec)
