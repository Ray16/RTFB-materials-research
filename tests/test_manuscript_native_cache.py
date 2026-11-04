import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[1]
definition = importlib.util.spec_from_file_location('paper_portable_cache', ROOT/'manuscript/check_native_cache.py')
checker = importlib.util.module_from_spec(definition)
definition.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source = ROOT/'manuscript/evidence/2026-10-30_empty_backbones'
    for name in ('polymer_native_replay.json', 'polymer_native_continuations.json',
                 'native_addition_five.json', 'native_addition_three.json'):
        shutil.copyfile(source/name, tmp_path/name)
    cache = ROOT/'validation/native_cache_2026-10-30'
    for a,b in [('replay.json','replay'),('isolation.json','isolation'),
                ('compression.json','compression'),('package_manifest.json','manifest')]:
        shutil.copyfile(cache/a, tmp_path/('portable_native_cache_'+b+'.json'))
    return tmp_path


def test_archive_is_tied_to_the_same_75_frozen_native_states(frozen):
    result = checker.audit(frozen)
    assert result['states'] == 75 and result['extracted_archive_fresh_numpy_replay']
    assert not result['full_historical_source_closure_supplied']
    assert not result['public_archive_deposited']


@pytest.mark.parametrize('mutation', ['energy','duplicate','nonfinite','scope','columns','archive','installer','maximum'])
def test_changed_raw_quantities_coverage_and_release_claims_cannot_pass(frozen, mutation):
    name = 'portable_native_cache_replay.json'
    if mutation in ('archive','installer','maximum'):name = 'portable_native_cache_isolation.json'
    path = frozen/name;data = json.loads(path.read_text())
    if mutation == 'energy':data['records'][0]['G_by_cutoff_eV']['50'] += .01
    if mutation == 'duplicate':data['records'][1] = data['records'][0]
    if mutation == 'nonfinite':data['records'][0]['nuclear_errors_eV']['50'] = float('nan')
    if mutation == 'scope':data['full_historical_source_closure_revalidated'] = True
    if mutation == 'columns':data['records'][0]['sampled_columns'] += 1
    if mutation == 'archive':data['archive_sha256'] = 'different'
    if mutation == 'installer':data['original_project_installed'] = True
    if mutation == 'maximum':data['maximum_frequency_error_cm'] = .001
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker.audit(frozen)
