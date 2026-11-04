import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'manuscript/check_polymer_replay.py'
spec = importlib.util.spec_from_file_location('paper_polymer_replay', SOURCE)
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source = ROOT/'manuscript/evidence/2026-10-30_polymer_replay'
    for name in ('polymer_native_replay.json', 'polymer_native_replay_summary.json',
                 'polymer_native_sampling_figure.json', 'polymer_native_sampling.pdf',
                 'polymer_native_replay_tests.xml'):
        shutil.copyfile(source/name, tmp_path/name)
    return tmp_path


def test_retained_scientific_quantities_reconstruct(frozen):
    result = checker.audit(frozen)
    assert result['retained_states'] == 60
    assert result['empty_ensemble_change_kJ_mol']['50'] == pytest.approx(-1.0313810086)
    assert result['occupied_observed_minimum_change_kJ_mol'] == pytest.approx(-3.8437660762)
    assert result['complete_occupied_batch_convergence_proven'] is False


@pytest.mark.parametrize('mutation', ['energy', 'missing', 'nan_frequency', 'duplicate',
                                    'scope', 'figure', 'configuration'])
def test_manuscript_rejects_changed_values_or_scope(frozen, mutation):
    path = frozen/'polymer_native_replay.json'
    data = json.loads(path.read_text())
    if mutation == 'energy': data['records'][0]['G_by_cutoff_eV']['50'] += .1
    if mutation == 'missing': data['records'].pop()
    if mutation == 'nan_frequency': data['records'][0]['minimum_frequency_cm'] = float('nan')
    if mutation == 'duplicate': data['records'][1]['key'] = data['records'][0]['key']
    if mutation == 'scope': data['global_ensemble_convergence_proven'] = True
    if mutation == 'configuration': data['records'][0]['configuration'] = 'alternating'
    if mutation == 'figure': (frozen/'polymer_native_sampling.pdf').write_bytes(b'changed')
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        checker.audit(frozen)
