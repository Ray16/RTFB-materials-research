import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('paper_empty_batch_claims', ROOT/'manuscript/check_empty_batches.py')
checker = importlib.util.module_from_spec(specification)
specification.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    inherited = ROOT/'manuscript/evidence/2026-10-30_polymer_continuations'
    for name in ('polymer_native_replay.json', 'polymer_native_continuations.json'):
        shutil.copyfile(inherited/name, tmp_path/name)
    additions = ROOT/'validation/polymer_empty_batches_2026-10-30'
    files = [(ROOT/'validation/native_working_collection_2026-10-30/replay.json', 'native_addition_five.json'),
             (additions/'native_additions.json', 'native_addition_three.json'),
             (additions/'analysis.json', 'empty_batch_analysis.json'),
             (additions/'figure.json', 'empty_batch_figure.json'),
             (additions/'figure.pdf', 'empty_batch_stability.pdf'),
             (additions/'values.tex', 'empty_batch_values.tex')]
    for source, target in files:
        shutil.copyfile(source, tmp_path/target)
    return tmp_path


def test_two_configuration_scientific_quantities_reconstruct(frozen):
    result = checker.audit(frozen)
    assert result['empty_states'] == 40
    assert result['unique_minima'] == {'all_same': 20, 'alternating': 19}
    assert result['nominal_ensemble_change_kJ_mol']['alternating'] == pytest.approx(-1.2311873927)
    assert result['occupied_ensemble_cycles_complete'] is False


@pytest.mark.parametrize('mutation', ['energy', 'duplicate', 'matrix', 'population', 'scope', 'macro', 'figure'])
def test_changed_values_or_inferences_cannot_pass(frozen, mutation):
    path = frozen/'empty_batch_analysis.json'
    data = json.loads(path.read_text())
    if mutation == 'energy':
        p = frozen/'polymer_native_replay.json'
        native = json.loads(p.read_text())
        native['records'][0]['G_by_cutoff_eV']['50'] += .1
        p.write_text(json.dumps(native))
    if mutation == 'duplicate':
        data['groups'][1]['partitions']['0.1']['50']['expanded']['duplicate_groups'].append(['not_a_state'])
    if mutation == 'matrix':
        data['groups'][1]['rmsd_A'][15][16] = data['groups'][1]['rmsd_A'][16][15] = .3
    if mutation == 'population':
        data['groups'][0]['partitions']['0.1']['50']['expanded']['conditional_populations'][0] += .01
    if mutation == 'scope': data['configurations_are_pooled'] = True
    if mutation == 'macro': (frozen/'empty_batch_values.tex').write_text('invented numbers')
    if mutation == 'figure': (frozen/'empty_batch_stability.pdf').write_bytes(b'changed')
    path.write_text(json.dumps(data))
    figure_path = frozen/'empty_batch_figure.json'
    figure = json.loads(figure_path.read_text())
    figure['analysis_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    figure_path.write_text(json.dumps(figure))
    with pytest.raises(ValueError):
        checker.audit(frozen)
