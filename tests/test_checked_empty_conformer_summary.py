import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/validation/summarize_checked_empty_conformers.py'
spec = importlib.util.spec_from_file_location('empty_summary', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def review(rank, energy):
    return {'status': 'alternate_native_state_independently_reconstructed',
            'primary_candidate_admission': False,
            'task': {'key': str(rank), 'occupied_sites': [], 'graph_sha256': 'graph',
                     'backbone_stereo_fingerprint': [[2, -1]], 'configuration_label': 'all_same',
                     'search_rank': rank},
            'experimental_model': {'solvent': 'acetonitrile'}, 'experimental_method_id': 'method',
            'G_eV': energy, 'G_by_cutoff_eV': {str(c): energy for c in (25, 50, 100)},
            'minimum_frequency_cm': 10.0}


def manifest(tmp_path, records):
    entries = []
    for i, record in enumerate(records):
        path = tmp_path / f'{i}.json'
        path.write_text(json.dumps(record))
        entries.append({'path': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    return {'reviews': entries}


def test_range_and_running_minimum_are_not_an_ensemble(tmp_path):
    result = module.summarize(manifest(tmp_path, [review(2, -10.2), review(1, -10.0)]), tmp_path)
    assert result['by_rotor_cutoff_cm']['50']['range_kJ_mol'] == pytest.approx(.2 * module.EV_KJ_MOL)
    assert result['by_rotor_cutoff_cm']['50']['running_minimum_eV'] == [-10.0, -10.2]
    assert not result['finite_ensemble_calculated'] and not result['wider_sampling_converged']


@pytest.mark.parametrize('mutation', ['method', 'graph', 'occupied', 'rank', 'nan', 'wrong_primary'])
def test_rejects_incompatible_or_invalid_review(tmp_path, mutation):
    first, second = review(1, -10), review(2, -10.2)
    if mutation == 'method': second['experimental_method_id'] = 'other'
    if mutation == 'graph': second['task']['graph_sha256'] = 'other'
    if mutation == 'occupied': second['task']['occupied_sites'] = [2]
    if mutation == 'rank': second['task']['search_rank'] = 1
    if mutation == 'nan': second['G_by_cutoff_eV']['50'] = float('nan')
    if mutation == 'wrong_primary': second['G_eV'] = -1
    with pytest.raises(ValueError):
        module.summarize(manifest(tmp_path, [first, second]), tmp_path)


def test_rejects_changed_review(tmp_path):
    m = manifest(tmp_path, [review(1, -10)])
    (tmp_path / '0.json').write_text('{}')
    with pytest.raises(ValueError, match='hash'):
        module.summarize(m, tmp_path)
