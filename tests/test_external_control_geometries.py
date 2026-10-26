import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


source = Path(__file__).resolve().parents[1] / 'scripts/validation/prepare_external_control_geometries.py'
spec = importlib.util.spec_from_file_location('external_control_seeds', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def challenge(tmp_path):
    mol = Chem.MolFromSmiles('c1ccccc1')
    path = tmp_path / 'challenge.json'
    path.write_text(json.dumps({'status': 'source_checked_NDI_external_challenge_not_scored',
        'computed_potentials': None, 'molecule': {'name': 'fixture', 'smiles': 'c1ccccc1', 'inchikey': Chem.MolToInchiKey(mol)},
        'proposed_states': [{'label': 'neutral', 'charge': 0, 'multiplicity': 1},
                            {'label': 'reduced_1', 'charge': -1, 'multiplicity': 2},
                            {'label': 'reduced_2', 'charge': -2, 'multiplicity': 1}]}))
    return path


def test_small_control_seeds_are_graph_checked_and_unqualified(tmp_path):
    result = module.prepare(challenge(tmp_path), tmp_path/'output', seed=23, attempts=3, top_k=2)
    assert result['selected'] >= 1
    assert result['quantum_jobs_submitted'] == 0 and not result['charged_state_minima_computed']
    assert (tmp_path/'output'/'manifest.csv').is_file()
    assert all(e['topology'] == 'ok' and (tmp_path/'output'/e['path']).is_file() for e in result['entries'])


def test_invalid_sampling_and_mismatched_identity_fail_before_output(tmp_path):
    path = challenge(tmp_path)
    with pytest.raises(ValueError, match='sampling'):
        module.prepare(path, tmp_path/'invalid', attempts=2, top_k=3)
    record = json.loads(path.read_text()); record['molecule']['inchikey'] = 'WRONG'; path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='identity'):
        module.prepare(path, tmp_path/'wrong', attempts=2, top_k=1)
    assert not (tmp_path/'invalid').exists() and not (tmp_path/'wrong').exists()


def test_independent_reviewer_rejects_changed_coordinates(tmp_path):
    path = source.with_name('review_external_control_geometries.py')
    review_spec = importlib.util.spec_from_file_location('external_control_review', path)
    reviewer = importlib.util.module_from_spec(review_spec)
    review_spec.loader.exec_module(reviewer)
    directory = tmp_path / 'reviewable'
    result = module.prepare(challenge(tmp_path), directory, seed=31, attempts=3, top_k=2)
    assert reviewer.review(directory)['selected'] == result['selected']
    xyz = directory / result['entries'][0]['path']
    xyz.write_text(xyz.read_text() + '\n')
    with pytest.raises(ValueError, match='hash changed'):
        reviewer.review(directory)
