import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_current_paper_and_all_six_packages_are_declared():
    config = json.loads((ROOT/'config/manuscript_reviewer_bundle_v4.json').read_text())
    assert {p['id'] for p in config['packages']} == {
        'primary_review', 'first_site', 'numerical_controls', 'neighbor_cycle',
        'tempo_addendum', 'local_motif'}
    pointer = json.loads((ROOT/'manuscript/evidence/current.json').read_text())
    assert str(Path(pointer['manifest']).parent) in config['source_paths']
    assert pointer['figures'] in config['source_paths']
    assert config['data_license'] is None and config['author_metadata'] is None
