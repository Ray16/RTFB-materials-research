import importlib.util
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts/validation'))
spec = importlib.util.spec_from_file_location('empty_batch', ROOT/'scripts/validation/compare_empty_conformer_batches.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(n):
    from redox.properties.conformer_ensemble import summarize
    keys = [str(i) for i in range(1, n+1)]
    rows = [dict(key=k, E_eV=-10., G_eV=-9.) for k in keys]
    part = summarize(rows, [[k] for k in keys], keys)
    finite = dict(keys=keys, partitions_by_rmsd_cutoff_A={
        rmsd: {rotor: part for rotor in ('25', '50', '100')} for rmsd in ('0.05', '0.1', '0.2')})
    review = dict(configuration='one', graph_sha256='graph', method_id='method', solvent='acetonitrile',
                  records=[dict(key=k, search_rank=int(k)) for k in keys])
    return finite, review


def test_batch_entropy_is_tested_even_when_minima_are_unchanged():
    first, fr = fixture(10)
    expanded, er = fixture(20)
    result = module.compare(first, expanded, fr, er)
    row = result['by_rmsd_A_and_rotor_cm']['0.1']['50']
    assert row['minimum_E_change_kJ_mol'] == 0
    assert row['ensemble_G_change_kJ_mol'] == pytest.approx(-1.7182820758)
    assert row['within_batch_tolerance'] is True
    assert result['global_ensemble_convergence_proven'] is False


@pytest.mark.parametrize('mutation', ['method', 'rank', 'count'])
def test_rejects_incompatible_or_postselected_batches(mutation):
    first, fr = fixture(10)
    expanded, er = fixture(19 if mutation == 'count' else 20)
    if mutation == 'method': er['method_id'] = 'another'
    if mutation == 'rank': er['records'][-1]['search_rank'] = 21
    with pytest.raises(ValueError):
        module.compare(first, expanded, fr, er)


def test_rejects_changed_stored_partition(tmp_path, monkeypatch):
    saved = tmp_path/'finite.json'
    saved.write_text('{"G_eV": -1}')
    monkeypatch.setattr(module, 'analyze', lambda *args: dict(G_eV=-2))
    with pytest.raises(ValueError, match='differs'):
        module.checked_batch(ROOT, tmp_path/'manifest.json', saved)
