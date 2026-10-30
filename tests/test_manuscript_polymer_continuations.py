import importlib.util
import hashlib
import json
from pathlib import Path
import shutil

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('current_polymer_claims', ROOT/'manuscript/check_polymer_continuations.py')
checker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checker)


@pytest.fixture
def frozen(tmp_path):
    source = ROOT/'manuscript/evidence/2026-10-30_polymer_continuations'
    for name in ('polymer_native_replay.json', 'polymer_native_continuations.json',
                 'polymer_native_sampling_current.json', 'polymer_native_sampling_current.pdf',
                 'nitro_repeatability.json'):
        shutil.copyfile(source/name, tmp_path/name)
    return tmp_path


def test_current_coverage_and_backend_limitation_reconstruct(frozen):
    result = checker.audit(frozen)
    assert result['checked_states'] == 67
    assert result['complete_occupied_batch_convergence_proven'] is False
    assert result['nitro_production_qualified'] is False


@pytest.mark.parametrize('mutation', ['energy', 'duplicate', 'model', 'nonfinite',
                                    'convergence', 'figure', 'panel', 'nitro'])
def test_changed_sampling_or_backend_claim_cannot_pass(frozen, mutation):
    path = frozen/'polymer_native_continuations.json'
    data = json.loads(path.read_text())
    if mutation == 'energy': data['records'][0]['G_by_cutoff_eV']['50'] += .1
    if mutation == 'duplicate': data['records'][0]['key'] = data['records'][1]['key']
    if mutation == 'model': data['records'][0]['method_id'] = 'a different solvent model'
    if mutation == 'nonfinite': data['records'][0]['minimum_frequency_cm'] = float('nan')
    if mutation == 'convergence': data['global_ensemble_convergence_proven'] = True
    if mutation == 'figure': (frozen/'polymer_native_sampling_current.pdf').write_bytes(b'changed')
    if mutation == 'panel':
        p = frozen/'polymer_native_sampling_current.json'
        fig = json.loads(p.read_text())
        fig['panels'][0] = fig['panels'][1]
        p.write_text(json.dumps(fig))
    if mutation == 'nitro':
        p = frozen/'nitro_repeatability.json'
        backend = json.loads(p.read_text())
        backend['production_calculator_qualified'] = True
        p.write_text(json.dumps(backend))
    path.write_text(json.dumps(data))
    # Rebind bytes so each negative check must reject the changed scientific
    # content, rather than passing merely because JSON formatting changed.
    figure_path = frozen/'polymer_native_sampling_current.json'
    figure = json.loads(figure_path.read_text())
    figure['continuation_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    figure_path.write_text(json.dumps(figure))
    with pytest.raises(ValueError):
        checker.audit(frozen)
