"""Relocated cached controls must reject altered chemistry and provenance."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

ROOT = Path(__file__).parents[1]
SOURCE = ROOT / 'examples/control_replay'
spec = importlib.util.spec_from_file_location('cached_controls_test', SOURCE / 'reproduce.py')
reader = importlib.util.module_from_spec(spec); spec.loader.exec_module(reader)


@pytest.fixture
def package(tmp_path):
    target = tmp_path / 'controls'
    shutil.copytree(SOURCE, target)
    return target


def save_manifest(package, manifest):
    (package / 'data/manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


def rehash(package, manifest, relative):
    manifest['files_sha256'][relative] = reader.digest(package / relative)
    save_manifest(package, manifest)


def test_full_cached_controls_reconstruct():
    result = reader.replay(SOURCE)
    assert result['quantum_states_reconstructed'] == 8
    assert result['Fc_functional_controls']['electronic_states'] == 6
    assert not result['polymer_G_qualified']


def test_electronic_replay_requires_no_site_packages():
    result = subprocess.run([sys.executable, '-I', '-S', str(SOURCE / 'reproduce.py'), '--electronic-only'],
                            capture_output=True, text=True, check=True)
    data = json.loads(result.stdout)
    assert data['quantum_states_reconstructed'] == 0
    assert data['Fc_functional_controls']['gaps_and_counterfactuals_reconstructed']
    assert not data['CUDA_or_xTB_loaded']


@pytest.mark.parametrize('change', ['checksum', 'scope', 'duplicate_state', 'escape',
                                  'displacement', 'unprobed_Hessian', 'borrowed_nuclear', 'symlink'])
def test_altered_cached_controls_fail_even_when_rehashed(package, change):
    manifest = reader.read_json(package / 'data/manifest.json')
    row = manifest['quantum_states'][0]
    if change == 'scope':
        manifest['polymer_free_energy_qualification'] = True
        save_manifest(package, manifest)
    elif change == 'duplicate_state':
        manifest['electronic_states'][0] = dict(manifest['electronic_states'][1])
        save_manifest(package, manifest)
    elif change == 'escape':
        manifest['files_sha256']['../outside.json'] = 'invalid'
        save_manifest(package, manifest)
    elif change == 'displacement':
        relative = row['evaluations'][1]
        point = reader.read_json(package / relative)
        point['coordinates_Bohr'][0][0] += .05
        (package / relative).write_text(json.dumps(point))
        rehash(package, manifest, relative)
    elif change == 'unprobed_Hessian':
        relative = row['files']['hessian']; path = package / relative
        matrix = np.load(path, allow_pickle=False)
        matrix[0, 0] += .1  # Outside the two finite-difference columns.
        np.save(path, matrix); rehash(package, manifest, relative)
    elif change == 'borrowed_nuclear':
        relative = row['files']['nuclear']
        other = next(r for r in manifest['quantum_states'] if r['species'] != row['species'])
        shutil.copyfile(package / other['files']['nuclear'], package / relative)
        rehash(package, manifest, relative)
    elif change == 'symlink':
        path = package / row['files']['geometry']; saved = path.with_suffix('.saved')
        path.rename(saved); path.symlink_to(saved.name)
    else:
        (package / row['files']['geometry']).write_text('changed raw bytes')
    with pytest.raises((ValueError, FileNotFoundError)):
        reader.replay(package)
