"""Exercise shared molecular checkpoint admission using actual CPU wavefunctions."""
import json
from pathlib import Path

from ase import Atoms
from ase.io import write
import pytest
from pyscf import gto, scf, grad, hessian

from redox.qm.scf_checkpoints import recording_builder
from redox.validation.primary_checkpoints import checkpoint_series


@pytest.fixture
def actual_checkpoint(tmp_path):
    molecule = gto.M(atom='H 0 0 0; H 0 0 .74', basis='sto-3g', verbose=0)
    root = tmp_path/'SCF_evaluations'
    mf = recording_builder(lambda: scf.RHF(molecule), root, {'test': 'actual_cpu_wavefunction'})()
    mf.kernel()
    geometry = tmp_path/'geometry.xyz'
    write(geometry, Atoms('H2', positions=molecule.atom_coords(unit='Angstrom')), format='xyz')
    return root, geometry, dict(charge=0, mult=1, basis='sto-3g', solvent_class=None)


def read(path):
    return json.loads(Path(path).read_text())


def test_actual_wavefunction_energy_nelec_spin_and_fixed_frame(actual_checkpoint):
    root, geometry, settings = actual_checkpoint
    result = checkpoint_series(root, geometry, settings, read, [], expected_count=1, fixed_geometry=True)
    assert result['actual_SCF_frames_checked'] == 1
    assert result['converged_energies_Ha'] == result['energies_Ha']
    assert abs(result['spins'][0]['spin_contamination']) < 1e-10
    assert result['spins'][0]['alpha_electrons'] == result['spins'][0]['beta_electrons'] == 1


@pytest.mark.parametrize('mutation', ['basis', 'charge', 'multiplicity', 'frame', 'energy', 'count'])
def test_swapped_state_frame_and_scalar_cannot_pass_actual_checkpoint_admission(actual_checkpoint, mutation):
    root, geometry, settings = actual_checkpoint
    count = 1
    if mutation == 'basis': settings['basis'] = 'def2-svp'
    if mutation == 'charge': settings['charge'] = 1
    if mutation == 'multiplicity': settings['mult'] = 3
    if mutation == 'frame': write(geometry, Atoms('H2', positions=[[0,0,0],[0,0,.8]]), format='xyz')
    if mutation == 'energy':
        path = next(root.glob('*/*/SCF_result.json'))
        result = read(path); result['energy_Ha'] += .01; path.write_text(json.dumps(result))
    if mutation == 'count': count = 9
    with pytest.raises(ValueError):
        checkpoint_series(root, geometry, settings, read, [], expected_count=count, fixed_geometry=True)


def test_unconverged_search_step_remains_unqualified(actual_checkpoint):
    root, geometry, settings = actual_checkpoint
    path = next(root.glob('*/*/SCF_result.json'))
    result = read(path); result['converged'] = False; path.write_text(json.dumps(result))
    with pytest.raises(ValueError): checkpoint_series(root, geometry, settings, read, [])
    report = checkpoint_series(root, geometry, settings, read, [], require_all_converged=False)
    assert report['converged_energies_Ha'] == []
    assert len(report['energies_Ha']) == 1
