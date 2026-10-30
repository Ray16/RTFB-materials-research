import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('measured_energy_followup_gate',
    ROOT/'scripts/pipeline/run_measured_energy_followup.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path, *, peak=100, witness_atoms=105):
    archive = tmp_path/'archive'
    archive.mkdir()
    source = archive/'source.py'
    source.write_text('source')
    sources = {'/old/source.py': hashlib.sha256(source.read_bytes()).hexdigest()}
    task = dict(atoms=105, active_basis_functions=2537, charge=0, mult=1, smiles='CCO',
                candidate_id='exp_phthalimide_control', role='interaction_method_control; not synthesis admission')
    witness = dict(status='alternate_electronic_geometry_and_resource_checked_native_G_pending',
                   native_Hessian_or_G_qualified=False, primary_candidate_admission=False,
                   portable_method_id='model', experimental_model=dict(solvent='acetonitrile'),
                   campaign='/old/reference', peak_host_RSS_bytes=peak*2**30,
                   observed_peak_GPU_MiB=20000, task=dict(task, atoms=witness_atoms), source_sha256=sources)
    path = archive/'review.json'
    path.write_text(json.dumps(witness))
    plan = dict(tasks=[task], maximum_host_RSS_GiB=160, prospective_resource_pilot=False,
                execution_driver='/old/scripts/pipeline/smooth_polymer_energy_bounded.py',
                native_G_computed=False, primary_candidate_admission=False, energy_only_execution=True,
                maximum_peak_GPU_MiB=31744, experimental_model=dict(solvent='acetonitrile'),
                experimental_method_id='model', electronic_RAM_class_review='/old/review.json',
                measured_alternate_model_electronic_reference='/old/reference', producer_sha256=sources,
                input_sha256={'/old/review.json': hashlib.sha256(path.read_bytes()).hexdigest()})
    return plan, archive


def test_measured_class_transfers_only_resource_budget(tmp_path):
    plan, archive = fixture(tmp_path)
    result = module.validate_measured(plan, Path('/old'), archive)
    assert result['bound_GiB'] == 160 and result['measured_peak_GiB'] == 100
    assert result['resource_budget_only_no_energy_transferred'] is True
    assert plan['prospective_resource_pilot'] is False


@pytest.mark.parametrize('peak,witness_atoms', [(129, 105), (100, 119)])
def test_insufficient_headroom_or_wrong_constitution_class_is_rejected(tmp_path, peak, witness_atoms):
    plan, archive = fixture(tmp_path, peak=peak, witness_atoms=witness_atoms)
    with pytest.raises(ValueError):
        module.validate_measured(plan, Path('/old'), archive)


def test_changed_reference_bytes_cannot_authorize_successor(tmp_path):
    plan, archive = fixture(tmp_path)
    (archive/'review.json').write_text('{}')
    with pytest.raises(ValueError, match='differs'):
        module.validate_measured(plan, Path('/old'), archive)
