import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('stage_archived_energy', ROOT/'scripts/pipeline/stage_checked_archived_energy.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path):
    archive = tmp_path/'archive'
    case = archive/'results/polymer_models/campaigns/batch/key'
    raw = archive/'calcs/polymer_models/batch/key'
    case.mkdir(parents=True)
    raw.mkdir(parents=True)
    prior = Path('/old')
    original_case = prior/'results/polymer_models/campaigns/batch/key'
    original_raw = prior/'calcs/polymer_models/batch/key'
    task, model = dict(key='key', charge=0, mult=1), dict(solvent='acetonitrile')
    (case/'input.xyz').write_text('checked input bytes')
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    (case/'plan.json').write_text(json.dumps(dict(tasks=[task], experimental_model=model,
        raw_root=str(original_raw), scored_geometry_sha256=sha(case/'input.xyz'))))
    (raw/'experimental_energy.json').write_text(json.dumps(dict(e_smd_eV=-10., e_gas_eV=-9.)))
    bindings = {str(original_case/name): sha(case/name) for name in ('plan.json', 'input.xyz')}
    bindings[str(original_raw/'experimental_energy.json')] = sha(raw/'experimental_energy.json')
    review = tmp_path/'review.json'
    review.write_text(json.dumps(dict(status='alternate_electronic_geometry_and_resource_checked_native_G_pending',
        primary_candidate_admission=False, native_Hessian_or_G_qualified=False, campaign=str(original_case),
        task=task, experimental_model=model, portable_method_id='model', E_smd_eV=-10., E_gas_eV=-9.,
        source_sha256=bindings)))
    return archive, prior, case, raw, review


def test_reviewed_evidence_is_staged_byte_identically_without_new_quantum_work(tmp_path):
    archive, prior, case, raw, review = fixture(tmp_path)
    destination = tmp_path/'staged'
    result = module.stage(archive, prior, review, destination)
    assert result['quantum_calculations_run'] is False
    assert (destination/'campaign/plan.json').read_bytes() == (case/'plan.json').read_bytes()
    assert (destination/'raw/key/experimental_energy.json').read_bytes() == (raw/'experimental_energy.json').read_bytes()
    assert (destination/'campaign/independent_energy_review.json').read_bytes() == review.read_bytes()


def test_changed_raw_values_cannot_be_staged_under_old_review(tmp_path):
    archive, prior, _, raw, review = fixture(tmp_path)
    (raw/'experimental_energy.json').write_text('{}')
    with pytest.raises(ValueError, match='differs'):
        module.stage(archive, prior, review, tmp_path/'staged')


def test_repeated_staging_preserves_previous_destination(tmp_path):
    archive, prior, _, _, review = fixture(tmp_path)
    destination = tmp_path/'staged'
    destination.mkdir()
    with pytest.raises(ValueError, match='fresh staged-evidence'):
        module.stage(archive, prior, review, destination)
