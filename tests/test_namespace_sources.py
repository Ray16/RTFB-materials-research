import hashlib
import json
from pathlib import Path

import pytest

from redox.validation.namespace_sources import NamespaceSources, read_checked_review


def setup(tmp_path):
    archive, attempt, raw = [tmp_path/name for name in ('archive', 'attempt', 'raw')]
    for path in (archive, attempt, raw):
        path.mkdir()
    return dict(prior_root='/old', archive=str(archive),
                original_campaign='/old/results/polymer_models/campaigns/batch/key', attempt=str(attempt),
                original_raw_parent='/old/calcs/polymer_models/batch', raw_parent=str(raw),
                original_raw_root='/old/calcs/polymer_models/batch/key', key='key', method_id='method')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_aliases_do_not_fall_back_to_historical_outputs(tmp_path):
    spec = setup(tmp_path)
    historical = Path(spec['archive'])/'calcs/polymer_models/batch/key/value'
    historical.parent.mkdir(parents=True)
    historical.write_text('previous failed attempt')
    reader = NamespaceSources(spec)
    assert reader.path('/old/calcs/polymer_models/batch/key/value') == Path(spec['raw_parent'])/'key/value'
    with pytest.raises(FileNotFoundError):
        reader.verify('/old/calcs/polymer_models/batch/key/value', digest(historical))


def test_raw_aliases_cannot_shadow_archived_sources(tmp_path):
    spec = setup(tmp_path)
    spec['raw_parent'] = str(Path(spec['archive'])/'mutable')
    Path(spec['raw_parent']).mkdir()
    with pytest.raises(ValueError):
        NamespaceSources(spec)


def energy_fixture(tmp_path):
    spec = setup(tmp_path)
    attempt, raw = Path(spec['attempt']), Path(spec['raw_parent'])/'key'
    raw.mkdir()
    task = dict(key='key', charge=0, mult=1)
    model = dict(solvent='acetonitrile')
    plan = attempt/'plan.json'
    plan.write_text(json.dumps(dict(tasks=[task], experimental_model=model)))
    energy = raw/'experimental_energy.json'
    energy.write_text(json.dumps(dict(e_smd_eV=-10., e_gas_eV=-9.)))
    report = dict(status='alternate_electronic_geometry_and_resource_checked_native_G_pending',
                  campaign=spec['original_campaign'], task=task, experimental_model=model,
                  portable_method_id='method', E_smd_eV=-10., E_gas_eV=-9.,
                  primary_candidate_admission=False, native_Hessian_or_G_qualified=False,
                  source_sha256={spec['original_campaign']+'/plan.json': digest(plan),
                                 spec['original_raw_root']+'/experimental_energy.json': digest(energy)})
    target = attempt/'independent_energy_review.json'
    target.write_text(json.dumps(report))
    path = tmp_path/'namespace.json'
    path.write_text(json.dumps(spec))
    return path, target, report


def test_existing_review_is_reused_only_after_all_raw_bindings_check(tmp_path):
    path, target, _ = energy_fixture(tmp_path)
    before = target.read_bytes()
    report, binding = read_checked_review(path, 'energy')
    assert report['E_smd_eV'] == -10.
    assert binding['verified_bindings'] == 2
    assert target.read_bytes() == before


@pytest.mark.parametrize('mutation', ['energy', 'method', 'scope', 'task'])
def test_wrong_reported_value_identity_or_scope_cannot_be_collected(tmp_path, mutation):
    path, target, report = energy_fixture(tmp_path)
    if mutation == 'energy': report['E_smd_eV'] = -11.
    if mutation == 'method': report['portable_method_id'] = 'other'
    if mutation == 'scope': report['native_Hessian_or_G_qualified'] = True
    if mutation == 'task': report['task']['charge'] = -1
    target.write_text(json.dumps(report))
    with pytest.raises(ValueError):
        read_checked_review(path, 'energy')


def test_missing_report_is_pending_without_creating_or_overwriting_one(tmp_path):
    path, target, _ = energy_fixture(tmp_path)
    target.unlink()
    with pytest.raises(FileNotFoundError):
        read_checked_review(path, 'energy', visibility_timeout=0)
    assert not target.exists()
