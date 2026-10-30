import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('checked_energy_native_preparation',
    ROOT/'scripts/pipeline/prepare_checked_energy_native_batch.py')
conveyor = importlib.util.module_from_spec(specification)
specification.loader.exec_module(conveyor)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixture(tmp_path, geometry='same coordinates', parent_energy=-100.):
    case, raw, archive = [tmp_path/name for name in ('case', 'raw', 'archive')]
    for path in (case, raw, archive):path.mkdir()
    (archive/'source.py').write_text('scientific source')
    (case/'input.xyz').write_text(geometry)
    parent = tmp_path/'electronic_review.json'
    parent.write_text(json.dumps(dict(E_smd_eV=parent_energy, E_gas_eV=-99.)))
    task = dict(key='exact', charge=0, mult=1)
    plan = dict(tasks=[task], experimental_model=dict(solvent='acetonitrile'), experimental_method_id='method',
        producer_sha256={'/old/source.py':digest(archive/'source.py')},
        input_sha256={str(parent):digest(parent)}, electronic_review=str(parent),
        scored_geometry_sha256=digest(case/'input.xyz'))
    (case/'plan.json').write_text(json.dumps(plan))
    spec = dict(key='exact', prior_root='/old', archive=str(archive), attempt=str(case),
        original_campaign='/old/results/polymer_models/campaigns/native/exact',
        original_raw_parent='/old/calcs/polymer_models/native', raw_parent=str(raw),
        plan_sha256=digest(case/'plan.json'), geometry_sha256=digest(case/'input.xyz'))
    path = tmp_path/'namespace.json'
    path.write_text(json.dumps(spec))
    report = dict(task=task, experimental_model=plan['experimental_model'], portable_method_id='method',
        campaign='/old/results/polymer_models/campaigns/energy/exact', E_smd_eV=-100., E_gas_eV=-99.,
        source_sha256={'/old/results/polymer_models/campaigns/energy/exact/input.xyz':
                       hashlib.sha256(b'same coordinates').hexdigest()})
    return path, report


def test_matching_existing_native_input_is_reused_without_rewriting(tmp_path):
    namespace, report = fixture(tmp_path)
    before = namespace.read_bytes()
    assert conveyor.match_native(namespace, report)['key'] == 'exact'
    assert namespace.read_bytes() == before


@pytest.mark.parametrize('geometry,energy', [('another pose', -100.), ('same coordinates', -99.5)])
def test_rehashed_wrong_pose_or_electronic_parent_cannot_be_reused(tmp_path, geometry, energy):
    namespace, report = fixture(tmp_path, geometry, energy)
    with pytest.raises(ValueError):
        conveyor.match_native(namespace, report)
