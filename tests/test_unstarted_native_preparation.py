import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('unstarted_native_preparation',
    ROOT/'scripts/pipeline/prepare_unstarted_archived_native.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def inputs(tmp_path):
    archive = tmp_path/'archive'
    source = archive/'results/polymer_models/campaigns/batch/key'
    source.mkdir(parents=True)
    plan = dict(tasks=[dict(key='key', charge=0, mult=1)],
                experimental_model=dict(solvent='acetonitrile'),
                execution_driver='/old/scripts/pipeline/smooth_polymer_native.py',
                maximum_host_RSS_GiB=80, maximum_peak_GPU_MiB=31744,
                raw_root='/old/calcs/polymer_models/batch/key', electronic_review='/old/parent.json',
                input_sha256={}, producer_sha256={})
    (source/'plan.json').write_text(json.dumps(plan))
    inventory = tmp_path/'inventory.json'
    inventory.write_text(json.dumps(dict(keys=['key'], nodes=[dict(host=str(i), processes=[]) for i in range(9)])))
    args = [archive, Path('/old'), Path('/old/results/polymer_models/campaigns/batch/key'),
            tmp_path/'new', inventory]
    return source, plan, args


@pytest.mark.parametrize('existing', ['stage', 'raw'])
def test_existing_scientific_attempt_cannot_be_reclassified_as_unstarted(tmp_path, existing):
    source, _, args = inputs(tmp_path)
    if existing == 'stage':
        (source/'experimental_native').mkdir()
    else:
        (args[0]/'calcs/polymer_models/batch/key').mkdir(parents=True)
    with pytest.raises(ValueError, match='existing native'):
        module.prepare(*args)
    assert not args[3].exists()


def test_live_or_unknown_predictor_keeps_preparation_closed(tmp_path):
    _, _, args = inputs(tmp_path)
    inventory = json.loads(args[-1].read_text())
    inventory['nodes'][0]['processes'] = [dict(pid=123, observed='live')]
    args[-1].write_text(json.dumps(inventory))
    with pytest.raises(ValueError, match='Stopped-predictor'):
        module.prepare(*args)


def test_missing_electronic_binding_cannot_create_a_native_plan(tmp_path):
    _, _, args = inputs(tmp_path)
    with pytest.raises(ValueError, match='electronic review is not bound'):
        module.prepare(*args)


def test_different_multiplicity_is_outside_the_validated_resource_class(tmp_path):
    source, plan, args = inputs(tmp_path)
    plan['tasks'][0]['mult'] = 3
    (source/'plan.json').write_text(json.dumps(plan))
    with pytest.raises(ValueError, match='neutral-singlet'):
        module.prepare(*args)
