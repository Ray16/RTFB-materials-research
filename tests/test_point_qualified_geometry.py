import copy
import json
from pathlib import Path

import pytest

from redox.core.common import ROOT
from redox.core.protocol import ACTIVE_SP, geom_sha1
from redox.workflows import polymer_campaign as api


@pytest.fixture
def point(tmp_path):
    pilot, review = tmp_path / 'point', tmp_path / 'review'
    xyz = pilot / 'live_inputs/graph/input.xyz'
    xyz.parent.mkdir(parents=True)
    xyz.write_text('2\ninput\nH 0 0 0\nH 0 0 0.74\n')
    review.mkdir()
    task = dict(key='graph', charge=0, mult=1)
    plan = dict(tasks=[task], producer_sha256={}, input_sha256={})
    (pilot / 'plan.json').write_text(json.dumps(plan))
    folder = pilot / 'energy/pilot'
    folder.mkdir(parents=True)
    record = dict(result=dict(xyz=str(xyz), energy_geom_sha1=geom_sha1(xyz), energy_protocol=ACTIVE_SP))
    (folder / 'energy_and_gradient.json').write_text(json.dumps(record))
    proof = dict(status='electronic_resource_point_independently_reviewed', source=str(pilot),
        task_key='graph', G_computed=False, optimized_geometry_qualified=False, source_sha256={})
    (review / 'summary.json').write_text(json.dumps(proof))
    module = api.load_module('point_geometry_fixture', ROOT / 'scripts/pipeline/run_point_qualified_geometry.py')
    return module, pilot, review


def test_exact_point_remains_a_starting_geometry(point):
    module, pilot, review = point
    plan, proof, xyz = module.point_evidence(pilot, review)
    assert plan['tasks'][0]['key'] == 'graph' and xyz.is_file()
    assert not proof['G_computed'] and not proof['optimized_geometry_qualified']


@pytest.mark.parametrize('change', ['charge', 'key', 'coordinates', 'protocol', 'false_minimum'])
def test_point_cannot_qualify_different_state_geometry_protocol_or_scope(point, change):
    module, pilot, review = point
    plan_path, record_path, proof_path = pilot / 'plan.json', pilot / 'energy/pilot/energy_and_gradient.json', review / 'summary.json'
    if change in ('charge', 'key'):
        value = json.loads(plan_path.read_text())
        value['tasks'][0][change] = -1 if change == 'charge' else 'other_graph'
        plan_path.write_text(json.dumps(value))
    elif change in ('coordinates', 'protocol'):
        value = json.loads(record_path.read_text())
        if change == 'coordinates':
            altered = pilot / 'altered.xyz'
            altered.write_text('2\naltered\nH 0 0 0\nH 0 0 0.75\n')
            value['result']['xyz'] = str(altered)
        else:
            value['result']['energy_protocol'] = copy.deepcopy(ACTIVE_SP)
            value['result']['energy_protocol']['solvent'] = 'water'
        record_path.write_text(json.dumps(value))
    else:
        value = json.loads(proof_path.read_text()); value['optimized_geometry_qualified'] = True
        proof_path.write_text(json.dumps(value))
    with pytest.raises(ValueError):
        module.point_evidence(pilot, review)
