import json

import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


@pytest.mark.parametrize('parent_status,expected', [
    (None, 'upstream_pending'),
    ({'status': 'GPU_geometry_stage_running'}, 'upstream_pending'),
    ({'status': 'geometry_stage_failed_preserved', 'finished': 1}, 'upstream_failed_preserved'),
])
def test_electronic_stage_waits_for_or_preserves_missing_failed_geometry(tmp_path, parent_status, expected):
    module = api.load_module('conformer_energy_dependency_test', ROOT / 'scripts/pipeline/prepare_conformer_electronic_batch.py')
    geometry = tmp_path / 'geometry'; geometry.mkdir()
    (geometry / 'plan.json').write_text('{}')
    if parent_status:
        (geometry / 'geometry').mkdir()
        (geometry / 'geometry/status.json').write_text(json.dumps(parent_status))
    out = tmp_path / 'electronic'
    entry = dict(geometry_campaign=str(geometry), campaign=str(out), geometry_plan_sha256=api.digest(geometry / 'plan.json'))
    assert module.prepare_ready(entry, {}) == expected
    assert not out.exists()
