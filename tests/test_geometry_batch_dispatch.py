import json

import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def module():
    return api.load_module('geometry_batch_dispatch_test', ROOT / 'scripts/pipeline/coordinate_geometry_batch.py')


def test_host_promises_reserve_future_growth_and_free_floor():
    m = module(); gib = 2**30
    assert m.can_fit(80 * gib, 32 * gib)
    assert not m.can_fit(79 * gib, 32 * gib)
    assert not m.can_fit(100 * gib, 90 * gib)
    with pytest.raises(ValueError):
        m.can_fit(float('nan'), 0)


def test_undeclared_or_changed_campaign_cannot_start_a_gpu_process(tmp_path):
    m = module(); campaign = tmp_path / 'campaign'; campaign.mkdir()
    (campaign / 'plan.json').write_text('{}')
    (tmp_path / 'batch.json').write_text(json.dumps(dict(source_sha256={}, entries=[])))
    with pytest.raises(ValueError, match='exact declared batch'):
        m.launch_local(tmp_path, campaign, 0)
    (tmp_path / 'batch.json').write_text(json.dumps(dict(source_sha256={},
        entries=[dict(campaign=str(campaign), plan_sha256='wrong')])) )
    with pytest.raises(ValueError, match='exact declared batch'):
        m.launch_local(tmp_path, campaign, 0)
