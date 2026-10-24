"""Device headroom is an explicit resource criterion, separate from physics."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('native_switch_resource_policy',ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_native_resource_ceiling_preserves_one_GiB_above_its_measured_peak_limit():
    policy=m.device_limit(31744,32768)
    assert policy['minimum_device_headroom_MiB']==1024
    assert policy['workspace_estimate_GiB']==28
    assert policy['resource_and_scientific_criteria_separate']


@pytest.mark.parametrize('maximum,minimum',[(32768,32768),(0,32768),(31744,32000)])
def test_no_native_profile_with_missing_or_insufficient_device_headroom(maximum,minimum):
    with pytest.raises(ValueError):m.device_limit(maximum,minimum)
