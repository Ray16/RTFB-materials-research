"""Smaller host admission needs an actual complete-pilot memory witness."""
import math
import pytest
from redox.workflows.counterion_thermochemistry import followup_host_bound


@pytest.mark.parametrize('peak,bound',[(2,16),(12.8,16),(13,40),(39,40)])
def test_measured_peak_and_headroom_choose_a_declared_host_envelope(peak,bound):
    assert followup_host_bound(dict(peak_host_RSS_bytes=peak*2**30))==bound


@pytest.mark.parametrize('peak',[None,True,0,-1,math.nan,math.inf,41*2**30])
def test_missing_or_out_of_class_memory_cannot_release_a_smaller_profile(peak):
    with pytest.raises(ValueError):followup_host_bound(dict(peak_host_RSS_bytes=peak))
