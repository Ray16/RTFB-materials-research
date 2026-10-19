from datetime import datetime
from zoneinfo import ZoneInfo

import pytest

from redox.workflows.publication import watch_deadline


def test_daily_window_uses_the_current_day_instead_of_retired_campaign_date():
    chicago=ZoneInfo('America/Chicago');now=datetime(2026,10,19,1,30,tzinfo=chicago)
    end=datetime.fromtimestamp(watch_deadline(now=now),chicago)
    assert end==datetime(2026,10,19,23,55,tzinfo=chicago)
    assert end>now


def test_explicit_window_counts_elapsed_hours_across_daylight_saving_change():
    now=datetime(2026,11,1,0,30,tzinfo=ZoneInfo('America/Chicago'))
    assert watch_deadline(3,now=now)-now.timestamp()==3*3600


@pytest.mark.parametrize('hours',[0,-1,float('nan')])
def test_invalid_watch_window_is_rejected(hours):
    with pytest.raises(ValueError):watch_deadline(hours)
