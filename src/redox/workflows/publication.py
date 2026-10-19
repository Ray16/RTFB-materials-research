"""Daily publication windows shared by result collectors and slide watchers."""
from datetime import datetime
import math
from zoneinfo import ZoneInfo


def watch_deadline(max_hours=None,*,now=None):
    """Stop at today's Chicago 23:55, or after explicit elapsed hours across midnight."""
    start=datetime.now(ZoneInfo('America/Chicago')) if now is None else now
    if start.tzinfo is None:raise ValueError('Timezone-aware watch start required')
    if max_hours is not None:
        if not math.isfinite(max_hours) or max_hours<=0:raise ValueError('max-hours must be positive and finite')
        return start.timestamp()+max_hours*3600
    local=start.astimezone(ZoneInfo('America/Chicago'))
    return local.replace(hour=23,minute=55,second=0,microsecond=0).timestamp()
