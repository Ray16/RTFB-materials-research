"""No GPU should be spent discovering diagnostic NumPy scalar serialization failures."""
import importlib.util,json
from pathlib import Path
import numpy as np
import pytest
from redox.core.common import ROOT


def module():
    from redox.core import metadata
    return metadata


def test_planned_allocation_shapes_and_metrics_roundtrip_as_plain_json():
    m=module();value=m.planned_allocation_metadata_check();roundtrip=json.loads(json.dumps(value,allow_nan=False))
    assert roundtrip['shape']==[632,632,714]
    assert type(roundtrip['bytes']) is int and type(roundtrip['relative_error']) is float
    assert roundtrip['elapsed_s']==.125 and roundtrip['native_gates_passed'] is True


def test_nested_numpy_diagnostic_metadata_is_checked_before_gpu_start():
    m=module();value=m.validate_metadata({'slots':(np.int64(0),np.int32(1)),'shape':np.array([3,3,7702,7702]),'path':Path('result.json')})
    assert json.loads(json.dumps(value))=={'slots':[0,1],'shape':[3,3,7702,7702],'path':'result.json'}
    with pytest.raises(ValueError):m.validate_metadata({'timing':np.float64(np.nan)})
    with pytest.raises(TypeError):m.validate_metadata({'metric':np.complex128(1+2j)})
