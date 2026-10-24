"""A memory forecast must cover both screened FP64 tensors without truncation."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('point_inventory', ROOT/'scripts/validation/inventory_polymer_point.py')
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_full_range_sum_matches_the_large_reference_inventory():
    shapes = [dict(shape=[7139,2857371],bytes=163190172552),
              dict(shape=[3018,2834461],bytes=68435226384)]
    assert MODULE.required_host_bytes(shapes) == 303720282499


@pytest.mark.parametrize('shapes', [[], [dict(shape=[2,3],bytes=48)],
    [dict(shape=[2,3],bytes=24)]*2, [dict(shape=[0,3],bytes=0)]*2,
    [dict(shape=[2.0,3],bytes=48)]*2])
def test_missing_arrays_or_inconsistent_storage_cannot_admit_a_job(shapes):
    with pytest.raises(ValueError):
        MODULE.required_host_bytes(shapes)
