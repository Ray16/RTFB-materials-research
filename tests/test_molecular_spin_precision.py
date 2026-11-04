import pytest

from redox.validation.molecular_native import spin_record_agreement


def test_four_decimal_native_record_matches_the_unrounded_density():
    result=spin_record_agreement(dict(s_squared=.7535607357117016),dict(s_squared=.7536,spin_contam=.0036),2)
    assert result['stored_decimal_places']==4
    assert result['full_precision_contamination']==pytest.approx(.0035607357117016)
    assert result['full_precision_limit']==.1


@pytest.mark.parametrize('value,saved',[(.7535607,.754),(.7535607,float('nan')),(float('nan'),.7536)])
def test_mismatched_or_nonfinite_spin_records_remain_rejected(value,saved):
    with pytest.raises(ValueError):spin_record_agreement(dict(s_squared=value),dict(s_squared=saved,spin_contam=.0036),2)


def test_rounding_cannot_hide_a_real_physical_spin_violation():
    with pytest.raises(ValueError,match='Full-precision'):
        spin_record_agreement(dict(s_squared=.850001),dict(s_squared=.85,spin_contam=.1),2)
