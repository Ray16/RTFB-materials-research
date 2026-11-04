import pytest
from redox.properties.contact_identity import publication_identity,validate_identities


def test_ids_are_stable_and_independent_of_roster_order():
    rows=[dict(original_start='N0_outer',branch='plus'),dict(original_start='N0_above',branch='minus')]
    assert [r[0] for r in validate_identities(rows,['N0_above','N0_outer'])]==['N0_outer_plus','N0_minus']
    assert publication_identity('separated_start',None)==('separated_start','separated start')


@pytest.mark.parametrize('rows,starts',[
    ([dict(original_start='N0_outer',branch='other')],['N0_outer']),
    ([dict(original_start='N1_above',branch=None)],['N0_outer']),
    ([dict(original_start='N0_outer',branch='plus')]*2,['N0_outer']),
    ([dict(original_start='N0_outer',branch=None)],['N0_outer','N0_outer']),
])
def test_undeclared_duplicate_or_invalid_identities_fail(rows,starts):
    with pytest.raises(ValueError):validate_identities(rows,starts)
