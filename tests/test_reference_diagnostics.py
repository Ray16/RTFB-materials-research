"""Offset diagnostics must exclude the complete held group and preserve missing models."""
import pytest
from redox.validation.reference_diagnostics import cross_validated_offsets,organic_reference


def row(identifier,wave,error):
    return dict(id=identifier,event=str(wave),doi=identifier,wave=wave,error_V=error,calc_V=error,exp_V=0.)


def test_molecule_holdout_does_not_learn_from_its_other_redox_step():
    rows=[row('a',1,.1),row('a',2,.3),row('b',1,.8)]
    result=cross_validated_offsets(rows,'id')
    assert [p['fitted_offset_V'] for p in result['predictions']]==pytest.approx([.8,.8,.2])
    assert [p['training_count'] for p in result['predictions']]==[1,1,2]


def test_missing_training_step_is_not_filled_from_a_different_step():
    result=cross_validated_offsets([row('a',1,.1),row('b',2,.3)],'id',True)
    assert result['coverage']==0
    assert all(p['corrected_error_V'] is None for p in result['predictions'])


def test_organic_reference_excludes_anchor_molecule_and_cancels_common_shift():
    rows=[row('a',1,.2),row('a',2,.3),row('b',1,.5)]
    before=organic_reference(rows,'a')['Fc_independent_relative_metrics']
    for r in rows:r['calc_V']+=.4;r['error_V']+=.4
    after=organic_reference(rows,'a')['Fc_independent_relative_metrics']
    assert before['n']==1 and after['MAE_V']==pytest.approx(before['MAE_V'])
