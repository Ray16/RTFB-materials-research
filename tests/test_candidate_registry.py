import pytest
from redox.screening.candidate_registry import merge_sources


def row(identifier='one', frag='[*:1]OC'):
    return dict(id=identifier, frag=frag, family='test', states=[['neu',0,1,0],['ox',1,2,1]])


def test_merge_preserves_aliases_and_sources_without_inheriting_properties():
    a=row(); b=row('two'); b['E1_vs_Fc_V']=99
    result=merge_sources([('enum','enumerated',[a]),('eval','evaluated_candidate',[b])])
    assert len(result)==1 and result[0]['aliases']==['one','two']
    assert len(result[0]['declarations'])==2 and result[0]['evaluations'][0]['source']=='eval'
    assert result[0]['evaluations'][0]['reported_values']['E1_vs_Fc_V']==99
    assert 'E1_vs_Fc_V' not in result[0] and not result[0]['evaluations'][0]['raw_records_revalidated']
    other=merge_sources([('enum','enumerated',[a,row('other','[*:1]OCC')])])
    assert len(other)==2 and all(not r['evaluations'] for r in other)


def test_conflicting_id_or_graph_and_parent_properties_rejected():
    with pytest.raises(ValueError,match='different chemical graphs'):
        merge_sources([('x','proposal',[row(),row(frag='[*:1]OCC')])])
    bad=row(); bad['product_smiles']='CCO'
    with pytest.raises(ValueError,match='differs'):
        merge_sources([('x','proposal',[bad])])
    bad=row(); bad['parent_metrics_inherited']=True
    with pytest.raises(ValueError,match='Parent properties'):
        merge_sources([('x','evaluated_candidate',[bad])])
