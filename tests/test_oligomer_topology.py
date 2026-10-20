"""Prevent cyclic/branched or incorrectly grafted structures masquerading as chains."""
import pytest
from redox.build.oligomers import decorate_sites
from redox.build.oligomer_topology import chain_topology


@pytest.mark.parametrize('n',[1,2,3])
def test_open_chain_has_two_n_backbone_carbons_and_para_sites(n):
    for smiles in decorate_sites('[*:1]Oc1ccccc1',n,range(1,n+1)):
        result=chain_topology(smiles,n)
        assert result['backbone_carbon_count']==2*n
        assert result['backbone_CC_bonds']==2*n-1
        assert result['chain_ends']==2 and not result['backbone_ring_closure']
        assert [s['backbone_position'] for s in result['para_graft_sites']]==list(range(2,2*n+1,2))


def test_meta_graft_is_not_accepted_as_the_para_model():
    with pytest.raises(ValueError,match='para'):chain_topology('CCc1cccc(COc2ccccc2)c1',1)


def test_a_cyclic_backbone_and_wrong_repeat_count_are_rejected():
    cyclic='C1C(c2ccc(COc3ccccc3)cc2)CC(c2ccc(COc3ccccc3)cc2)CC1c2ccc(COc3ccccc3)cc2'
    with pytest.raises(ValueError):chain_topology(cyclic,3)
    with pytest.raises(ValueError):chain_topology(decorate_sites('[*:1]Oc1ccccc1',2,(1,2))[0],3)
