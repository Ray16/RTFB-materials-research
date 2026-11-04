"""A symmetric viologen query must select the benzyl carbon, not N-methyl."""
import pytest
from rdkit import Chem

from redox.workflows.viologen_motif_fragments import benzyl_query

FRAGMENT='C[n+]1ccc(-c2cc[n+]([*:1])cc2)cc1'
GRAFT='Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1'


def test_constrained_attachment_query_excludes_the_other_n_methyl():
    parent=Chem.AddHs(Chem.MolFromSmiles(GRAFT))
    query=Chem.MolFromSmarts(benzyl_query(FRAGMENT))
    wildcard=next(a.GetIdx() for a in query.GetAtoms() if a.GetAtomicNum()==0)
    anchors={match[wildcard] for match in parent.GetSubstructMatches(query,uniquify=False)}
    assert len(anchors)==1
    atom=parent.GetAtomWithIdx(next(iter(anchors)))
    assert sum(a.GetAtomicNum()==1 for a in atom.GetNeighbors())==2
    assert any(a.GetIsAromatic() for a in atom.GetNeighbors())
    assert Chem.GetFormalCharge(parent)==2


@pytest.mark.parametrize('bad',['CC','[*:1]C[*:1]','[*]C'])
def test_ambiguous_or_missing_attachment_declarations_are_rejected(bad):
    with pytest.raises(ValueError):benzyl_query(bad)
