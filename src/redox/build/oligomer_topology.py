"""Independently verify an open, H-capped styrene backbone and para graft sites."""
from rdkit import Chem


def chain_topology(smiles, repeats):
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol))!=1:raise ValueError('One intact covalent chain required')
    nodes={a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==6 and not a.IsInRing()}
    remaining=set(nodes);components=[]
    while remaining:
        todo=[min(remaining)];part=set()
        while todo:
            i=todo.pop()
            if i in part:continue
            part.add(i);todo.extend(a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() in nodes and a.GetIdx() not in part)
        remaining-=part;components.append(part)
    candidates=[part for part in components if len(part)==2*repeats and
                sum(any(a.GetIsAromatic() and a.GetAtomicNum()==6 for a in mol.GetAtomWithIdx(i).GetNeighbors()) for i in part)==repeats]
    if len(candidates)!=1:raise ValueError('Cannot identify a unique 2n-carbon styrene backbone')
    backbone=candidates[0]
    adjacent={i:[a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() in backbone] for i in backbone}
    ends=[i for i in backbone if len(adjacent[i])==1]
    edges=sum(map(len,adjacent.values()))//2
    if edges!=2*repeats-1 or len(ends)!=2 or max(map(len,adjacent.values()))>2:
        raise ValueError('Backbone must be an open unbranched carbon path')
    methyl=[i for i in ends if mol.GetAtomWithIdx(i).GetTotalNumHs()==3]
    if len(methyl)!=1:raise ValueError('Expected one terminal methyl H-cap')
    order=[methyl[0]]
    while len(order)<len(backbone):order.append(next(i for i in adjacent[order[-1]] if i not in order))
    sites=[];rings=mol.GetRingInfo().AtomRings()
    for position,index in enumerate(order,1):
        aromatic=[a.GetIdx() for a in mol.GetAtomWithIdx(index).GetNeighbors() if a.GetIsAromatic()]
        if not aromatic:continue
        if len(aromatic)!=1:raise ValueError('Backbone carbon has ambiguous aryl substitution')
        ring=next((set(r) for r in rings if aromatic[0] in r and len(r)==6 and all(mol.GetAtomWithIdx(i).GetAtomicNum()==6 for i in r)),None)
        if ring is None:raise ValueError('Styrene side group must contain a six-carbon phenyl ring')
        linkers=[(i,a.GetIdx()) for i in ring for a in mol.GetAtomWithIdx(i).GetNeighbors()
                 if a.GetAtomicNum()==6 and a.GetIdx() not in ring and a.GetIdx() not in backbone and not a.IsInRing()]
        if len(linkers)!=1:raise ValueError('Expected one benzyl graft linker on each phenyl ring')
        aryl,linker=linkers[0]
        if len(Chem.GetShortestPath(mol,aromatic[0],aryl))-1!=3:raise ValueError('Graft linker is not para to the backbone')
        sites.append(dict(backbone_position=position,backbone_atom=index,phenyl_atoms=sorted(ring),benzyl_atom=linker))
    if [s['backbone_position'] for s in sites]!=list(range(2,2*repeats+1,2)):
        raise ValueError('Aryl groups must alternate along the 2n-carbon backbone')
    return dict(topology='open_linear_chain',repeat_units=repeats,backbone_atoms=order,
                backbone_carbon_count=len(order),backbone_CC_bonds=edges,chain_ends=2,
                backbone_ring_closure=False,para_graft_sites=sites)
