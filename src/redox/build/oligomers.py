"""H-capped polystyrene fragments with explicit site occupancy and stereochemistry."""
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers,StereoEnumerationOptions
from redox.build.structures import unassigned_stereo


def oligomer_scaffold(n):
    if type(n) is not int or not 1<=n<=5:raise ValueError('Use one to five declared repeat units')
    aryl=lambda i:f'c1ccc(C[*:{i}])cc1'
    return ''.join('CC('+aryl(i)+')' for i in range(1,n))+'CC'+aryl(n)


def decorate_sites(fragment,n,occupied):
    sites=set(occupied)
    if not sites<=set(range(1,n+1)):raise ValueError('Site index outside oligomer')
    scaffold=Chem.RWMol(Chem.MolFromSmiles(oligomer_scaffold(n)))
    for atom in scaffold.GetAtoms():
        if atom.GetAtomicNum()==0 and atom.GetAtomMapNum() not in sites:
            atom.SetAtomicNum(17);atom.SetAtomMapNum(0);atom.SetNoImplicit(True)
    combined=scaffold.GetMol()
    for site in sorted(sites):
        frag=Chem.MolFromSmiles(fragment)
        if frag is None:raise ValueError('Invalid fragment')
        dummy=[a for a in frag.GetAtoms() if a.GetAtomicNum()==0]
        if len(dummy)!=1:raise ValueError('One motif attachment dummy required')
        dummy[0].SetAtomMapNum(site);combined=Chem.CombineMols(combined,frag)
    product=Chem.molzip(combined) if sites else combined
    Chem.SanitizeMol(product)
    if any(a.GetAtomicNum()==0 for a in product.GetAtoms()):raise ValueError('Unresolved attachment site')
    if len(Chem.GetMolFrags(product))!=1:raise ValueError('Oligomer must be one covalent fragment')
    options=StereoEnumerationOptions(onlyUnassigned=True,unique=True,maxIsomers=16)
    isomers=list(EnumerateStereoisomers(product,options=options))
    if not isomers:isomers=[product]
    if any(unassigned_stereo(mol) for mol in isomers):raise ValueError('Oligomer stereochemistry remains unspecified')
    return sorted({Chem.MolToSmiles(mol,isomericSmiles=True) for mol in isomers})
