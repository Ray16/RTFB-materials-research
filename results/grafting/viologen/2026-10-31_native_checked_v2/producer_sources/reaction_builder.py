"""Balanced local Merrifield grafting reactions; no assumed base or kinetics."""
from collections import Counter
from rdkit import Chem

SITE_SMILES='Cc1ccc(CCl)cc1'


def atom_inventory(molecule):
    return Counter(a.GetAtomicNum() for a in Chem.AddHs(molecule).GetAtoms())


def precursor(fragment):
    """Undo attachment: an N-H/O-H precursor, or the free quaternizing N."""
    mol=Chem.MolFromSmiles(fragment)
    if mol is None:raise ValueError('Invalid motif fragment')
    dummies=[a for a in mol.GetAtoms() if a.GetAtomicNum()==0]
    if len(dummies)!=1 or dummies[0].GetDegree()!=1:raise ValueError('One terminal attachment dummy required')
    dummy=dummies[0];anchor=dummy.GetNeighbors()[0]
    if anchor.GetSymbol() not in ('N','O'):raise ValueError('Unreviewed reaction class')
    quaternization=anchor.GetSymbol()=='N' and anchor.GetFormalCharge()==1
    index=anchor.GetIdx();hydrogens=anchor.GetTotalNumHs();charge=anchor.GetFormalCharge()
    rw=Chem.RWMol(mol);target=rw.GetAtomWithIdx(index)
    target.SetFormalCharge(charge-1 if quaternization else charge)
    target.SetNumExplicitHs(hydrogens if quaternization else hydrogens+1);target.SetNoImplicit(True)
    rw.RemoveAtom(dummy.GetIdx());result=rw.GetMol();Chem.SanitizeMol(result)
    return result,('N_quaternization' if quaternization else 'NH_OH_substitution')


def grafting_reaction(fragment,product_smiles):
    reactant,kind=precursor(fragment);product=Chem.MolFromSmiles(product_smiles)
    if product is None:raise ValueError('Invalid grafted product')
    from redox.build.structures import decorate
    expected=decorate('Cc1ccc(C[*:1])cc1',fragment)
    if Chem.MolToSmiles(expected)!=Chem.MolToSmiles(product):
        raise ValueError('Product does not preserve the declared motif and attachment site')
    site=Chem.MolFromSmiles(SITE_SMILES)
    leaving=Chem.MolFromSmiles('[Cl-]' if kind=='N_quaternization' else '[H]Cl')
    if atom_inventory(reactant)+atom_inventory(site)!=atom_inventory(product)+atom_inventory(leaving):
        raise ValueError('Grafting reaction is not atom balanced')
    q=lambda mol:Chem.GetFormalCharge(mol)
    if q(reactant)+q(site)!=q(product)+q(leaving):raise ValueError('Grafting reaction is not charge balanced')
    return dict(reaction_class=kind,precursor_smiles=Chem.MolToSmiles(reactant),site_smiles=SITE_SMILES,product_smiles=Chem.MolToSmiles(product),byproduct_smiles=Chem.MolToSmiles(leaving),precursor_charge=q(reactant),product_charge=q(product),byproduct_charge=q(leaving),base_correction_required=kind!='N_quaternization',atom_balanced=True,charge_balanced=True,model='capped-site substitution; explicit byproduct; no base correction')
