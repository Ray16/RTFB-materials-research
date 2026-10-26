#!/usr/bin/env python3
"""Build an exact local N-linked 4-amino-TEMPO graft; no route claim."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from redox.build.attachment_identity import exact_attachment
from redox.build.grafting import grafting_reaction
from redox.build.structures import decorate


def build(identity):
    if identity['CID'] != 550942 or identity['InChIKey'] != 'XUXUHDYTLNCYQQ-UHFFFAOYSA-N':
        raise ValueError('Unexpected source precursor identity')
    precursor = Chem.MolFromSmiles(identity['ConnectivitySMILES'])
    if precursor is None or Chem.MolToInchiKey(precursor) != identity['InChIKey']:
        raise ValueError('Precursor graph does not match source identity')
    matches = precursor.GetSubstructMatches(Chem.MolFromSmarts('[N;H2;X3]-[C;R]'))
    if len(matches) != 1:
        raise ValueError('Expected one exocyclic primary amino handle')
    atom_index = matches[0][0]
    editable = Chem.RWMol(precursor)
    anchor = editable.GetAtomWithIdx(atom_index)
    anchor.SetNumExplicitHs(1)
    anchor.SetNoImplicit(True)
    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(1)
    editable.AddBond(atom_index, editable.AddAtom(dummy), Chem.BondType.SINGLE)
    fragment = editable.GetMol()
    Chem.SanitizeMol(fragment)
    frag = Chem.MolToSmiles(fragment)
    product = decorate('Cc1ccc(C[*:1])cc1', frag)
    product_smiles = Chem.MolToSmiles(product)
    reaction = grafting_reaction(frag, product_smiles)
    if Chem.MolToInchiKey(Chem.MolFromSmiles(reaction['precursor_smiles'])) != identity['InChIKey']:
        raise ValueError('Mapped precursor identity changed')
    def radical_count(mol):
        return sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
    if radical_count(precursor) != 1 or radical_count(product) != 1:
        raise ValueError('The precursor or product lost the nitroxide radical')
    states = [['rad', 0, 2, 0], ['ox', 1, 1, 1]]
    task = dict(identity=dict(scaffold_smiles='Cc1ccc(C[*:1])cc1', attachment_fragment=frag,
                              precursor_smiles=reaction['precursor_smiles'], product_smiles=product_smiles),
                payload=dict(attachment_fragment=frag, family='nitroxide', states=states,
                             precursor_smiles=reaction['precursor_smiles'], product_smiles=product_smiles))
    binding = exact_attachment(task)
    return dict(id='amino_tempo_capped', name='4-amino-TEMPO N-benzyl capped-site proposal',
                frag=frag, precursor_inchikey=identity['InChIKey'],
                product_inchikey=Chem.MolToInchiKey(product),
                product_formula=rdMolDescriptors.CalcMolFormula(product),
                product_molecular_weight=Descriptors.MolWt(product),
                proposed_charge_spin_states=states,
                state_assignment_scope='proposed states; physical spin and redox properties unverified',
                exact_design_key=binding['exact_design_key'], route_binding=binding['route_binding'],
                reaction=reaction)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--identity', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    identity = json.loads(args.identity.read_text())['PropertyTable']['Properties'][0]
    with args.out.open('x') as stream:
        json.dump([build(identity)], stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
