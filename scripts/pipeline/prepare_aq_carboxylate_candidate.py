#!/usr/bin/env python
"""Construct the AQ-2-carboxylate capped graft from a source-identified acid."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from redox.build.grafting import grafting_reaction
from redox.build.structures import decorate
from redox.build.attachment_identity import exact_attachment


def build(acid_smiles):
    acid = Chem.MolFromSmiles(acid_smiles)
    if acid is None:
        raise ValueError('Invalid precursor structure')
    matches = acid.GetSubstructMatches(Chem.MolFromSmarts('[C](=[O])[O;H1]'))
    if len(matches) != 1:
        raise ValueError('Expected exactly one carboxylic acid handle')
    anchor = matches[0][2]
    graph = Chem.RWMol(acid)
    oxygen = graph.GetAtomWithIdx(anchor)
    oxygen.SetNumExplicitHs(0)
    oxygen.SetNoImplicit(True)
    dummy = Chem.Atom(0)
    dummy.SetAtomMapNum(1)
    graph.AddBond(anchor, graph.AddAtom(dummy), Chem.BondType.SINGLE)
    motif = graph.GetMol()
    Chem.SanitizeMol(motif)
    fragment = Chem.MolToSmiles(motif)
    product = decorate('Cc1ccc(C[*:1])cc1', fragment)
    reaction = grafting_reaction(fragment, Chem.MolToSmiles(product))
    if Chem.MolToInchiKey(Chem.MolFromSmiles(reaction['precursor_smiles'])) != Chem.MolToInchiKey(acid):
        raise ValueError('Precursor identity was not preserved')
    carbonyl = Chem.MolFromSmarts('[#6]=[O]')
    if len(acid.GetSubstructMatches(carbonyl)) != 3 or len(product.GetSubstructMatches(carbonyl)) != 3:
        raise ValueError('Expected both quinone carbonyls and the ester carbonyl')
    identity = dict(scaffold_smiles='Cc1ccc(C[*:1])cc1', attachment_fragment=fragment,
                    precursor_smiles=reaction['precursor_smiles'], product_smiles=reaction['product_smiles'])
    task = dict(identity=identity, payload=dict(attachment_fragment=fragment, family='anthraquinone',
                states=[['neutral', 0, 1, 0], ['reduced', -1, 2, -1]],
                precursor_smiles=reaction['precursor_smiles'], product_smiles=reaction['product_smiles']))
    mapping = exact_attachment(task)
    return dict(id='aq2_carboxylate_capped', name='AQ-2-carboxylate benzyl ester capped-site proposal',
                frag=fragment, precursor_inchikey=Chem.MolToInchiKey(acid),
                product_inchikey=Chem.MolToInchiKey(product),
                product_formula=rdMolDescriptors.CalcMolFormula(product),
                product_molecular_weight=Descriptors.MolWt(product),
                proposed_charge_spin_states=task['payload']['states'],
                state_assignment_scope='proposed parity-compatible states; physical spin and redox properties unverified',
                exact_design_key=mapping['exact_design_key'], route_binding=mapping['route_binding'],
                reaction=reaction)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--identity', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    identity = json.loads(args.identity.read_text())['PropertyTable']['Properties'][0]
    if identity['CID'] != 67030 or identity['InChIKey'] != 'ASDLSKCKYGVMAI-UHFFFAOYSA-N':
        raise ValueError('Unexpected source identity')
    result = build(identity['ConnectivitySMILES'])
    if result['precursor_inchikey'] != identity['InChIKey']:
        raise ValueError('Computed identity differs from source')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump([result], stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
