#!/usr/bin/env python3
"""Construct a proposed O-linked graft of source-identified hydroxyethyl PTZ."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdMolDescriptors

from redox.build.attachment_identity import exact_attachment
from redox.build.grafting import grafting_reaction
from redox.build.structures import decorate
from redox.core.common import load_config


def build(identity):
    if identity['CID']!=352306 or identity['InChIKey']!='FMENRHNHBUWRMN-UHFFFAOYSA-N':
        raise ValueError('Unexpected source precursor identity')
    precursor=Chem.MolFromSmiles(identity['ConnectivitySMILES'])
    if precursor is None or Chem.MolToInchiKey(precursor)!=identity['InChIKey']:
        raise ValueError('Precursor graph does not match its source identity')
    if rdMolDescriptors.CalcMolFormula(precursor)!='C14H13NOS':
        raise ValueError('Unexpected precursor formula')
    alcohol=precursor.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]'))
    if len(alcohol)!=1:raise ValueError('Expected one alcohol attachment handle')
    anchor=alcohol[0][0];edit=Chem.RWMol(precursor);oxygen=edit.GetAtomWithIdx(anchor)
    oxygen.SetNumExplicitHs(0);oxygen.SetNoImplicit(True)
    dummy=Chem.Atom(0);dummy.SetAtomMapNum(1);edit.AddBond(anchor,edit.AddAtom(dummy),Chem.BondType.SINGLE)
    motif=edit.GetMol();Chem.SanitizeMol(motif);fragment=Chem.MolToSmiles(motif,isomericSmiles=True)
    scaffold=load_config('family_search').SCAFFOLD
    product=decorate(scaffold,fragment);smiles=Chem.MolToSmiles(product,isomericSmiles=True)
    reaction=grafting_reaction(fragment,smiles)
    if Chem.MolToInchiKey(Chem.MolFromSmiles(reaction['precursor_smiles']))!=identity['InChIKey']:
        raise ValueError('Attachment mapping altered the precursor')
    states=[['neutral',0,1,0],['oxidized',1,2,1]]
    binding_identity=dict(scaffold_smiles=scaffold,attachment_fragment=fragment,
                          precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles)
    task=dict(identity=binding_identity,payload=dict(attachment_fragment=fragment,family='phenothiazine p-type',
        states=states,precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles))
    binding=exact_attachment(task)
    return dict(id='source_ptz_hydroxyethyl_'+hashlib.sha256(smiles.encode()).hexdigest()[:12],
        name='hydroxyethyl phenothiazine O-linked capped-site proposal',frag=fragment,
        family='phenothiazine p-type',states=states,precursor_inchikey=identity['InChIKey'],
        precursor_formula=rdMolDescriptors.CalcMolFormula(precursor),product_inchikey=Chem.MolToInchiKey(product),
        product_smiles=smiles,product_formula=rdMolDescriptors.CalcMolFormula(product),
        product_molecular_weight=Descriptors.MolWt(product),exact_design_key=binding['exact_design_key'],
        route_binding=binding['route_binding'],reaction=reaction,rdkit_version=rdBase.rdkitVersion,
        state_scope='one oxidation proposed; charge/spin and reversible operation unverified',
        properties_computed=False,synthesis_eligibility_established=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--identity',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    row=build(json.loads(a.identity.read_text())['PropertyTable']['Properties'][0])
    with a.out.open('x') as stream:json.dump([row],stream,indent=2);stream.write('\n')
