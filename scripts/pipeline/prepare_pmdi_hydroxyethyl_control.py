#!/usr/bin/env python3
"""Construct the local O-graft of source-reported mono(hydroxyethyl) PMDI.

The precursor retains a second, N-H imide handle. No attachment selectivity or
material property is inferred from the source precursor report.
"""
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


FRAGMENT = '[*:1]OCCn1c(=O)c2cc3c(=O)[nH]c(=O)c3cc2c1=O'
PRECURSOR_KEY = 'HOUHMLHIFBPDGF-UHFFFAOYSA-N'


def build():
    config = load_config('family_search')
    fragment = Chem.MolToSmiles(Chem.MolFromSmiles(FRAGMENT), isomericSmiles=True)
    product = decorate(config.SCAFFOLD,fragment)
    smiles = Chem.MolToSmiles(product,isomericSmiles=True)
    reaction = grafting_reaction(fragment,smiles)
    precursor = Chem.MolFromSmiles(reaction['precursor_smiles'])
    if (Chem.MolToInchiKey(precursor)!=PRECURSOR_KEY
            or rdMolDescriptors.CalcMolFormula(precursor)!='C12H8N2O5'):
        raise ValueError('Precursor graph differs from the named mono(hydroxyethyl) PMDI')
    if len(precursor.GetSubstructMatches(Chem.MolFromSmarts('[nH]'))) != 1:
        raise ValueError('Expected a competing imide N-H handle')
    if len(precursor.GetSubstructMatches(Chem.MolFromSmarts('[O;H1]'))) != 1:
        raise ValueError('Expected one alcohol handle')
    if len(product.GetSubstructMatches(Chem.MolFromSmarts('[nH]'))) != 1:
        raise ValueError('O-linked local product must retain the imide N-H')
    spec=config.CORES['pmdi'];identity=dict(scaffold_smiles=config.SCAFFOLD,attachment_fragment=fragment,
        precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles)
    task=dict(identity=identity,payload=dict(attachment_fragment=fragment,family=spec['family'],
        states=spec['states'],precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles))
    binding=exact_attachment(task)
    return dict(id='source_pmdi_hydroxyethyl_nh_'+hashlib.sha256(smiles.encode()).hexdigest()[:12],
        name='mono(hydroxyethyl) PMDI O-linked capped-site proposal', frag=fragment,
        family=spec['family'], states=spec['states'], precursor_inchikey=PRECURSOR_KEY,
        precursor_formula=rdMolDescriptors.CalcMolFormula(precursor),
        product_inchikey=Chem.MolToInchiKey(product),product_smiles=smiles,
        product_formula=rdMolDescriptors.CalcMolFormula(product),
        product_molecular_weight=Descriptors.MolWt(product),
        exact_design_key=binding['exact_design_key'],route_binding=binding['route_binding'],reaction=reaction,
        competing_handle='precursor contains both primary alcohol and imide N-H',
        precursor_source='https://patents.google.com/patent/US4320129A/en',
        source_scope='Preparation 4 reports isolation of the named precursor; source claim not independent replication',
        rdkit_version=rdBase.rdkitVersion,properties_computed=False,synthesis_eligibility_established=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with a.out.open('x') as stream:
        json.dump([build()],stream,indent=2);stream.write('\n')
