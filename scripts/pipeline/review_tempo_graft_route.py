#!/usr/bin/env python
"""Bind an existing TEMPO local graft to an exact-support source precedent."""
import argparse
from collections import Counter
import csv
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path

from rdkit import Chem

from redox.build.attachment_identity import exact_attachment
from redox.build.grafting import atom_inventory, grafting_reaction
from redox.build.structures import decorate
from redox.core.common import load_config


class Text(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts = []

    def handle_data(self, data):
        self.parts.append(data)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def balanced_smiles(left, right):
    sides = [[Chem.MolFromSmiles(value) for value in side] for side in (left,right)]
    if any(m is None for side in sides for m in side):
        raise ValueError('Invalid species in declared reaction')
    atoms = [sum((atom_inventory(m) for m in side),Counter()) for side in sides]
    charges = [sum(Chem.GetFormalCharge(m) for m in side) for side in sides]
    radicals = [sum(a.GetNumRadicalElectrons() for m in side for a in m.GetAtoms()) for side in sides]
    if atoms[0] != atoms[1] or charges[0] != charges[1] or radicals[0] != radicals[1]:
        raise ValueError('Reaction does not preserve atom inventory, charge and radical count')
    return dict(atom_balanced=True, charge_balanced=True, radical_electrons_balanced=True,
                charge=charges[0], radical_electrons=radicals[0])


def build(patent, identity, manifest):
    parser = Text()
    parser.feed(Path(patent).read_text())
    text = ' '.join(' '.join(parser.parts).split())
    start = text.find('Example A4 Preparation of Compound 12')
    end = text.find('Examples B:',start)
    section = text[start:end] if start >= 0 and end > start else ''
    required = ('4-hydroxy-TEMPO','NaH','Merrifield polymer','1% divinylbenzene',
                '1.7 mmol of Cl/g','0.9 mmol of N-oxyl/g','DMF or THF','1 to 4 days')
    if not all(term in section for term in required) or '4-oxy-bound to a Merrifield polymer' not in text:
        raise ValueError('Patent lacks the expected exact-support route evidence')
    pubchem = json.loads(Path(identity).read_text())['PropertyTable']['Properties'][0]
    if pubchem['CID'] != 137994 or pubchem['InChIKey'] != 'UZFMOKQJFYMBGY-UHFFFAOYSA-N':
        raise ValueError('Unexpected 4-hydroxy-TEMPO identity source')
    config = load_config('redox_groups')
    candidates = [v for v in vars(config).values() if isinstance(v,list)
                  for g in v if isinstance(g,dict) and g.get('id') == 'tempo']
    if len(candidates) != 1:
        raise ValueError('Expected one existing TEMPO candidate definition')
    group = next(g for g in candidates[0] if isinstance(g,dict) and g.get('id') == 'tempo')
    scaffold = 'Cc1ccc(C[*:1])cc1'
    product = Chem.MolToSmiles(decorate(scaffold,group['frag']))
    reaction = grafting_reaction(group['frag'],product)
    precursor = Chem.MolFromSmiles(reaction['precursor_smiles'])
    if (Chem.MolToInchiKey(precursor) != pubchem['InChIKey']
            or sum(a.GetNumRadicalElectrons() for a in precursor.GetAtoms()) != 1):
        raise ValueError('Precursor differs from the source nitroxide')
    with Path(manifest).open(newline='') as stream:
        existing = [row for row in csv.DictReader(stream) if row['id'] == 'tempo']
    if len(existing) != 3 or any(Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles'])) != product for row in existing):
        raise ValueError('Existing state manifest does not match the exact local graft')
    identity_graph = dict(scaffold_smiles=scaffold,attachment_fragment=group['frag'],
                          precursor_smiles=reaction['precursor_smiles'],product_smiles=product)
    payload = dict(attachment_fragment=group['frag'],family=group['family'],states=group['states'],
                   precursor_smiles=reaction['precursor_smiles'],product_smiles=product)
    mapped = exact_attachment(dict(identity=identity_graph,payload=payload))
    left = [reaction['precursor_smiles'],reaction['site_smiles'],'[Na+].[H-]']
    right = [product,'[H][H]','[Na+].[Cl-]']
    balance = balanced_smiles(left,right)
    evidence = dict(route_binding=mapped['route_binding'],exact_design_key=mapped['exact_design_key'],
        exact_precursor_access_reviewed=True,attachment_selectivity_reviewed=True,conditions_defined=True,
        balanced_route_reviewed=True,redox_core_retention_reviewed=True,
        conditions=dict(solvent='DMF or THF (source)',reactive_form='4-hydroxy-TEMPO with NaH',
                        base='NaH',counterions='Na+ and Cl- in overall stoichiometric bookkeeping',
                        temperature='0 C addition, then room temperature',time='1 to 4 days after resin addition'),
        source_urls=['https://patents.google.com/patent/US6441243B1/en','https://pubchem.ncbi.nlm.nih.gov/compound/137994'],
        evidence_scope='Patent Example A4 reports 4-oxy-bound TEMPO on 1% DVB Merrifield resin and 0.9 mmol N-oxyl/g after washing; the capped product is its local graph model.',
        limitations='One patent report; N-oxyl/g and starting chloride/g have different mass bases and do not define conversion. Selectivity and radical retention are source claims, not independent replication. Electrochemical cycling, loading in our beads, permeability, grafting free energy, motif interaction and device performance remain unresolved.',
        balance_scope='Overall NaH stoichiometry is atom/charge/radical balanced; ion pairing, phases, kinetics and reaction free energy are not modeled.')
    candidate = dict(id='tempo',name=group['name'],family=group['family'],frag=group['frag'],
                     product_smiles=product,precursor_smiles=reaction['precursor_smiles'],
                     source_scope='existing local capped model, not new molecule',states=group['states'])
    return dict(candidate=candidate,evidence=evidence,neutral_HCl_reference=reaction,
                NaH_overall=dict(reactants=left,products=right,**balance),
                source_sha256={str(p):sha(p) for p in (patent,identity,manifest,Path(__file__).resolve())},
                full_material_eligibility=False,experimental_electrochemical_cycling_inferred=False)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('patent','identity','manifest','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();result=build(a.patent.resolve(),a.identity.resolve(),a.manifest.resolve())
    a.out.mkdir(parents=True,exist_ok=False)
    for name,data in [('candidate.json',[result['candidate']]),('evidence.json',{'candidates':{'tempo':result['evidence']}}),('route_review.json',result)]:
        (a.out/name).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'status':'source_route_review_prepared','candidate':'tempo'}))


if __name__=='__main__':main()
