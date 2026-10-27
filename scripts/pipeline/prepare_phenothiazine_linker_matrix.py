#!/usr/bin/env python3
"""Build exact, route-reviewable phenothiazine linker/substitution hypotheses."""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem,rdBase
from rdkit.Chem import Descriptors,rdMolDescriptors

from redox.build.attachment_identity import exact_attachment
from redox.build.grafting import grafting_reaction
from redox.build.structures import decorate
from redox.core.common import load_config


def build(spec,baseline,ledger):
    if (spec.get('schema')!=1 or spec.get('core')!='phenothiazine'
            or set(spec['anchors'])!={'flexible','cyclic'} or set(spec['substitution_cores'])!={'H','3_7_dimethoxy'}):
        raise ValueError('Unexpected declared phenothiazine design factors')
    if len(baseline)!=1 or baseline[0]['precursor_inchikey']!='FMENRHNHBUWRMN-UHFFFAOYSA-N':
        raise ValueError('Source-identified hydroxyethyl baseline required')
    known={}
    for entry in ledger:
        key=entry['product_key']
        if key in known:raise ValueError('Duplicate product identity in frozen ledger')
        known[key]=entry
    scaffold=load_config('family_search').SCAFFOLD;rows=[]
    for linker,anchor in spec['anchors'].items():
        for substituent,core in spec['substitution_cores'].items():
            motif=Chem.MolFromSmiles(anchor+core)
            if motif is None:raise ValueError('Invalid declared motif')
            fragment=Chem.MolToSmiles(motif,isomericSmiles=True)
            product=decorate(scaffold,fragment);smiles=Chem.MolToSmiles(product,isomericSmiles=True)
            reaction=grafting_reaction(fragment,smiles);precursor=Chem.MolFromSmiles(reaction['precursor_smiles']);key=Chem.MolToInchiKey(product)
            identity=dict(scaffold_smiles=scaffold,attachment_fragment=fragment,precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles)
            task=dict(identity=identity,payload=dict(attachment_fragment=fragment,family='phenothiazine p-type',states=spec['states'],
                precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles))
            binding=exact_attachment(task);centers=Chem.FindMolChiralCenters(precursor,includeUnassigned=True)
            if linker=='cyclic' and (len(centers)!=2 or any(c[1]=='?' for c in centers)):
                raise ValueError('The cyclic linker needs explicit stereochemistry')
            if linker=='flexible' and centers:raise ValueError('Unexpected flexible-precursor stereochemistry')
            prior=known.get(key)
            if prior:
                if Chem.MolToSmiles(Chem.MolFromSmiles(prior['grafted_smiles']),isomericSmiles=True)!=smiles:
                    raise ValueError('Frozen ledger identity conflicts with product graph')
                identifier=prior['id'];origin='reused_frozen_ledger'
            elif key==baseline[0]['product_inchikey']:
                if binding['exact_design_key']!=baseline[0]['exact_design_key']:
                    raise ValueError('Baseline graft route differs from prior source assessment')
                identifier=baseline[0]['id'];origin='reused_source_proposal'
            else:
                identifier='hyp_ptz_linker_'+hashlib.sha256(smiles.encode()).hexdigest()[:12];origin='new_exact_graph'
            rows.append(dict(id=identifier,name=f'phenothiazine / {linker} / {substituent}',frag=fragment,
                family='phenothiazine p-type',states=spec['states'],design_factors=dict(linker=linker,substituent=substituent),
                precursor_inchikey=Chem.MolToInchiKey(precursor),product_inchikey=key,product_smiles=smiles,
                product_formula=rdMolDescriptors.CalcMolFormula(product),product_molecular_weight=Descriptors.MolWt(product),
                product_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(product),precursor_stereocenters=centers,
                exact_design_key=binding['exact_design_key'],route_binding=binding['route_binding'],reaction=reaction,
                origin=origin,properties_computed=False,synthesis_eligibility_established=False,
                state_scope='One oxidation is a proposal; no reversible behavior or usable capacity inferred'))
    if len(rows)!=4 or len({r['product_inchikey'] for r in rows})!=4:raise ValueError('Four distinct products required')
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('spec','baseline','ledger','out','candidates-out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.out.exists() or a.candidates_out.exists():raise FileExistsError('Preserve existing candidate outputs')
    rows=build(json.loads(a.spec.read_text()),json.loads(a.baseline.read_text()),json.loads(a.ledger.read_text())['candidates'])
    with a.out.open('x') as stream:json.dump(dict(status='exact_designs_prepared_route_review_required',rdkit_version=rdBase.rdkitVersion,candidates=rows),stream,indent=2);stream.write('\n')
    with a.candidates_out.open('x') as stream:json.dump(rows,stream,indent=2);stream.write('\n')
