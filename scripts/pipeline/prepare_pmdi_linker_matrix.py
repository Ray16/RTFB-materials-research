#!/usr/bin/env python3
"""Prepare a route-gated 2x2 PMDI linker/edge-substituent hypothesis panel.

The methyl group is a proposed steric perturbation, not proof of shielding.
The cyclic linker has two explicit stereocenters; precursor access and grafting
conditions remain independent questions.
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
from redox.screening.family_search import ring_variants


def build(library):
    config = load_config('family_search')
    spec = config.CORES['pmdi']
    known = {row['inchikey']: row for row in library}
    if len(known) != len(library):
        raise ValueError('Duplicate InChIKey in source library')
    anchors = {'flexible': config.ANCHORS['ether_C2'],
               'cyclic': '[*:1]O[C@H]1C[C@H]1'}
    result = []
    for linker, anchor in anchors.items():
        template = spec['template'].format(anchor=anchor, cap='CC')
        variants = {}
        for label, site, frag in ring_variants(template, {'H':None,'Me':'C'}):
            variants.setdefault((label, frag), []).append(site)
        for edge in ('H','Me'):
            selected = [(frag, sites) for (label,frag),sites in variants.items() if label == edge]
            if len(selected) != 1:
                raise ValueError('PMDI edge variant is ambiguous under the chosen template')
            frag, sites = selected[0]
            product = decorate(config.SCAFFOLD, frag)
            smiles = Chem.MolToSmiles(product, isomericSmiles=True)
            reaction = grafting_reaction(frag, smiles)
            precursor = Chem.MolFromSmiles(reaction['precursor_smiles'])
            inchikey = Chem.MolToInchiKey(product)
            previous = known.get(inchikey)
            if previous is not None and previous['grafted_smiles'] != smiles:
                raise ValueError('Library graph differs despite matching identity')
            candidate_id = previous['group']['id'] if previous else 'search_pmdi_'+hashlib.sha256(smiles.encode()).hexdigest()[:12]
            identity = dict(scaffold_smiles=config.SCAFFOLD, attachment_fragment=frag,
                            precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles)
            task = dict(identity=identity,payload=dict(attachment_fragment=frag,family=spec['family'],
                        states=spec['states'],precursor_smiles=reaction['precursor_smiles'],product_smiles=smiles))
            binding = exact_attachment(task)
            centers = Chem.FindMolChiralCenters(precursor, includeUnassigned=True)
            if linker == 'cyclic' and (len(centers) != 2 or any(c[1]=='?' for c in centers)):
                raise ValueError('Cyclic linker stereochemistry is not explicit')
            if linker == 'flexible' and centers:
                raise ValueError('Unexpected stereochemistry in flexible precursor')
            result.append(dict(id=candidate_id, name=f'PMDI ethyl / {linker} / {edge}',
                frag=frag, family=spec['family'], states=spec['states'],
                design_factors=dict(linker=linker, edge_substituent=edge,
                    interpretation='methyl is a hypothesized steric perturbation, not demonstrated shielding'),
                template_atom_indices=sites, precursor_inchikey=Chem.MolToInchiKey(precursor),
                product_inchikey=inchikey, product_smiles=smiles,
                product_formula=rdMolDescriptors.CalcMolFormula(product),
                product_molecular_weight=Descriptors.MolWt(product),
                product_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(product),
                precursor_stereocenters=centers, exact_design_key=binding['exact_design_key'],
                route_binding=binding['route_binding'], reaction=reaction,
                source_library_id=previous['group']['id'] if previous else None,
                origin='reused_existing_graph' if previous else 'new_exact_graph',
                properties_computed=False, synthesis_eligibility_established=False))
    if len(result)!=4 or len({row['product_inchikey'] for row in result})!=4:
        raise ValueError('Panel must have four distinct exact graphs')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--candidates-out',type=Path,
                   help='Optional route-gate input containing just the candidate list')
    args=p.parse_args()
    if args.out.exists() or (args.candidates_out and args.candidates_out.exists()):
        raise FileExistsError('Preserve existing candidate outputs')
    rows=build(json.loads(args.library.read_text()))
    with args.out.open('x') as stream:
        json.dump(dict(status='exact_designs_prepared_route_review_required', rdkit_version=rdBase.rdkitVersion,
                       candidates=rows),stream,indent=2);stream.write('\n')
    if args.candidates_out:
        with args.candidates_out.open('x') as stream:
            json.dump(rows,stream,indent=2);stream.write('\n')


if __name__ == '__main__':
    main()
