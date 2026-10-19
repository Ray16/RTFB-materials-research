#!/usr/bin/env python
"""Prepare isolated vinyl-monomer graphs without inheriting capped-model properties.

The synthesis route is redox functionalization of 4-vinylbenzyl chloride followed
by polymerization. No shared manifest, configuration, conformers or QM is changed.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from redox.build.contracts import validate_group
from redox.core.common import ROOT
from redox.properties.capacity_and_proxies import F_C_PER_MOL, _sa_prep, _sascorer

CAPPED_SCAFFOLD='Cc1ccc(C[*:1])cc1'
VINYL_SCAFFOLD='C=Cc1ccc(C[*:1])cc1'
VINYL_HANDLE=Chem.MolFromSmarts('[CH2]=[CH]-[c]1[cH][cH][c]([CH2]-[!#0])[cH][cH]1')
PILOT_IDS={'design_ndi_ethyl','design_pmdi_ethyl','ethylviologen'}


def canonical(mol):
    return Chem.MolToSmiles(mol,isomericSmiles=True)


def prepare(row):
    states=json.loads(row['states']) if isinstance(row['states'],str) else row['states']
    group=dict(id='vm_'+row['id'],name='Vinyl counterpart of '+row['name'],
               family=row['family'],frag=row['frag'],states=states)
    capped=validate_group(CAPPED_SCAFFOLD,group)
    if canonical(capped)!=canonical(Chem.MolFromSmiles(row['exact_grafted_smiles'])):
        raise ValueError('Snapshot fragment does not reconstruct its stated capped model: '+row['id'])
    monomer=validate_group(VINYL_SCAFFOLD,group)
    matches=monomer.GetSubstructMatches(VINYL_HANDLE)
    if len(matches)!=1:raise ValueError('Require exactly one specified para-vinylbenzyl handle')
    # Independent inverse graph operation: removing only the terminal vinyl carbon
    # must regenerate the complete original capped model, including stereochemistry.
    inverse=Chem.RWMol(monomer);inverse.RemoveAtom(matches[0][0]);inverse=inverse.GetMol()
    Chem.SanitizeMol(inverse)
    if canonical(inverse)!=canonical(capped):raise ValueError('Vinyl-to-cap inverse graph check failed')
    before=Counter(a.GetSymbol() for a in Chem.AddHs(capped).GetAtoms())
    after=Counter(a.GetSymbol() for a in Chem.AddHs(monomer).GetAtoms())
    expected=before.copy();expected['C']+=1
    if after!=expected:raise ValueError('Transformation must add C1 with no other elemental changes')
    mw=Descriptors.MolWt(monomer);cap_mw=Descriptors.MolWt(capped)
    if abs(mw-cap_mw-12.011)>1e-7:raise ValueError('Unexpected mass change')
    smiles=canonical(monomer);sa_mol,_,sa_note=_sa_prep(smiles)
    electrons=max(s[3] for s in states)-min(s[3] for s in states)
    is_quinone='quinone' in row['family']
    return dict(id=group['id'],name=group['name'],source_capped_id=row['id'],family=row['family'],
        smiles=smiles,inchikey=Chem.MolToInchiKey(monomer),formula=rdMolDescriptors.CalcMolFormula(monomer),
        states=states,declared_ladder_status='proposal inherited from design; vinyl-state properties uncomputed',
        MW_g_mol=mw,n_e_nominal=electrons,capacity_nominal_mAh_g=electrons*F_C_PER_MOL/(3.6*mw),
        SA_score=_sascorer().calculateScore(sa_mol) if sa_mol is not None else None,sa_note=sa_note,
        source_capped_smiles=canonical(capped),source_capped_MW_g_mol=cap_mw,
        carbon_added=1,other_element_counts_unchanged=True,cap_inverse_graph_verified=True,
        exact_grafted_experimental_validation=False,vinyl_monomer_synthesis_verified=False,
        polymerization_verified=False,polymer_redox_performance_verified=False,
        initial_modeling_pilot=row['id'] in PILOT_IDS,
        initial_pilot_reason=('Conventional diimide / viologen route review and model-fidelity controls'
                              if row['id'] in PILOT_IDS else None),
        polymerization_review=('Quinone radical trapping / inhibition requires explicit review'
                               if is_quinone else 'Initiation conditions, conversion and redox-center retention unverified'),
        attachment_notes=row.get('chemistry_notes'),
        E1_vs_Fc_V=None,E2_vs_Fc_V=None,lambda_solution_eV=None,
        inherited_computed_properties=False,counterion='not assigned; charge-bearing component only')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--feed',type=Path,default=ROOT/'results/candidates/current.json')
    parser.add_argument('--out',type=Path,default=ROOT/'results/discovery/vinyl_monomers')
    args=parser.parse_args();payload=json.loads(args.feed.read_text())
    source=payload['original39']+payload['additional_candidates']+payload.get('comparison_controls',[])
    manifest=ROOT/'library/manifest.csv';manifest_before=hashlib.sha256(manifest.read_bytes()).hexdigest()
    rows=[prepare(row) for row in source]
    if len({r['id'] for r in rows})!=len(rows):raise ValueError('Duplicate proposal IDs')
    molecules={}
    for row in rows:molecules.setdefault(row['smiles'],[]).append(row['id'])
    out=args.out;out.mkdir(parents=True,exist_ok=False)
    with (out/'proposals.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    proposal=dict(status='graph_and_descriptor_preparation_only',
        synthesis_route='Attach redox center to para-vinylbenzyl monomer, then polymerize',
        scaffold=VINYL_SCAFFOLD,assumption='Pure para isomer; other substitution isomers require separate proposals',
        source_snapshot=payload['snapshot_path'],source_as_of_Chicago=payload['as_of_Chicago'],
        designs=len(payload['original39'])+len(payload['additional_candidates']),controls=len(payload.get('comparison_controls',[])),
        proposals=len(rows),unique_molecular_graphs=len(molecules),aliases=[ids for ids in molecules.values() if len(ids)>1],
        initial_modeling_pilot=[r['id'] for r in rows if r['initial_modeling_pilot']],
        all_inverse_graph_checks_passed=True,all_element_and_charge_ladder_checks_passed=True,
        shared_library_manifest_unchanged=hashlib.sha256(manifest.read_bytes()).hexdigest()==manifest_before,
        production_QM_registered=False,quantum_jobs_submitted=0,
        scope='Nominal capacities and SA are graph descriptors of ionic components where charged. '
              'Counterions, comonomer fraction, polymer loading and cycling accessibility are unqualified. '
              'No capped-model redox or reorganization value is assigned to a vinyl monomer.',
        source_sha256={str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       str(manifest):manifest_before},records=rows)
    assert proposal['shared_library_manifest_unchanged']
    (out/'proposals.json').write_text(json.dumps(proposal,indent=2)+'\n')
    print(json.dumps({k:v for k,v in proposal.items() if k not in ('records','source_sha256','scope')}))


if __name__=='__main__':main()
