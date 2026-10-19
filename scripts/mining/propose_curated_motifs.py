#!/usr/bin/env python
"""Review prior leads and generate chemically explicit, uncomputed expansion proposals.

Writes only a fresh output directory. Does not register, build library files or launch QM.
RDKit graph checks establish a valid representation, not synthesis/cycling feasibility.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import pprint
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts' / 'mining'))
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from redox.build.contracts import validate_group
from redox.core.common import load_config, read_manifest
from redox.screening.discovery import describe, families, graft_products, identify, precursor_from_frag
from register_candidates import frag_from_grafted

SCAFFOLD = 'Cc1ccc(C[*:1])cc1'
RED2 = [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)]
RED1 = RED2[:2]
VIO2 = [('ox2', 2, 1, 0), ('ox1', 1, 2, -1), ('neu', 0, 1, -2)]
OX2 = [('neu', 0, 1, 0), ('ox1', 1, 2, 1), ('ox2', 2, 1, 2)]
SOURCES = {
 'pmdi': {'url':'https://pubs.acs.org/doi/10.1021/acsaem.1c01490', 'doi':'10.1021/acsaem.1c01490',
          'title':'Imide-Based Multielectron Anolytes as High-Performance Materials in Nonaqueous Redox Flow Batteries',
          'evidence_scope':'Parent monoimides support one reversible reduction; PMDI diimides support two. Grafted asymmetric caps are new designs.', 'access':'Primary publisher abstract inspected'},
 'ndi': {'url':'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011', 'doi':'10.1002/cssc.202501011',
         'title':'Ethylene Glycol-Functionalization of Naphthalene Diimides in Nonaqueous Redox Flow Batteries: Balance Between Solubility, Diffusivity, and Crossover',
         'evidence_scope':'Parent NDI two-electron chemistry and EG-chain trade-offs. Abstract describes nucleophilic substitution; no claim that asymmetric N-caps reproduce measured compounds.', 'access':'Primary abstract inspected; SI retrieval returned HTTP 403'},
 'quinone': {'url':'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632', 'doi':'10.1021/acs.jpcc.0c03632',
             'title':'Effect of Cycling Ion and Solvent on the Redox Chemistry of Substituted Quinones and Solvent-Induced Breakdown of the Correlation between Redox Potential and Electron-Withdrawing Power of Substituents',
             'evidence_scope':'Parent quinone reduction chemistry and strong ion/solvent dependence, not measurements of proposed grafted structures.', 'access':'Primary publisher abstract inspected'},
 'viologen': {'url':'https://www.sciencedirect.com/science/article/pii/S2095495617311403', 'doi':'10.1016/j.jechem.2018.02.014',
              'title':'Two electron utilization of methyl viologen anolyte in nonaqueous organic redox flow battery',
              'evidence_scope':'Methyl viologen two-electron utilization with TFSI in MeCN. No transfer of demonstrated performance to PF6 or asymmetric grafted ether derivatives.', 'access':'Primary publisher abstract inspected'},
 'glycol_viologen': {'url':'https://pubs.rsc.org/en/content/articlelanding/2022/ta/d2ta07876g', 'doi':'10.1039/D2TA07876G',
              'title':'A prototype of high-performance two-electron non-aqueous organic redox flow battery operated at -40 C',
              'evidence_scope':'Glycolated viologen/phenothiazine parent design rationale; exact proposed asymmetric benzyl derivatives unmeasured.', 'access':'Primary publisher abstract inspected'},
 'phenothiazine': {'url':'https://pubs.rsc.org/en/content/articlelanding/2017/ta/c7ta05883g', 'doi':'10.1039/C7TA05883G',
                  'title':'A stable two-electron-donating phenothiazine for application in nonaqueous redox flow batteries',
                  'evidence_scope':'3,7-methoxy donor substitution stabilizes a parent N-ethyl phenothiazine dication; proposed N-benzyl derivative still requires validation.', 'access':'Primary publisher abstract inspected'},
 'glycol_phenothiazine': {'url':'https://pubs.acs.org/doi/10.1021/acs.chemmater.8b04770', 'doi':'10.1021/acs.chemmater.8b04770',
                  'title':'Tailoring Two-Electron-Donating Phenothiazines To Enable High-Concentration Redox Electrolytes for Use in Nonaqueous Redox Flow Batteries',
                  'evidence_scope':'Specifically substituted 3,7-dialkoxy phenothiazines give two oxidation events. Does not justify two electrons for unsubstituted phenothiazines.', 'access':'Primary publisher abstract inspected'},
 'd3tales': {'url':'https://doi.org/10.18126/v5sj-6q93', 'doi':'10.18126/v5sj-6q93',
             'title':'D3TaLES public computational dataset', 'evidence_scope':'Local D3TaLES parent computed properties; no experimental evidence or production grafted energies.', 'access':'Previously downloaded local public dataset/pool'},
}


def group(gid, name, family, frag, source, notes, priority='series', states=RED2, origin='new expansion'):
    return dict(id=gid, name=name, family=family, frag=frag, states=states,
                provenance=dict(source_keys=source if isinstance(source,list) else [source],
                                type='source-inspired grafted design', measured_for_this_structure=False,
                                production_qm=False, origin=origin),
                chemistry_notes=notes, review_priority=priority)


def designs():
    out = []
    templates = {'pmdi':'[*:1]n1c(=O)c2cc3c(=O)n({cap})c(=O)c3cc2c1=O',
                 'ndi':'[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N({cap})C3=O'}
    caps = [('n_butyl','CCCC','butyl'),
            ('trioxyethyl','CCOCCOCCOC','2-[2-(2-methoxyethoxy)ethoxy]ethyl')]
    for core,template in templates.items():
        for capid,cap,name in caps:
            out.append(group(f'exp_{core}_{capid}', f'N-(4-methylbenzyl)-N\'-{name} {core.upper()}',
                             'imide (n-type)',template.format(cap=cap),core,
                             'Exactly one imide N-H precursor; isolate asymmetric mono-capped precursor before N-benzylation. Both imide units/four carbonyls retained. No residual donor proton. Longer chains add mass and conformers; no assumed solubility advantage in polymer phase.'))
    for stereo in ('@','@@'):
        out.append(group('exp_pmdi_pent2_'+('at' if stereo=='@' else 'atat'),
                         'N-(4-methylbenzyl)-N\'-(pentan-2-yl) PMDI, explicit stereoisomer '+stereo,
                         'imide (n-type)',templates['pmdi'].format(cap=f'[C{stereo}H](C)CCC'),'pmdi',
                         'Explicit pentan-2-yl stereoisomer; CIP identity exported from RDKit. Source parent uses branched pentyl substitution but proposed asymmetric tethered enantiomer is unmeasured. Isolate mono-N-H precursor; four imide carbonyls retained.'))
    out += [
        group('exp_aq_1oxy', '1-(4-methylbenzyloxy)anthracene-9,10-dione', 'quinone (n-type)',
              '[*:1]Oc1cccc2c1C(=O)c1ccccc1C2=O','quinone',
              'Single phenolic handle, two quinone carbonyls retained. O-alkylation removes peri OH hydrogen bonding; potential cannot be inherited from 1-hydroxy parent.', 'next chemistry smoke'),
        group('exp_aq_2oxy_6methoxy', '2-(4-methylbenzyloxy)-6-methoxyanthracene-9,10-dione', 'quinone (n-type)',
              '[*:1]Oc1ccc2c(c1)C(=O)c1ccc(OC)cc1C2=O','quinone',
              'Single non-peri phenolic handle and capped second oxygen. Preserve quinone carbonyls and avoid residual OH/PCET. Rational ether derivative; exact parent is not claimed as a measured paper compound.', 'next chemistry smoke'),
        group('exp_bq_oxy_methoxy', '2-(4-methylbenzyloxy)-5-methoxy-1,4-benzoquinone', 'quinone (n-type)',
              '[*:1]OC1=CC(=O)C(OC)=CC1=O','quinone',
              'Enolic hydroxyquinone precursor needs O/C alkylation selectivity check. Two carbonyls retained and no residual OH. Unshielded ring positions may undergo nucleophilic addition.', 'defer route review'),
        group('exp_bq_oxy_trimethyl', '2-(4-methylbenzyloxy)-3,5,6-trimethyl-1,4-benzoquinone', 'quinone (n-type)',
              '[*:1]OC1=C(C)C(=O)C(C)=C(C)C1=O','quinone',
              'Rational methyl shielding series, not a verified measured compound. Enolic OH precursor requires O/C selectivity review; methyl substitution does not prove cycling stability.', 'defer route review'),
        group('exp_viologen_methoxyethyl', 'N-(4-methylbenzyl)-N\'-(2-methoxyethyl)-4,4\'-bipyridinium', 'pyridine-multi-e',
              '[*:1][n+]1ccc(-c2cc[n+](CCOC)cc2)cc1',['viologen','glycol_viologen'],
              'Mono-quaternized bipyridine precursor has one remaining pyridine N; Menshutkin graft completes +2 viologen. Counterions excluded from structure. Second-wave accessibility and counterion dependence remain uncomputed.', 'next chemistry smoke',VIO2),
        group('exp_viologen_methoxyethoxyethyl', 'N-(4-methylbenzyl)-N\'-[2-(2-methoxyethoxy)ethyl]-4,4\'-bipyridinium', 'pyridine-multi-e',
              '[*:1][n+]1ccc(-c2cc[n+](CCOCCOC)cc2)cc1',['viologen','glycol_viologen'],
              'One remaining pyridine N after mono-quaternization; +2/+1/0 ladder. Ether chain is a deliberate polarity/mass control, not a solubility prediction.', 'series',VIO2),
        group('exp_phthalimide_control', 'N-(4-methylbenzyl)phthalimide', 'mono-imide (n-type)',
              '[*:1]N1C(=O)c2ccccc2C1=O','pmdi',
              'Single imide N-H gives one defined graft; two carbonyls retained. One-electron control only: do not assign two electrons merely from two carbonyls.', 'one-electron control',RED1),
        group('exp_naphthalimide_control', 'N-(4-methylbenzyl)-1,8-naphthalimide', 'mono-imide (n-type)',
              '[*:1]N1C(=O)c2cccc3cccc(c23)C1=O','pmdi',
              'Rational mono-imide core-size control, not a structure measured in cited PMDI paper. One-electron proposed ladder; reversibility for this enlarged mono-imide remains to test.', 'one-electron control',RED1),
        group('exp_ptz_37dimethoxy', 'N-(4-methylbenzyl)-3,7-dimethoxyphenothiazine', 'phenothiazine (p-type, two-electron proposal)',
              '[*:1]N1c2ccc(OC)cc2Sc2cc(OC)ccc21',['phenothiazine','glycol_phenothiazine'],
              'N-H alkylation preserves the two methoxy donors para to N that support parent dication stabilization. Explicit 0/+1/+2 oxidative ladder; second oxidation is proposed for N-benzyl analog. Generic discovery labels phenothiazine 1e, so this source-specific exception must remain explicit.', 'separate catholyte branch',OX2),
    ]
    return out


def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def all_proposals():
    prior = json.loads((ROOT/'results/discovery/imide_proposals_2026-10-10/proposals.json').read_text())
    for g in prior:
        core='pmdi' if '_pmdi_' in g['id'] else 'ndi'
        g['provenance']=dict(source_keys=[core], type='source-inspired grafted design', measured_for_this_structure=False,
                             production_qm=False, origin='prior nine-imide proposal review')
        g['review_priority']='next chemistry smoke' if g['id'] in ('design_pmdi_methoxyethyl','design_ndi_methyl','design_ndi_methoxyethyl') else 'series'
        g['chemistry_notes']='Single imide N-H precursor must be selectively prepared and isolated; indiscriminate mixed-amine dianhydride condensation can give symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N capped; no residual imide NH. Parent paper supports scaffold, not measurements of this benzyl structure.'
    pool=read_csv(ROOT/'results/discovery/review_2026-10-10/d3tales_pool.csv')
    ids={r['id'].removeprefix('d3_').upper() for r in read_csv(ROOT/'results/discovery/next_candidates_2026-10-10.csv') if r['id'].startswith('d3_')}
    ids |= {'80PRFF','80EEOZ','80GXWG','80TBWR','80NECN'}
    # One exact preliminary-front graft per source, not every competing graft site.
    front={r['source_id']:r['grafted_smiles'] for r in read_csv(ROOT/'results/discovery/review_2026-10-10/d3tales_prefilter_front.csv')}
    for r in pool:
        if r['source_id'] not in ids or r['status']!='ok' or (r['source_id'] in front and r['grafted_smiles'] != front[r['source_id']]):
            continue
        if not r['graft'].startswith('O-benzyl') and r['source_id'] in {'80PRFF','80EEOZ','80GXWG','80TBWR','80NECN'}:
            continue
        notes=[]
        pre=Chem.MolFromSmiles(r['parent_smiles'])
        product=Chem.MolFromSmiles(r['grafted_smiles'])
        if any(a.GetSymbol() in ('Cl','Br') for a in pre.GetAtoms()): notes.append('Ring halogen: competing nucleophilic substitution must be checked.')
        if r['graft'].startswith('N-benzyl (amine'): notes.append('Conjugated aminoquinone N is a weak nucleophile; modeled graft is not a verified direct Merrifield substitution route. Resin aminomethyl addition to quinone is an alternate architecture.')
        if pre.HasSubstructMatch(Chem.MolFromSmarts('[nX2;H0;+0]')): notes.append('Aza ring N competes through Menshutkin alkylation and can change core charge; attachment selectivity required.')
        if r['graft'].startswith('benzyl ester'): notes.append('Ester linkage cleavage and reactive o-quinone chemistry require stability review.')
        if r['graft'].startswith('O-benzyl'):
            if pre.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]-[C]=[C]')): notes.append('Enolic OH: O/C alkylation selectivity required; both carbonyls must survive.')
            else: notes.append('Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact.')
        if Descriptors.NumHDonors(product): notes.append('Residual OH/NH donor remains; proton transfer and competitive attachment require review.')
        if r['source_id']=='80NECN': notes.append('Conjugated tertiary amine may compete by quaternization; single phenolic handle is not proof of selectivity.')
        if r['source_id']=='80PRFF': notes.append('Peri phenolic O graft removes OH-to-carbonyl H-bond; recompute grafted properties.')
        notes.append('D3TaLES lambda is parent computed evidence only; no second-step grafted lambda or potential is known.')
        priority='next chemistry smoke' if r['source_id'] in {'80PRFF','80EEOZ'} else 'defer route review'
        if r['source_id'] in {'80GXWG','80TBWR'}: priority='series'
        g=group('d3_'+r['source_id'].lower(),'D3TaLES '+r['source_id']+' '+r['graft']+' product (canonical structure supplied)',
                'quinone (n-type)',frag_from_grafted(r['grafted_smiles']),['d3tales','quinone'],' '.join(notes),priority,
                origin='prior ten-front review' if r['source_id'] in front else 'additional D3TaLES pool candidate')
        g['provenance'].update(source_id=r['source_id'],parent_smiles=r['parent_smiles'],d3_parent_electron_reorg_eV=r.get('d3_electron_reorg_eV'))
        prior.append(g)
    return prior+designs()


def validate(g, known, seen, faraday):
    m=validate_group(SCAFFOLD,g)
    smi=Chem.MolToSmiles(m)
    key=Chem.MolToInchiKey(m)
    if key in known or key in seen:
        return None
    seen.add(key)
    pre=precursor_from_frag(g['frag'])
    p=Chem.MolFromSmiles(pre)
    prods=graft_products(p)
    matches=[chem for gsmi,chem in prods if Chem.MolToInchiKey(Chem.MolFromSmiles(gsmi))==key]
    if not matches: raise ValueError(f"{g['id']}: precursor does not regenerate intended graft")
    carbonyl=Chem.MolFromSmarts('[#6]=[O]')
    before=len(p.GetSubstructMatches(carbonyl)); after=len(m.GetSubstructMatches(carbonyl))
    if before != after: raise ValueError(f"{g['id']}: carbonyl lost during graft")
    fam=families(m)
    n_e=max(abs(t[3]) for t in g['states'])
    stats=describe(smi,n_e,faraday)
    exact_nominal=n_e*faraday/(3.6*Descriptors.MolWt(m))
    assert abs(stats['capacity_nominal_mAh_g']-exact_nominal)<=0.051
    if 'imide' in g['family'] and after != (2 if n_e==1 else 4): raise ValueError(f"{g['id']}: wrong imide carbonyl count")
    if 'quinone' in g['family'] and not any('quinone' in f for f in fam): raise ValueError(f"{g['id']}: quinone core missing")
    if 'pyridine-multi-e' in g['family'] and 'viologen' not in fam: raise ValueError(f"{g['id']}: viologen core missing")
    mh=Chem.AddHs(m)
    ep=AllChem.ETKDGv3(); ep.randomSeed=20261010; ep.numThreads=1
    embedded=AllChem.EmbedMolecule(mh,ep)==0
    if not embedded: raise ValueError(f"{g['id']}: embedding failed")
    stereo=Chem.FindMolChiralCenters(m,includeUnassigned=True,useLegacyImplementation=False)
    rule_rows=identify([dict(source_id=g['id'],smiles=pre)],faraday)
    statuses=[r['status'] for r in rule_rows if r.get('grafted_smiles') and Chem.MolToInchiKey(Chem.MolFromSmiles(r['grafted_smiles']))==key]
    return dict(id=g['id'],name=g['name'],family=g['family'],origin=g['provenance']['origin'],
                priority=g['review_priority'],frag=g['frag'],precursor_smiles=pre,grafted_smiles=smi,inchikey=key,
                formula=rdMolDescriptors.CalcMolFormula(m),states=json.dumps(g['states']),
                source_urls='; '.join(SOURCES[s]['url'] for s in g['provenance']['source_keys']),
                measured_for_grafted_structure=False,production_qm=False,**stats,
                HBD=Descriptors.NumHDonors(m), carbonyl_before=before,carbonyl_after=after,
                recognized_families='+'.join(fam),graft='; '.join(matches),distinct_precursor_graft_products=len(prods),
                exact_graft_roundtrip=True,fully_assigned_stereo=True,stereocenters=json.dumps(stereo),
                embedded=True,default_discovery_status='; '.join(statuses),chemistry_notes=g['chemistry_notes'],
                d3_parent_electron_reorg_eV=g['provenance'].get('d3_parent_electron_reorg_eV',''))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',required=True,type=Path)
    a=ap.parse_args()
    if a.out.exists(): ap.error('choose a fresh output directory; prior proposals preserved')
    faraday=load_config('electrolyte').FARADAY
    known={Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles'])) for r in read_manifest()}
    for conf in ('redox_groups','starting_candidates','merrifield_multielectron','discovered_candidates','custom_candidates'):
        for g in load_config(conf).GROUPS:
            known.add(Chem.MolToInchiKey(validate_group(SCAFFOLD,g)))
    groups=[]; rows=[]; excluded=[]; seen=set()
    for g in all_proposals():
        r=validate(g,known,seen,faraday)
        if r is None: excluded.append(dict(id=g['id'],reason='existing library/config or duplicate proposal')); continue
        if g['id'].startswith('exp_pmdi_pent2_'):
            cip=json.loads(r['stereocenters'])[0][1]
            g['id']=r['id']='exp_pmdi_pent2_'+cip.lower()
            g['name']=r['name']=f"N-(4-methylbenzyl)-N'-({cip})-pentan-2-yl PMDI"
        groups.append(g); rows.append(r)
    a.out.mkdir(parents=True)
    (a.out/'proposals.json').write_text(json.dumps(groups,indent=2)+'\n')
    with (a.out/'curated_candidates.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (a.out/'proposals.py').write_text('"""Uncomputed curated designs; chemistry review before registration/production."""\n'+f'SCAFFOLD = {SCAFFOLD!r}\nGROUPS = '+pprint.pformat(groups,width=110,sort_dicts=False)+'\n')
    (a.out/'sources.json').write_text(json.dumps(SOURCES,indent=2)+'\n')
    counts={origin:sum(r['origin']==origin for r in rows) for origin in sorted({r['origin'] for r in rows})}
    inputs=['results/discovery/imide_proposals_2026-10-10/proposals.json','results/discovery/review_2026-10-10/d3tales_pool.csv','results/discovery/review_2026-10-10/d3tales_prefilter_front.csv','results/discovery/next_candidates_2026-10-10.csv','library/manifest.csv']
    summary=dict(count=len(rows),origins=counts,excluded=excluded,all_graph_graft_stereo_embedding_checks_pass=True,
                 rdkit_version=rdBase.rdkitVersion,faraday_C_mol=faraday,
                 input_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
                 limitation='Graph/stereo/embedding/arithmetic validation only. No property predictions or synthetic/stability verification.',
                 counterion_convention='Nominal capacity uses grafted active species mass; no counterions or backbone-repeat correction.',
                 quantum_calculations_submitted=0)
    (a.out/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
