#!/usr/bin/env python
"""Prepare a small candidate batch with explicit route/provenance gates; never submit QM."""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import pprint
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from rdkit import Chem
from redox.build.contracts import validate_group
from redox.core.common import load_config,read_manifest
from redox.screening.discovery import describe,graft_products,precursor_from_frag
SCAFFOLD='Cc1ccc(C[*:1])cc1'
SELECTED=['exp_viologen_methoxyethyl','design_pmdi_methoxyethyl','exp_ndi_hexyl_pentyloxy',
          'd3_80prff','exp_aq_1oxy','exp_phthalimide_control']
DEFERRED=['design_ndi_methyl','design_ndi_methoxyethyl','d3_80eeoz','exp_aq_2oxy_6methoxy']
# Each evidence record separates actual experimental source structures from new grafted designs.
EVIDENCE={
'exp_viologen_methoxyethyl':dict(
 precursor_name="1-(2-methoxyethyl)-4,4'-bipyridinium (monoquaternized)",
 experimental_parent_name='methyl viologen (exact redox parent); branched-alkyl styryl viologen (attachment precedent)',
 experimental_parent_smiles='C[n+]1ccc(-c2cc[n+](C)cc2)cc1',
 experimental_parent_reference='10.1016/j.jechem.2018.02.014',
 route_source='https://patents.google.com/patent/US8773745B2/en',route_locus='Synthetic Example2: exact monoquaternized precursor and HBPS-EC-10',
 route_additional_source='https://onlinelibrary.wiley.com/doi/10.1002/macp.202100373',
 proposed_route=['Prepare 2-methoxyethyl tosylate, mono-quaternize 4,4-prime-bipyridine, and chromatographically isolate the mono-cation as reported in US8773745B2 Synthetic Example2.',
                 'Quaternize its single remaining pyridine N using the Merrifield benzylic chloride or 4-methylbenzyl chloride proxy.',
                 'Perform explicit counterion exchange for the target PF6 electrolyte; do not inherit TFSI-cell behavior.'],
 route_status='exact ether-cap precursor and benzylic-bromide polymer graft published; chloride transfer proposed',
 direct_merrifield_route_verified=False,precursor_exact_primary_synthesis_verified=True,
 route_note='US8773745B2 publishes the exact monoquaternized ether-cap precursor and its attachment to benzylic-bromide polystyrene. Borchers2022 demonstrates chloromethylstyrene attachment with a different cap. Ether-cap transfer to Merrifield chloride remains to establish.',
 positions='4,4-prime inter-ring bond; one tether on each pyridinium N'),
'design_pmdi_methoxyethyl':dict(
 precursor_name="N-(2-methoxyethyl)-N'-H pyromellitic diimide",
 experimental_parent_name="N,N'-bis(2-methoxyethyl) PMDI, Daub compound2e",
 experimental_parent_smiles='COCCn1c(=O)c2cc3c(=O)n(CCOC)c(=O)c3cc2c1=O',
 experimental_parent_reference='10.1021/acsaem.1c01490',
 route_source='https://acs.figshare.com/articles/journal_contribution/16684876',route_locus='SI S5-S6 compound2e; S7 asymmetric2i',
 proposed_route=['A direct chloride route requires preparing and isolating the N-(2-methoxyethyl)-N-prime-H diimide, then deprotonating its only imide NH and benzylating with resin.',
                 'The inspected SI does not establish that mono-N-H precursor route; do not assume selective half-alkylation.',
                 'Alternative to the same product connectivity: prepare/isolate a mono-methoxyethyl imide/anhydride, then imidize the remaining anhydride with a Merrifield-derived aminomethyl support; use 4-methylbenzylamine for the molecular proxy.'],
 route_status='direct mono-N-H chloride route unresolved; staged imidization has an experimental asymmetric analog',
 direct_merrifield_route_verified=False,precursor_exact_primary_synthesis_verified=False,
 route_note='Daub2i uses staged primary-amine addition and chromatography, not an isolated mono-N-H diimide. Its yield cannot be transferred to this structure. Aminomethyl support preparation is an additional route step.',
 positions='Pyromellitic benzene-1,2:4,5 bis-imide; attachment through one imide N and methoxyethyl on the other'),
'exp_ndi_hexyl_pentyloxy':dict(
 precursor_name="N-(5-hydroxypentyl)-N'-n-hexyl naphthalene-1,8:4,5-bis(dicarboximide)",
 experimental_parent_name="N-(5-hydroxypentyl)-N'-n-hexyl NDI, hydroxyl-functionalized intermediate of Al Kurdi2021",
 experimental_parent_smiles='OCCCCCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCCCCC)C3=O',
 experimental_parent_reference='10.1039/D0QM00685H',
 route_source='https://pubs.rsc.org/en/content/articlehtml/2021/qm/d0qm00685h',route_locus='Abstract and Results/Synthesis Scheme1',
 proposed_route=['Use the published asymmetric-imidization route: NTDA to the N-hexyl monoimide dicarboxylic acid, then imidize with 5-amino-1-pentanol.',
                 'Isolate the fully N-capped NDI bearing one terminal primary alcohol.',
                 'O-benzylate that alcohol using Merrifield chloride or the molecular 4-methylbenzyl chloride proxy; establish conditions that avoid imide hydrolysis.'],
 route_status='exact isolated hydroxyalkyl precursor published; terminal-alcohol chloride graft proposed',
 direct_merrifield_route_verified=False,precursor_exact_primary_synthesis_verified=True,
 route_note='This removes the unestablished mono-N-H diimide requirement. The paper esterifies the OH with norbornenecarboxylic acid; it does not report this ether graft. Its redox data are in DCM and are not MeCN property anchors.',
 positions='Naphthalene-1,8:4,5 diimide; one N-n-hexyl, other N-(CH2)5-O-CH2-aryl tether; no free imide NH'),
'd3_80prff':dict(
 precursor_name='5-hydroxy-1,4-naphthoquinone (juglone)',
 experimental_parent_name='juglone; 5-benzyloxy-1,4-naphthoquinone, Brimble2001 intermediate15',
 experimental_parent_smiles='O=C1C=CC(=O)c2c(O)cccc21',
 experimental_parent_reference='10.1039/B103080A',
 route_source='https://pubs.rsc.org/en/content/articlelanding/2001/p1/b103080a',route_locus='Scheme2 and discussion of juglone14 to benzyl ether15',
 proposed_route=['Start from isolated juglone; preserve both quinone C=O groups.',
                 'Use the published mild Ag2O/benzyl bromide phenolic benzylation as soluble-proxy precedent.',
                 'For actual resin chloride, establish comparable mild graft conditions experimentally rather than assuming a strong-base Williamson procedure.'],
 route_status='exact parent phenolic benzyl ether has primary synthesis precedent; p-methyl/resin-chloride transfer proposed',
 direct_merrifield_route_verified=False,precursor_exact_primary_synthesis_verified=True,
 route_note='Parent benzylation uses BnBr/Ag2O/CHCl3; changing Br to resin-bound chloride can change conversion. Removal of peri OH H-bond changes redox potentials. Local juglone electrochemical benchmark contains peak potentials, not reliable E1/2 anchors.',
 positions='5-phenolic O tether on 1,4-naphthoquinone fused benzene ring; not lawsone2-enolic O'),
'exp_aq_1oxy':dict(
 precursor_name='1-hydroxyanthracene-9,10-dione',
 experimental_parent_name='1-hydroxyanthraquinone (benchmark parent); nitro-substituted hydroxy-AQ benzyl ether route analog',
 experimental_parent_smiles='Oc1cccc2c1C(=O)c1ccccc1C2=O',
 experimental_parent_reference='10.1021/acs.jpcc.2c03129',
 route_source='https://patents.google.com/patent/WO2004084891A1/en',route_locus='Example1F hydroxyanthraquinone benzyl ether starting material',
 proposed_route=['Use isolated 1-hydroxyanthraquinone (one phenol, no amine).',
                 'Screen phenolic O-benzylation with a benzylic halide; the primary patent reports K2CO3/DMF with benzyl bromide for a nitro-substituted hydroxy-AQ analog.',
                 'Adapt to 4-methylbenzyl chloride/resin only after verifying intact quinone carbonyls and O-selectivity.'],
 route_status='hydroxy-AQ benzyl ether reaction-class precedent; exact 1-O resin route proposed',
 direct_merrifield_route_verified=False,precursor_exact_primary_synthesis_verified=True,
 route_note='The cited reaction analog has a nitro substituent absent from this target. 1-O is peri to a carbonyl; this changes intramolecular H-bonding relative to parent. Parent benchmark half-peak potentials are not exact thermodynamic E values.',
 positions='1-O tether; central 9,10 anthraquinone carbonyls retained'),
'exp_phthalimide_control':dict(
 precursor_name='phthalimide',experimental_parent_name='N-methyl and N-(2-methoxyethyl)phthalimides, Daub1a/1e',
 experimental_parent_smiles='Cn1c(=O)c2ccccc2c1=O',experimental_parent_reference='10.1021/acsaem.1c01490',
 route_source='https://pubs.rsc.org/en/content/articlehtml/2019/ra/c9ra01131e',route_locus='Section2.2: phthalimide substitution onto poly(4-chloromethylstyrene)',
 route_additional_source='https://acs.figshare.com/articles/journal_contribution/16684876',
 proposed_route=['Start with phthalimide and deprotonate its sole NH to the imide nucleophile.',
                 'Displace the benzylic chloride on Merrifield or 4-methylbenzyl chloride to give one N-benzyl imide.',
                 'Keep the control at one declared electron; a second carbonyl does not justify a second electron.'],
 route_status='direct phthalimide substitution onto chloromethyl-polystyrene published; exact truncated proxy properties uncomputed',
 direct_merrifield_route_verified=True,precursor_exact_primary_synthesis_verified=True,
 route_note='RSC2019 C9RA01131E reports potassium phthalimide substitution onto poly(4-chloromethylstyrene). This flag means a published direct support-reaction precedent, not experimental execution here or validated electrochemical performance. Daub provides monoimide one-electron evidence.',
 positions='One benzene-1,2 dicarboximide; graft on its sole imide N'),
}


def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--out',required=True,type=Path)
 ap.add_argument('--curated',type=Path,default=ROOT/'results/discovery/parallel_expansion_2026-10-10/proposals.json')
 ap.add_argument('--only',help='comma-separated IDs from the balanced six-motif selection')
 a=ap.parse_args();out=a.out.resolve()
 for name in ('batch_config.py','batch_config.json','screening_inputs.csv','route_evidence.json','preparation_summary.json'):
  if (out/name).exists():ap.error(f'output already exists: {out/name}')
 byid={g['id']:g for g in json.loads(a.curated.read_text())}
 byid['exp_ndi_hexyl_pentyloxy']=dict(id='exp_ndi_hexyl_pentyloxy',
    name="N-[5-(4-methylbenzyloxy)pentyl]-N'-n-hexyl naphthalene diimide",family='imide (n-type)',
    frag='[*:1]OCCCCCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCCCCC)C3=O',
    states=[('neu',0,1,0),('red1',-1,2,-1),('red2',-2,1,-2)])
 groups=[];rows=[];already_registered=[];faraday=load_config('electrolyte').FARADAY
 known={Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles'])) for r in read_manifest()}
 for conf in ('redox_groups','starting_candidates','merrifield_multielectron',
              'discovered_candidates','custom_candidates'):
  for registered in load_config(conf).GROUPS:
   known.add(Chem.MolToInchiKey(validate_group(SCAFFOLD,registered)))
 selected=list(SELECTED)
 if a.only:
  requested=set(a.only.split(','))
  if requested-set(SELECTED):ap.error(f'unknown shortlist IDs: {sorted(requested-set(SELECTED))}')
  selected=[gid for gid in SELECTED if gid in requested]
 seen=set()
 for gid in selected:
  g=copy.deepcopy(byid[gid]);e=EVIDENCE[gid]
  g['provenance']=dict(type='uncomputed grafted proposal',source=e['route_source'],
     exact_grafted_experiment=False,production_qm=False,route_status=e['route_status'],
     direct_merrifield_route_verified=e['direct_merrifield_route_verified'],
     experimental_parent_reference=e['experimental_parent_reference'])
  g['chemistry_notes']=e['route_note'];g['review_priority']='balanced first calculation batch; physical protocol gate pending'
  m=validate_group(SCAFFOLD,g);key=Chem.MolToInchiKey(m)
  if key in known:
   already_registered.append(dict(id=gid,inchikey=key,reason='exact structure already registered; no duplicate proposal'))
   continue
  if key in seen:raise ValueError(f'{gid}: duplicate candidate')
  seen.add(key);smi=Chem.MolToSmiles(m);pre=precursor_from_frag(g['frag']);p=Chem.MolFromSmiles(pre)
  prods=graft_products(p)
  chem=next((chem for sm,chem in prods if Chem.MolToInchiKey(Chem.MolFromSmiles(sm))==key),None)
  if chem is None:raise ValueError(f'{gid}: exact graft round trip failed')
  # Keep witness structures separate; no source energies, kinetics or solubility copied into inputs.
  witness=Chem.MolFromSmiles(e['experimental_parent_smiles'])
  if witness is None:raise ValueError(f'{gid}: invalid experimental witness structure')
  n=max(abs(s[3]) for s in g['states'])
  groups.append(g);rows.append(dict(id=gid,name=g['name'],family=g['family'],frag=g['frag'],
     precursor_smiles=pre,grafted_smiles=smi,inchikey=key,states=json.dumps(g['states']),
     graft=chem,graph_roundtrip=True,**describe(smi,n,faraday),
     route_status=e['route_status'],exact_grafted_experiment=False,production_qm=False,
     experimental_parent_smiles=Chem.MolToSmiles(witness),experimental_parent_reference=e['experimental_parent_reference']))
 if not groups:ap.error('all selected structures are registered; no new candidate inputs to prepare')
 out.mkdir(parents=True,exist_ok=True)
 (out/'batch_config.py').write_text('"""Isolated screening inputs. New grafted properties uncomputed; route gates remain explicit."""\n'+f'SCAFFOLD = {SCAFFOLD!r}\nGROUPS = '+pprint.pformat(groups,width=100,sort_dicts=False)+'\n')
 (out/'batch_config.json').write_text(json.dumps(groups,indent=2)+'\n')
 with (out/'screening_inputs.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (out/'route_evidence.json').write_text(json.dumps(EVIDENCE,indent=2)+'\n')
 summary=dict(ids=[g['id'] for g in groups],molecules=len(groups),declared_states=sum(len(g['states']) for g in groups),
     already_registered=already_registered,
     excluded_from_first_batch=DEFERRED,all_computed_electrochemical_metrics_missing=True,
     exact_grafted_experimental_validation=False,all_direct_resin_routes_unverified=not any(g['provenance']['direct_merrifield_route_verified'] for g in groups),
     published_direct_resin_precedent_ids=[g['id'] for g in groups if g['provenance']['direct_merrifield_route_verified']],
     input_sha256=hashlib.sha256(a.curated.read_bytes()).hexdigest(),quantum_jobs_submitted=0,
     isolated_library='library',isolated_builder_command=f'python -m redox.build.candidates --config {out}/batch_config.py --library {out}/library',
     dryrun_command=f'python -m redox.build.candidates --config {out}/batch_config.py --library {out}/library --dry-run')
 (out/'preparation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
