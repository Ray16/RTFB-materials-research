#!/usr/bin/env python
"""Independent read-only audit of candidate registration, structures and pilot semantics.

Writes only a fresh report folder. No production builds, energy calculations or registry edits.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds, rdMolDescriptors
from redox.build.contracts import validate_group
from redox.core.common import load_config,read_manifest,read_result,UMA,DFT,_skeleton
from redox.properties.potentials import uma_gas_diagnostic
from redox.screening.discovery import describe,precursor_from_frag,graft_products
from redox.screening.scorecard import _candidate_chem_flags
from redox.validation.readiness import inspect_states
CONFIGS=('redox_groups','starting_candidates','merrifield_multielectron','discovered_candidates','custom_candidates')
PILOT_IDS=('exp_viologen_methoxyethyl','exp_phthalimide_control')


def registry_audit():
    rows=read_manifest();seen_ids={};seen_keys={};checks=[];defs=[]
    for conf in CONFIGS:
        cfg=load_config(conf)
        for g in cfg.GROUPS:
            m=validate_group(cfg.SCAFFOLD,g);smi=Chem.MolToSmiles(m);key=Chem.MolToInchiKey(m)
            if g['id'] in seen_ids and seen_ids[g['id']]!=smi:raise ValueError(f"conflicting registered ID: {g['id']}")
            seen_ids[g['id']]=smi; seen_keys.setdefault(key,[]).append((conf,g['id']))
            actual=[r for r in rows if r['id']==g['id']]
            declared={tuple(s) for s in g['states']}
            recorded={(r['state'],int(r['charge']),int(r['mult_hint']),int(r['n_e'])) for r in actual}
            check=dict(id=g['id'],config=conf,manifest_present=bool(actual),states_match=recorded==declared,
               smiles_match=bool(actual) and all(Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles']))==key for r in actual),
               formula=rdMolDescriptors.CalcMolFormula(m),inchikey=key,declared_states=len(declared),
               metadata_source=g.get('provenance',{}).get('source',''),
               chemistry_notes=bool(g.get('chemistry_notes')),scorecard_route_note_present=bool(_candidate_chem_flags().get(g['id'])))
            checks.append(check);defs.append((g,m))
    aliases=[dict(inchikey=k,entries=v) for k,v in seen_keys.items() if len({i for _,i in v})>1]
    return rows,checks,defs,aliases


def pilot_audit(rows):
    out=[];faraday=load_config('electrolyte').FARADAY
    for gid in PILOT_IDS:
        mr=[r for r in rows if r['id']==gid]
        for r in mr:
            st=r['state'];sd=UMA/gid/st;d=json.loads((sd/'result.json').read_text());xyz=sd/'relaxed.xyz'
            mol=Chem.AddHs(Chem.MolFromSmiles(r['smiles']))
            # ASE's warm-start XYZ may append force columns. Independently normalize
            # coordinates to plain XYZ in memory before RDKit graph perception.
            lines=xyz.read_text().splitlines()
            plain='\n'.join(lines[:2]+[' '.join(line.split()[:4]) for line in lines[2:]])+'\n'
            actual=Chem.MolFromXYZBlock(plain)
            if actual is None:raise ValueError(f'{gid}/{st}: unreadable coordinates')
            rdDetermineBonds.DetermineConnectivity(actual,charge=0)
            ref=_skeleton(Chem.RemoveHs(mol));observed=_skeleton(Chem.RemoveHs(actual,sanitize=False))
            nuclei=sum(a.GetAtomicNum() for a in actual.GetAtoms());q=int(r['charge']);mult=d['chosen_mult']
            delta_nuclei=nuclei-q
            p=DFT/gid/st
            check=dict(id=gid,state=st,model=d.get('model'),converged=d.get('converged') is True,
                record_smiles_matches=d.get('smiles')==r['smiles'],charge_matches=d.get('charge')==q,
                declared_n_e_matches=d.get('n_e')==int(r['n_e']),electron_count_matches=d.get('n_elec')==delta_nuclei,
                spin_parity_valid=(delta_nuclei-(mult-1))%2==0,
                composition_matches=sorted(a.GetAtomicNum() for a in actual.GetAtoms())==sorted(a.GetAtomicNum() for a in mol.GetAtoms()),
                connectivity_matches=observed.HasSubstructMatch(ref) and ref.HasSubstructMatch(observed),
                geometry_sha256=hashlib.sha256(xyz.read_bytes()).hexdigest(),
                result_sha256=hashlib.sha256((sd/'result.json').read_bytes()).hexdigest(),
                checkpoint_hash_recorded=bool(d.get('checkpoint_sha256')),geometry_hash_recorded=bool(d.get('geom_sha256')),
                production_dft_result_present=(p/'result.json').exists(),
                production_energy_available=(read_result(gid,st) or {}).get('e_smd_eV') is not None)
            out.append(check)
    return out


def descriptor_audit(defs):
    f=load_config('electrolyte').FARADAY;out=[]
    for g,m in defs:
        if g['id'] not in PILOT_IDS:continue
        smi=Chem.MolToSmiles(m);n=max(abs(s[3]) for s in g['states']);stats=describe(smi,n,f)
        atom_mass=sum(a.GetMass() for a in Chem.AddHs(m).GetAtoms());pre=precursor_from_frag(g['frag'])
        regenerated=[s for s,_ in graft_products(Chem.MolFromSmiles(pre))]
        out.append(dict(id=g['id'],grafted_smiles=smi,precursor_smiles=pre,declared_nominal_electrons=n,
            MW_atomic_sum=round(atom_mass,6),MW_matches=abs(stats['MW']-atom_mass)<.00051,
            nominal_capacity_atomic_sum=round(n*f/(3.6*atom_mass),6),
            nominal_capacity_matches=abs(stats['capacity_nominal_mAh_g']-n*f/(3.6*atom_mass))<.051,
            precursor_graft_matches=any(Chem.MolToInchiKey(Chem.MolFromSmiles(s))==Chem.MolToInchiKey(m) for s in regenerated),
            **stats,convention='exact grafted active-species mass; no counterions or polymer-repeat correction'))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args()
    if a.out.exists():ap.error('choose a fresh audit directory')
    rows,registry,defs,aliases=registry_audit();pilot=pilot_audit(rows);desc=descriptor_audit(defs)
    ready=inspect_states([r for r in rows if r['id'] in PILOT_IDS])
    gas=[]
    for gid in PILOT_IDS:
        mr=sorted([r for r in rows if r['id']==gid],key=lambda r:-int(r['charge']))
        for o,r in zip(mr,mr[1:]):
            delta,E,qc=uma_gas_diagnostic(gid,o['state'],r['state'],int(o['charge']),int(r['charge']))
            gas.append(dict(id=gid,event=o['state']+'->'+r['state'],dE_uma_gas_eV=delta,E_vs_UMA_Fc_gas_V=E,
                            diagnostic_qc=qc,production_solvated_property=False))
    state_fields=('converged','record_smiles_matches','charge_matches','declared_n_e_matches','electron_count_matches',
                  'spin_parity_valid','composition_matches','connectivity_matches')
    failed=[f"{r['id']}/{r['state']}:{key}" for r in pilot for key in state_fields if not r[key]]
    summary=dict(pilot_molecules=len(PILOT_IDS),pilot_states=len(pilot),pilot_structure_checks=len(pilot)*len(state_fields),
                 failures=failed,registry_entries=len(registry),manifest_states=len(rows),
                 registry_manifest_failures=[r['id'] for r in registry if not all(r[k] for k in ('manifest_present','states_match','smiles_match'))],
                 duplicate_structure_aliases=aliases,
                 pilot_dft_properties_available=sum(r['production_energy_available'] for r in pilot),
                 report_semantics='UMA gas geometry/energy evidence only. No solution properties or usable capacity established.',
                 source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'config/custom_candidates.py',ROOT/'library/manifest.csv')},
                 route_notes_missing_from_scorecard=[r['id'] for r in registry if r['chemistry_notes'] and not r['scorecard_route_note_present']])
    a.out.mkdir(parents=True)
    for name,records in [('registry_manifest.csv',registry),('pilot_states.csv',pilot),('pilot_descriptors.csv',desc),('pilot_readiness.csv',ready),('pilot_gas_diagnostics.csv',gas)]:
        with (a.out/name).open('w',newline='') as fp:
            w=csv.DictWriter(fp,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (a.out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
    if failed or summary['registry_manifest_failures']:raise SystemExit(1)
if __name__=='__main__':main()
