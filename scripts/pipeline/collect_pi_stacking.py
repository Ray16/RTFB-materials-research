#!/usr/bin/env python
"""Collect stationary electronic structures without substituting them for deltaG."""
import csv
import fcntl
import io
import json
import math
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT,EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint

BASE=ROOT/'results/polymer_models/campaigns/2026-10-19_pi_stacking'
OUT=ROOT/'results/polymer_models'


def collect():
    with (OUT/'.collection.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        plan=json.loads((BASE/'plan.json').read_text())
        records={};invalid=[]
        for receipt in (BASE/'done').glob('*.json'):
            try:
                record=verify_state_checkpoint(receipt)
                if not math.isfinite(record['E_smd_eV']):raise ValueError('Nonfinite energy')
                records[record['task']['key']]=record
            except (OSError,KeyError,ValueError) as exc:invalid.append(dict(receipt=str(receipt),error=repr(exc)))
        rows=[];by_model={}
        for task in plan['tasks']:
            ready=BASE/'live_inputs'/task['key']/'task.json'
            if ready.exists():task={**task,**json.loads(ready.read_text())}
            key=task['key'];record=records.get(key);failed=(BASE/'failures'/f'{key}.json').exists()
            row=dict(task,stage='stationary_DFT_electronic_structure' if record else ('failed_preserved' if failed else 'DFT_pending'),
                     E_smd_eV=record['E_smd_eV'] if record else None,
                     deltaG_kJ_mol=None,Hessian_minimum_tested=False,
                     final_xyz=str(Path(record['state_dir'])/'opt.xyz') if record else None,
                     contact_metrics=record['contact_metrics'] if record else None,
                     stacked_pairs=record['stacked_pairs'] if record else None)
            rows.append(row);by_model.setdefault(task['model_key'],{})[task['target']]=row
        folding=[]
        for key,states in by_model.items():
            if 'open' not in states or 'stacked' not in states:continue
            a,b=states['open'],states['stacked']
            ready=a['E_smd_eV'] is not None and b['E_smd_eV'] is not None
            distinct=ready and a['stacked_pairs']==0 and b['stacked_pairs']>0
            folding.append(dict(model_key=key,candidate_id=a['candidate_id'],repeats=a['repeats'],stereoisomer=a['stereoisomer'],
                stacked_minus_open_deltaE_kJ_mol=(b['E_smd_eV']-a['E_smd_eV'])*EV_KJ if distinct else None,
                status='sampled_stationary_conformer_comparison' if distinct else ('end_states_not_distinct_stack_open' if ready else 'pending'),
                reference='Same composition, charge, spin and backbone stereochemistry',
                conformational_population_qualified=False,deltaG_kJ_mol=None))
        associations=[]
        monomers={r['candidate_id']:r for r in rows if r['target']=='monomer'}
        from rdkit import Chem
        from collections import Counter
        def inventory(smiles):return Counter(a.GetAtomicNum() for a in Chem.AddHs(Chem.MolFromSmiles(smiles)).GetAtoms())
        for row in rows:
            if row['model_type']!='two_nonbonded_segments':continue
            mono=monomers[row['candidate_id']]
            assert inventory(row['smiles'])==inventory(mono['smiles'])+inventory(mono['smiles'])
            ready=row['E_smd_eV'] is not None and mono['E_smd_eV'] is not None
            associations.append(dict(model_key=row['model_key'],candidate_id=row['candidate_id'],orientation=row['target'],
                association_deltaE_kJ_mol=(row['E_smd_eV']-2*mono['E_smd_eV'])*EV_KJ if ready else None,
                status='relaxed_electronic_association_screen' if ready else 'pending',
                atom_balanced=True,includes_monomer_deformation=True,BSSE_corrected=False,
                deltaG_kJ_mol=None,association_entropy_qualified=False))
        pilot=BASE/'pilot/status.json'
        status=dict(updated=time.time(),scope='Neutral interaction-method controls; no synthesis admission',
            expected_DFT_cases=len(rows),stationary_electronic_cases=len(records),
            sampled_stack_open_comparisons=sum(f['stacked_minus_open_deltaE_kJ_mol'] is not None for f in folding),
            electronic_association_comparisons=sum(f['association_deltaE_kJ_mol'] is not None for f in associations),
            failures=len(list((BASE/'failures').glob('*.json'))),invalid_receipts=invalid,
            active_claims=[p.name for p in (BASE/'claims').iterdir() if p.is_dir()],
            memory_pilot=json.loads(pilot.read_text())['status'] if pilot.exists() else 'awaiting_free_H100',
            completed_interaction_free_energies=0,charged_state_interactions_computed=False,
            whole_bead_interactions_modeled=False,raw_root=plan['raw_root'],campaign=str(BASE.relative_to(ROOT)),
            viewer='results/polymer_models/index.html',property_scope=plan['property_scope'])
        atomic_write_json(OUT/'interaction_results.json',dict(status=status,models=rows,folding=folding,association=associations))
        atomic_write_json(OUT/'interaction_status.json',status)
        text=io.StringIO();writer=csv.DictWriter(text,fieldnames=['model_key','candidate_id','repeats','stereoisomer','stacked_minus_open_deltaE_kJ_mol','status','deltaG_kJ_mol'],extrasaction='ignore')
        writer.writeheader();writer.writerows(folding)
        tmp=OUT/'interaction_comparisons.csv.tmp';tmp.write_text(text.getvalue());tmp.replace(OUT/'interaction_comparisons.csv')
        subprocess.run([str(ROOT/'.venv-figures/bin/python'),str(ROOT/'scripts/plotting/polymer_models/update_interaction_figures.py')],check=True)
        print(json.dumps({k:status[k] for k in ['stationary_electronic_cases','expected_DFT_cases','sampled_stack_open_comparisons','failures','memory_pilot']}))


if __name__=='__main__':collect()
