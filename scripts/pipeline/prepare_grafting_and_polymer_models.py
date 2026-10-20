#!/usr/bin/env python
"""Prepare all local graft reactions and matched oligomer inventories, without QM submission."""
import csv,hashlib,json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from redox.build.grafting import grafting_reaction
from redox.build.oligomers import decorate_sites
from redox.core.common import ROOT,read_result,free_energy
from redox.core.io import atomic_write_json

PILOT={'exp_bq_oxy_methoxy','exp_bq_oxy_trimethyl','exp_aq_1oxy','design_ndi_ethyl','design_pmdi_ethyl','hyp_cata_ndi_ethyl','exp_viologen_methoxyethyl'}

def save_csv(path,rows):
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    feed=json.loads((ROOT/'results/candidates/current.json').read_text());rows=feed['original39']+feed['additional_candidates']
    graft=ROOT/'results/grafting';polymer=ROOT/'results/polymer_models'
    graft.mkdir(exist_ok=True);polymer.mkdir(exist_ok=True)
    assert not (graft/'reactions.json').exists(),'Existing reaction inventory must be reviewed, not overwritten'
    reactions=[];results=[];models=[];species={}
    for row in rows:
        r=grafting_reaction(row['frag'],row['exact_grafted_smiles']);r.update(id=row['id'],name=row['name'],family=row['family'])
        source=row.get('precursor_smiles');r['earlier_precursor_matches']=bool(source and Chem.MolToSmiles(Chem.MolFromSmiles(source))==r['precursor_smiles'])
        key=hashlib.sha256(r['precursor_smiles'].encode()).hexdigest()[:12];r['precursor_id']='gp_'+key
        species[r['precursor_id']]=dict(id=r['precursor_id'],smiles=r['precursor_smiles'],charge=r['precursor_charge'],mult=1,role='grafting_precursor')
        states=json.loads(row['states']) if isinstance(row['states'],str) else row['states'];state=next(s for s in states if s[3]==0)
        gid=row.get('registry_id') or row['id'];record=read_result(gid,state[0]);G,qc=free_energy(record)
        ready=G is not None and qc=='ok' and record.get('thermal_definition')=='native_solution_dft_qrrho_spin'
        r.update(product_id=gid,product_state=state[0],existing_product_G_eV=G if ready else None,product_QC=qc,existing_product_usable=ready)
        reactions.append(r)
        results.append(dict(id=row['id'],name=row['name'],reaction_class=r['reaction_class'],precursor_id=r['precursor_id'],solvent='acetonitrile',temperature_K=298.15,deltaG_reference_eV=None,deltaG_reference_kJ_mol=None,deltaG_with_base_kJ_mol=None,base_correction_required=r['base_correction_required'],status='awaiting_precursor_and_reference_QM' if ready else 'awaiting_precursor_reference_and_product_QM',atom_balanced=True,charge_balanced=True,product_QC=qc,product_state_path=str(ROOT/'calcs/dft'/gid/state[0]),result_record=None))
        for n in (1,2,3):
            for i,smi in enumerate(decorate_sites(row['frag'],n,range(1,n+1))):
                mol=Chem.MolFromSmiles(smi);models.append(dict(id=f"{row['id']}__n{n}__s{i}",candidate_id=row['id'],repeat_units=n,occupied_sites=','.join(map(str,range(1,n+1))),smiles=smi,charge=Chem.GetFormalCharge(mol),formula=rdMolDescriptors.CalcMolFormula(mol),model_type='covalent_H_capped_oligomer',pilot=row['id'] in PILOT,QM_status='not_submitted',deltaG_site_cooperativity_kJ_mol=None))
        if row['id'] in PILOT:
            for occupied in ((1,),(2,)):
                for i,smi in enumerate(decorate_sites(row['frag'],2,occupied)):
                    mol=Chem.MolFromSmiles(smi);models.append(dict(id=f"{row['id']}__n2__sites{occupied[0]}__s{i}",candidate_id=row['id'],repeat_units=2,occupied_sites=str(occupied[0]),smiles=smi,charge=Chem.GetFormalCharge(mol),formula=rdMolDescriptors.CalcMolFormula(mol),model_type='matched_single_site_control',pilot=True,QM_status='not_submitted',deltaG_site_cooperativity_kJ_mol=None))
    for n in (2,3):
        for i,smi in enumerate(decorate_sites(rows[0]['frag'],n,())):
            mol=Chem.MolFromSmiles(smi);models.append(dict(id=f'chlorinated_backbone__n{n}__s{i}',candidate_id='shared_backbone',repeat_units=n,occupied_sites='',smiles=smi,charge=0,formula=rdMolDescriptors.CalcMolFormula(mol),model_type='matched_ungrafted_backbone',pilot=True,QM_status='not_submitted',deltaG_site_cooperativity_kJ_mol=None))
    atomic_write_json(graft/'reactions.json',dict(status='balanced_reference_reactions_prepared',source_snapshot=feed['snapshot_path'],solvent='acetonitrile',reactions=reactions))
    atomic_write_json(graft/'species.json',dict(precursors=list(species.values()),shared_references=[dict(id='merrifield_local_chloride',smiles='Cc1ccc(CCl)cc1',charge=0,mult=1),dict(id='hydrogen_chloride',smiles='Cl',charge=0,mult=1),dict(id='chloride',smiles='[Cl-]',charge=-1,mult=1)]))
    save_csv(graft/'current.csv',results);save_csv(polymer/'current.csv',models)
    ready=sum(r['existing_product_usable'] for r in reactions)
    status=dict(status='prepared_not_submitted',candidate_reactions=len(reactions),balanced_reactions=len(reactions),unique_precursors=len(species),shared_references=3,existing_product_free_energies=ready,products_needing_QM=len(reactions)-ready,completed_grafting_deltaG=0,solvent='acetonitrile',base='not specified; neutral HCl-reference reaction for NH/OH motifs',base_coupled_results_pending=True,new_quantum_jobs_submitted=0)
    atomic_write_json(graft/'status.json',status)
    atomic_write_json(polymer/'status.json',dict(status='graph_inventory_prepared',models=len(models),candidate_count=len(rows),pilot_candidates=sorted(PILOT),covalent_repeats=[1,2,3],whole_bead_interactions_modeled=False,new_quantum_jobs_submitted=0,interaction_cycle='G(P_AA)+G(P_ClCl)-G(P_ACl)-G(P_ClA)',statistical_site_degeneracy_must_be_separated=True))
    print(json.dumps(status));print('polymer graphs',len(models),'earlier precursor mismatches',[r['id'] for r in reactions if not r['earlier_precursor_matches']])


if __name__=='__main__':main()
