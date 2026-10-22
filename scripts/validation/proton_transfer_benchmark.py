#!/usr/bin/env python
"""Check relative base energetics without HCl, chloride or a free proton reference."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
CAMPAIGN = ROOT/'results/validation/proton_transfer_2026-10-22'
OUT = CAMPAIGN
BASE_PILOT = ROOT/'results/grafting/campaigns/2026-10-19_base_assisted_pilot'
TABLE = ROOT/'calcs/validation/acid_base_external_2026-10-22/sources/basicity.html'
BASES = {
    'pyridine':dict(neutral='c1ccncc1',protonated='c1cc[nH+]cc1'),
    'triethylamine':dict(neutral='CCN(CC)CC',protonated='CC[NH+](CC)CC'),
}


def engine():
    module=api.load_module('proton_transfer_component_runner',ROOT/'scripts/pipeline/base_assisted_grafting.py')
    module.BASE=CAMPAIGN;module.OUT=OUT;module.SCRIPT=SCRIPT
    module.G.BASE=CAMPAIGN;module.G.SCRIPT=SCRIPT
    original=module.G.run_state
    checkpoint=api.load_module('proton_transfer_checkpoints',ROOT/'scripts/validation/complete_benchmark_spin.py')

    def evaluate(task,plan):
        checkpoint.OUT=CAMPAIGN/'tasks'/task['key']
        with checkpoint.checkpoint_optimizers():
            return original(task,plan)

    module.G.run_state=evaluate
    return module


def launch_with_health_checks():
    """An unavailable device must not prevent other suitable devices from starting."""
    import socket,time
    from redox.workflows import resources
    original=resources.gpu_suitability

    def checked(gpu):
        try:return original(gpu)
        except (subprocess.TimeoutExpired,subprocess.CalledProcessError) as error:
            with (CAMPAIGN/('resource_probe_'+socket.gethostname()+'.jsonl')).open('a') as stream:
                stream.write(json.dumps(dict(gpu=gpu,error=repr(error),checked=time.time(),
                    interpretation='Unavailable for this launch; no molecular attempt consumed'))+'\n')
            raise RuntimeError('GPU health/identity probe unavailable') from error

    resources.gpu_suitability=checked
    try:engine().launch()
    finally:resources.gpu_suitability=original


def balanced_transfer(neutral,protonated,reference_neutral,reference_protonated):
    from rdkit import Chem
    from redox.build.grafting import atom_inventory
    molecules=[Chem.MolFromSmiles(s) for s in [neutral,protonated,reference_neutral,reference_protonated]]
    if any(m is None for m in molecules):raise ValueError('Invalid proton-transfer graph')
    a,b,c,d=molecules
    return (atom_inventory(a)+atom_inventory(d)==atom_inventory(b)+atom_inventory(c)
            and Chem.GetFormalCharge(a)+Chem.GetFormalCharge(d)==Chem.GetFormalCharge(b)+Chem.GetFormalCharge(c))


def prepare():
    from rdkit import Chem
    from redox.build.structures import conformer_ensemble,to_xyz
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os,time
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    parser=api.load_module('proton_transfer_author_table',ROOT/'scripts/validation/benchmark_acid_base.py')
    table=parser.BasicityTable();table.feed(TABLE.read_text())
    reference_smiles=dict(neutral='C1CCC2=NCCCN2CC1',protonated='C1CCC2=[NH+]CCCN2CC1')
    anchor=parser.exact_base(table,reference_smiles['neutral'])
    if anchor is None:raise ValueError('Exact reference base missing from author table')
    references={};reference_sources={}
    for form in ('neutral','protonated'):
        receipt=BASE_PILOT/'done'/('dbu_'+form+'.json');report=verify_state_checkpoint(receipt)
        if report['status']!='state_records_verified':raise ValueError('DBU reference does not qualify')
        references[form]=dict(G_eV=report['G_eV'],receipt=str(receipt),receipt_sha256=api.digest(receipt))
        reference_sources.update(report['output_sha256']);reference_sources[str(receipt)]=api.digest(receipt)
    CAMPAIGN.mkdir()
    for name in ('claims','done','failures','tasks','workers','workspace/library','workspace/calcs'):
        (CAMPAIGN/name).mkdir(parents=True,exist_ok=True)
    raw=ROOT/'calcs/validation'/CAMPAIGN.name;raw.mkdir()
    (CAMPAIGN/'workspace/calcs/dft').symlink_to(raw,target_is_directory=True)
    tasks=[];comparisons=[]
    for name,forms in BASES.items():
        source=parser.exact_base(table,forms['neutral'])
        if source is None:raise ValueError('No exact external base anchor: '+name)
        if not balanced_transfer(forms['neutral'],forms['protonated'],reference_smiles['neutral'],reference_smiles['protonated']):
            raise ValueError('Proton-transfer atom/charge balance failed')
        comparisons.append(dict(base=name,source=source,reference=anchor,
            reaction='B + DBUH+ -> BH+ + DBU',delta_n_solutes=0,
            external_deltaG_kJ_mol=parser.capture_deltaG(anchor['pKaH_MeCN'],source['pKaH_MeCN']),
            source_condition_qualification='Author-revised MeCN relative scale; exact source temperature/activity review remains open'))
        for form,smiles in forms.items():
            key=name+'_'+form;mol=Chem.MolFromSmiles(smiles);charge=Chem.GetFormalCharge(mol)
            directory=CAMPAIGN/'workspace/library'/key;conf=directory/'conformers';conf.mkdir(parents=True)
            molecule,cids,forcefield,energies=conformer_ensemble(mol,net_charge=charge,top_k=3)
            for index,(cid,energy) in enumerate(zip(cids,energies)):
                one=Chem.Mol(molecule);one.RemoveAllConformers();one.AddConformer(molecule.GetConformer(cid),assignId=True)
                (conf/f'conf_{index:02d}.xyz').write_text(to_xyz(one,f'{forcefield} E={energy}; proton-transfer method control'))
            shutil.copyfile(conf/'conf_00.xyz',directory/'rest.xyz')
            row=dict(id=key,name=key,family='proton_transfer_method_control',state='rest',charge=charge,
                     mult_hint=1,n_e=0,n_conf=len(cids),solv_preopt=1,smiles=smiles)
            tasks.append(dict(key=key,id=key,state='rest',role='base_method_control_not_synthesis_admission',
                atoms=Chem.AddHs(mol).GetNumAtoms(),manifest_row=row,atom=False,preopt='alpb',
                input_sha256={str(p):api.digest(p) for p in directory.rglob('*') if p.is_file()}))
    snapshot=CAMPAIGN/'reference_snapshot.json'
    atomic_write_json(snapshot,dict(references=references,reference_sources=reference_sources,comparisons=comparisons,
        table_sha256=api.digest(TABLE),table_url='https://analytical.chem.ut.ee/databases/pka-values/'))
    sources=[SCRIPT,ROOT/'scripts/pipeline/grafting_campaign.py',ROOT/'scripts/pipeline/base_assisted_grafting.py',
        ROOT/'scripts/validation/complete_benchmark_spin.py',ROOT/'scripts/validation/benchmark_acid_base.py']
    sources += list((ROOT/'src/redox/qm').glob('*.py'))+list((ROOT/'src/redox/core').glob('*.py'))
    sources += [ROOT/'src/redox'/p for p in ['build/structures.py','build/geometry.py','validation/protonation.py',
        'validation/readiness.py','validation/spin_ensemble.py','workflows/geometry_checkpoints.py']]
    sources += [ROOT/'config/project.json',ROOT/'config/electrolyte.py',*list((ROOT/'config').glob('native_*.json'))]
    original=api.read_plan(BASE_PILOT)
    atomic_write_json(CAMPAIGN/'plan.json',dict(schema=1,created=time.time(),tasks=tasks,
        raw_root=str(raw),attempts_per_item=1,worker_consecutive_failure_limit=3,
        optimizer_stopping_tolerances=original['optimizer_stopping_tolerances'],
        producer_sha256={str(p):api.digest(p) for p in sorted(set(sources))},
        input_sha256={str(TABLE):api.digest(TABLE),str(snapshot):api.digest(snapshot),**reference_sources},
        purpose='Exact-base proton-transfer method check; HCl/Cl and isolated-proton terms cancel or are absent',
        temperature_K=298.15,solvent='acetonitrile',full_grafting_accuracy_validated=False))
    collect()


def collect():
    plan=api.read_plan(CAMPAIGN);api.verify_sources(plan['input_sha256'])
    source=json.loads((CAMPAIGN/'reference_snapshot.json').read_text());values={}
    for path in (CAMPAIGN/'done').glob('*.json'):
        report=verify_state_checkpoint(path);values[report['task']['key']]=report['G_eV']
    rows=[]
    for comparison in source['comparisons']:
        name=comparison['base'];row=dict(comparison,DFT_deltaG_kJ_mol=None,DFT_minus_external_kJ_mol=None)
        if name+'_neutral' in values and name+'_protonated' in values:
            delta=(values[name+'_protonated']+source['references']['neutral']['G_eV']
                   -values[name+'_neutral']-source['references']['protonated']['G_eV'])*EV_KJ
            row.update(DFT_deltaG_kJ_mol=delta,DFT_minus_external_kJ_mol=delta-row['external_deltaG_kJ_mol'])
        rows.append(row)
    atomic_write_json(OUT/'current.json',dict(status='components_complete' if len(values)==4 else 'components_pending',
        components_qualified=len(values),components_expected=4,comparisons=rows,source_plan_sha256=api.digest(CAMPAIGN/'plan.json'),
        raw_root=plan['raw_root'],temperature_K=298.15,solvent='acetonitrile',
        interpretation='Relative basicity check only; not calibrated full grafting, ion pairing, activation barriers or yield',
        source_conditions_fully_qualified=False,full_grafting_accuracy_validated=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=('prepare','worker','launch','collect'))
    args=parser.parse_args()
    if args.action=='prepare':prepare()
    elif args.action=='collect':collect()
    elif args.action=='launch':launch_with_health_checks()
    else:engine().worker()
