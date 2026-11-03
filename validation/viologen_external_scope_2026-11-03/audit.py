"""Recheck exact identity and benchmark overlap before admitting external evidence."""
import hashlib
import json
import os
from pathlib import Path

from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ROOT=Path(__file__).resolve().parents[2]
epoch=ROOT/'manuscript/evidence/2026-11-03_editorial_pass'
benchmark=epoch/'benchmark.json';candidate=epoch/'viologen_native_candidate.json'
doi='10.1002/chem.202200149'
rows=[r for r in json.loads(benchmark.read_text())['outliers'] if r['doi']==doi]
assert len(rows)==2 and len({r['id'] for r in rows})==1
assert {r['wave']:r['exp_V'] for r in rows}=={1:-.785,2:-1.204}
measured=Chem.MolFromSmiles(rows[0]['parent_smiles'])
model=Chem.MolFromSmiles(json.loads(candidate.read_text())['product_smiles'])
assert measured is not None and model is not None
measured_smiles=Chem.MolToSmiles(measured,isomericSmiles=True)
model_smiles=Chem.MolToSmiles(model,isomericSmiles=True)
assert measured_smiles!=model_smiles and Chem.GetFormalCharge(measured)==Chem.GetFormalCharge(model)==2
points=[]
for r in sorted(rows,key=lambda r:r['wave']):
    error=r['calc_V']-r['exp_V'];assert abs(error-r['error_V'])<1e-12
    points.append(dict(wave=r['wave'],experimental_V_vs_Fc=r['exp_V'],
        existing_same_parent_benchmark_calculated_V_vs_Fc=r['calc_V'],
        existing_benchmark_signed_error_V=error,new_exact_graft_error_V=None))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
result=dict(status='closest_experimental_analogue_identity_and_overlap_audited',
    doi=doi,primary_url='https://pmc.ncbi.nlm.nih.gov/articles/PMC9310624/',
    source_location='Results and Discussion, Figure 1 discussion',
    experimental_context=dict(object='soluble benzyl/methyl viologen mediator',
        solvent='acetonitrile',reference='Fc+/0',supporting_electrolyte='tetrabutylammonium hexafluorophosphate',
        supporting_electrolyte_M=.50,mediator_counterions='two hexafluorophosphate ions'),
    measured_identity=dict(smiles=measured_smiles,formula=rdMolDescriptors.CalcMolFormula(measured)),
    local_graft_model_identity=dict(smiles=model_smiles,formula=rdMolDescriptors.CalcMolFormula(model)),
    exact_graph_match=False,matched_solvent_label=True,matched_reference_label=True,
    electrolyte_and_ion_environment_matched=False,
    identity_difference='The local model adds a para-methyl substituent to the benzyl ring',
    already_present_in_retrospective_benchmark=True,overlapping_endpoints=2,
    independent_held_out_endpoints_added=0,points=points,
    exact_graft_experimental_accuracy_established=False,polymer_potential_validated=False,
    existing_benchmark_modified=False,manuscript_numerical_results_changed=False,
    source_sha256={str(p):sha(p) for p in [benchmark,candidate,Path(__file__)]})
out=Path(__file__).with_name('review.json')
if out.exists():raise FileExistsError('Preserve prior scope audits')
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','measured_identity','local_graft_model_identity','independent_held_out_endpoints_added']}))
