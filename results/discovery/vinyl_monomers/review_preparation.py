"""Independent graph, mass and source review; no production calculations are touched."""
import csv
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdChemReactions, rdMolDescriptors

BASE=Path(__file__).resolve().parent
proposal=json.loads((BASE/'proposals.json').read_text())
snapshot=Path(proposal['source_snapshot'])
sources={}
source_hashes={}
for name in ('original39.csv','additional_candidates.csv','comparison_controls.csv'):
    path=snapshot/name;source_hashes[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    with path.open() as stream:
        for row in csv.DictReader(stream):sources[row['id']]=row
checks=[]


def check(name, passed):
    checks.append(dict(check=name,passed=bool(passed)))


reaction=rdChemReactions.ReactionFromSmarts('[CH2:1]=[CH:2]-[c:3]>>[CH3:2]-[c:3]')
check('All source designs and the comparison control covered exactly once',
      {r['source_capped_id'] for r in proposal['records']}==set(sources)
      and len(proposal['records'])==len(sources))
for row in proposal['records']:
    source=sources[row['source_capped_id']];mol=Chem.MolFromSmiles(row['smiles'])
    original=Chem.MolFromSmiles(source['exact_grafted_smiles'])
    products=reaction.RunReactants((mol,));restored=set()
    for product in products:
        Chem.SanitizeMol(product[0]);restored.add(Chem.MolToSmiles(product[0],isomericSmiles=True))
    prefix=row['id']+' '
    check(prefix+'independent reaction restores source capped graph',
          restored=={Chem.MolToSmiles(original,isomericSmiles=True)})
    check(prefix+'one additional carbon and unchanged formal charge',
          mol.GetNumAtoms()==original.GetNumAtoms()+1 and Chem.GetFormalCharge(mol)==Chem.GetFormalCharge(original))
    check(prefix+'monoisotopic mass increases by exactly C12',
          abs(rdMolDescriptors.CalcExactMolWt(mol)-rdMolDescriptors.CalcExactMolWt(original)-12.)<1e-7)
    masses=Chem.GetPeriodicTable()
    independent_mass=sum(masses.GetAtomicWeight(atom.GetAtomicNum()) for atom in Chem.AddHs(mol).GetAtoms())
    check(prefix+'independent explicit-atom molecular mass',abs(independent_mass-row['MW_g_mol'])<1e-7)
    states=json.loads(source['states'])
    check(prefix+'declared ladder preserved and not computed',row['states']==states and not row['inherited_computed_properties'])
    electrons=max(state[3] for state in states)-min(state[3] for state in states)
    capacity=electrons*96485.33212331002/(3.6*independent_mass)
    check(prefix+'independent nominal capacity',abs(capacity-row['capacity_nominal_mAh_g'])<1e-8)
    check(prefix+'pending properties remain blank',all(row[k] is None for k in ('E1_vs_Fc_V','E2_vs_Fc_V','lambda_solution_eV')))
    check(prefix+'synthetic and polymer validation remain false',
          not row['vinyl_monomer_synthesis_verified'] and not row['polymerization_verified'] and not row['polymer_redox_performance_verified'])
for path,digest in proposal['source_sha256'].items():
    check('Unchanged producer/library '+path,hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest)
result=dict(status='checked_vinyl_graph_and_descriptor_preparation',proposals=len(proposal['records']),
            checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),
            source_record_sha256=source_hashes,producer_sha256=proposal['source_sha256'],
            quantum_properties_verified=False,polymerization_verified=False,details=checks)
(BASE/'independent_review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('details','source_record_sha256','producer_sha256')}))
assert result['failed_checks']==0
