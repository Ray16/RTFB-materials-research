"""Reconstruct this study's graphs from CIF bonds; prepare two uncomputed grafts.

The 7-Ph asymmetric unit is exactly half a molecule. Only its deposited 3_656
operation is supported; any changed symmetry is refused instead of guessed.
8-Ph supplies two independent full molecules plus separately retained CH2Cl2.
"""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from ase.io.cif import parse_cif
from rdkit import Chem
from rdkit.Chem import Descriptors,Draw,rdDetermineBonds,rdMolDescriptors
from redox.core.common import ROOT,load_config,read_manifest
from redox.core.io import atomic_write_json
from redox.build.contracts import validate_group
from redox.screening.discovery import precursor_from_frag
from redox.screening.expansion import registry_definitions
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

base=Path(__file__).resolve().parent
raw=ROOT/'data/raw/discovery_sources/2026-10-19/cata_diimides'
checks=[]
def check(label,passed):
 checks.append(dict(check=label,passed=bool(passed)))
 if not passed:raise ValueError(label)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canonical(m):return Chem.MolToSmiles(Chem.RemoveHs(m),isomericSmiles=True)

def reconstruct(name):
 path=raw/f'1860-5397-20-155-{name}.cif'
 b=next(parse_cif(str(path),reader='ase'))
 labels=b['_atom_site_label'];symbols=b['_atom_site_type_symbol'];lookup={label:i for i,label in enumerate(labels)}
 xyz=np.array([b['_atom_site_fract_'+axis] for axis in 'xyz']).T
 cell=np.asarray(b.get_cell());cart=xyz@cell;n=len(labels)
 symmetry=b['_geom_bond_site_symmetry_2'];has_mate=any(s!='.' for s in symmetry)
 if has_mate:
  check(name+' known half-molecule symmetry',set(symmetry)=={'.','3_656'} and b['_space_group_symop_operation_xyz'][2].replace(' ','')=='-x,y,-z+1/2')
  mate=xyz.copy();mate[:,0]=1-xyz[:,0];mate[:,2]=1.5-xyz[:,2]
  cart=np.concatenate((cart,mate@cell));symbols=symbols*2;labels=labels+[label+'_mate' for label in labels]
 rw=Chem.RWMol()
 for label,symbol in zip(labels,symbols):
  a=Chem.Atom(symbol);a.SetNoImplicit(True);a.SetProp('source_label',label);rw.AddAtom(a)
 edges={}
 for a,c,length,sym in zip(b['_geom_bond_atom_site_label_1'],b['_geom_bond_atom_site_label_2'],b['_geom_bond_distance'],symmetry):
  i,j=lookup[a],lookup[c]
  edge=(i,j if sym=='.' else j+n)
  candidates=[edge]
  if has_mate:candidates.append((i+n,j+n if sym=='.' else j))
  for i,j in candidates:
   pair=tuple(sorted((i,j)));check(name+f' bond distance {labels[i]}-{labels[j]}',abs(np.linalg.norm(cart[i]-cart[j])-length)<.003)
   edges[pair]=length
 for i,j in edges:rw.AddBond(i,j,Chem.BondType.SINGLE)
 mol=rw.GetMol();conf=Chem.Conformer(len(symbols))
 for i,p in enumerate(cart):conf.SetAtomPosition(i,p)
 mol.AddConformer(conf)
 fragments=[]
 for indices in Chem.GetMolFrags(mol,sanitizeFrags=False):
  part=Chem.RWMol();mapping={}
  for i in indices:mapping[i]=part.AddAtom(Chem.Atom(mol.GetAtomWithIdx(i)))
  for bond in mol.GetBonds():
   i,j=bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
   if i in mapping and j in mapping:part.AddBond(mapping[i],mapping[j],Chem.BondType.SINGLE)
  frag=part.GetMol();rdDetermineBonds.DetermineBondOrders(frag,charge=0,embedChiral=False)
  Chem.SanitizeMol(frag)
  check(name+' neutral closed-shell component',Chem.GetFormalCharge(frag)==0 and sum(a.GetNumRadicalElectrons() for a in frag.GetAtoms())==0)
  fragments.append(frag)
 return b,fragments


def uncapped_core(parent):
 molecule=Chem.RWMol(Chem.RemoveHs(parent));nitrogens=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==7]
 check('two imide nitrogens',len(nitrogens)==2)
 for i in nitrogens:
  atom=molecule.GetAtomWithIdx(i)
  carbonyl=[a for a in atom.GetNeighbors() if any(b.GetBondType()==Chem.BondType.DOUBLE and b.GetOtherAtom(a).GetAtomicNum()==8 for b in a.GetBonds())]
  check('imide carbonyl connectivity',len(carbonyl)==2)
  carbonyl_ids={a.GetIdx() for a in carbonyl}
  outside=[a for a in atom.GetNeighbors() if a.GetIdx() not in carbonyl_ids]
  check('one N-phenyl cap',len(outside)==1 and outside[0].GetAtomicNum()==6 and outside[0].GetIsAromatic())
  molecule.RemoveBond(i,outside[0].GetIdx());atom.SetNumExplicitHs(1);atom.SetNoImplicit(True)
 fragments=Chem.GetMolFrags(molecule.GetMol(),asMols=True,sanitizeFrags=True)
 core=[m for m in fragments if sum(a.GetAtomicNum()==7 for a in m.GetAtoms())==2]
 check('one intact bis-imide core',len(core)==1)
 return core[0]


def capped(core,caps):
 molecule=Chem.RWMol(core);ns=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()==7]
 for index,cap in zip(ns,caps):
  atom=molecule.GetAtomWithIdx(index);atom.SetNumExplicitHs(0)
  branch=Chem.MolFromSmiles(cap);offset=molecule.GetNumAtoms()
  molecule=Chem.RWMol(Chem.CombineMols(molecule,branch));molecule.AddBond(index,offset,Chem.BondType.SINGLE)
 Chem.SanitizeMol(molecule)
 return molecule.GetMol()


def main():
 with allocation_guard():
  resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
 sources=json.loads((raw/'sources.json').read_text());original=json.loads((ROOT/'data/curated/discovery/cata_diimides_2024.json').read_text())
 for source in sources['files'].values():check('download hash '+source['path'],sha(ROOT/source['path'])==source['sha256'])
 for source,digest in original['source_sha256'].items():check('original curation hash '+source,sha(ROOT/source)==digest)
 parents={};cores={};components={}
 for name,label,expected,count in [('S2','7-Ph','C26H14N2O4',1),('S3','8-Ph','C30H16N2O4',2)]:
  block,fragments=reconstruct(name)
  components[label]=[dict(formula=rdMolDescriptors.CalcMolFormula(f),smiles=canonical(f)) for f in fragments]
  selected=[f for f in fragments if rdMolDescriptors.CalcMolFormula(f)==expected]
  check(label+' formula and component count',len(selected)==count)
  check(label+' independent component graph identity',len({canonical(f) for f in selected})==1)
  parents[label]=Chem.RemoveHs(selected[0]);cores[label]=uncapped_core(parents[label])
  check(label+' carbonyl count',len(parents[label].GetSubstructMatches(Chem.MolFromSmarts('[C,c](=[O])')))==4)
  hex_label=label.replace('Ph','Hex');parents[hex_label]=capped(cores[label],['CCCCCC','CCCCCC'])
  check(hex_label+' NMR total H count',sum(a.GetTotalNumHs() for a in parents[hex_label].GetAtoms())==(30 if label=='7-Ph' else 32))
 rows=[]
 for row in original['records']:
  r=dict(row);mol=parents[r['source_label']]
  r.update(smiles=canonical(mol),formula=rdMolDescriptors.CalcMolFormula(mol),inchikey=Chem.MolToInchiKey(mol),structure_identity_status='CIF bond graph with source symmetry and formula checks; hexyl analogues use the source-declared N-cap replacement',graph_provenance='S2 7-Ph; S3 two independent 8-Ph molecules; SI S4-S6 NMR integrates 30/32 H for 7/8-Hex')
  check(r['source_label']+' no MeCN label transfer',r['solvent']=='CH2Cl2' and not r['MeCN_calibration_admitted'] and not r['exact_graft_property'])
  rows.append(r)
 scaffold=load_config('family_search').SCAFFOLD;known={Chem.MolToInchiKey(validate_group(s,g)) for s,g in registry_definitions()}
 known.update(Chem.MolToInchiKey(Chem.MolFromSmiles(r['smiles'])) for r in read_manifest())
 proposals=[];graft_mols=[]
 for label,core in cores.items():
  variants={canonical(capped(core,caps)) for caps in [['[*:1]','CC'],['CC','[*:1]']]}
  check(label+' symmetry-equivalent N attachment',len(variants)==1)
  fragment=next(iter(variants));topology='cata_ndi' if label=='7-Ph' else 'cata_adi'
  group=dict(id='hyp_'+topology+'_ethyl',name='N-(4-methylbenzyl)-N-ethyl 1,2,3,4-'+('naphthalene' if label=='7-Ph' else 'anthracene')+' diimide',family='imide (n-type)',frag=fragment,states=[['neu',0,1,0],['red1',-1,2,-1],['red2',-2,1,-2]],hypothesis=('Test cata vs conventional NDI as a core-topology isomer pair with the same formula, ethyl cap and N-graft anchor.' if label=='7-Ph' else 'Test adding one fused benzene to the cata-NDI core, with the same ethyl cap and N-graft anchor; core size changes explicitly.'),comparison_baseline=('design_ndi_ethyl' if label=='7-Ph' else 'hyp_cata_ndi_ethyl'),provenance=dict(type='source-directed alternative imide core',publication=original['publication'],source_parent=label,exact_grafted_experiment=False,parent_metrics_inherited=False,production_qm=False,route_status='Source parents are bis-phenyl/bis-hexyl. Mixed N-ethyl/N-graft precursor isolation and selectivity are unverified.'))
  mol=validate_group(scaffold,group);key=Chem.MolToInchiKey(mol)
  check(topology+' new exact graph',key not in known)
  group.update(exact_grafted_smiles=canonical(mol),inchikey=key,formula=rdMolDescriptors.CalcMolFormula(mol),total_atoms_including_H=Chem.AddHs(mol).GetNumAtoms(),MW=Descriptors.MolWt(mol),precursor_smiles=precursor_from_frag(fragment))
  query=Chem.Mol(core)
  for atom in query.GetAtoms():
   if atom.GetAtomicNum()==7:atom.SetNumExplicitHs(0)
  check(topology+' unchanged full core after graft',mol.HasSubstructMatch(query))
  proposals.append(group);graft_mols.append(mol)
 atomic_write_json(base/'source_parents.json',dict(publication=original['publication'],license=original['license'],records=rows,CIF_components=components,source_sha256={v['path']:v['sha256'] for v in sources['files'].values()},parent_measurements_separate_from_grafts=True))
 atomic_write_json(base/'proposals.json',dict(SCAFFOLD=scaffold,GROUPS=proposals))
 Draw.MolsToGridImage([parents[r['source_label']] for r in rows],molsPerRow=2,subImgSize=(440,280),legends=[r['source_label'] for r in rows]).save(str(base/'parents.png'))
 Draw.MolsToGridImage(graft_mols,molsPerRow=2,subImgSize=(480,320),legends=[p['id'] for p in proposals]).save(str(base/'prepared_grafts.png'))
 atomic_write_json(base/'summary.json',dict(status='source_graphs_curated_and_two_exact_grafts_prepared',source_parents=4,prepared_exact_grafts=2,checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),registered_candidates=0,quantum_jobs_submitted=0,resources=resources,CPU_affinity=sorted(os.sched_getaffinity(0)),model_properties_computed=False,source_adapter_sha256=sha(Path(__file__)),details=checks,created=time.time()))
 print(json.dumps(dict(checks=len(checks),prepared=[dict(id=p['id'],atoms=p['total_atoms_including_H'],smiles=p['exact_grafted_smiles']) for p in proposals])))
if __name__=='__main__':main()
