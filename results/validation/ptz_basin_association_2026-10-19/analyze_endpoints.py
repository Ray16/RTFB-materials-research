"""Read-only tether geometry diagnostics; angle proximity alone does not prove a basin."""
import os,json,hashlib
from pathlib import Path
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from rdkit import Chem
from rdkit.Chem import Lipinski
from ase.io import read
from redox.core.common import ROOT,read_manifest,smiles_to_xyz_map
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent;gid='hyp_ptz_37difluoro'
row=next(r for r in read_manifest() if r['id']==gid);mol=Chem.MolFromSmiles(row['smiles'])
paths={s:ROOT/'calcs/dft'/gid/s/'opt.xyz' for s in ['neu','ox1']}
angles=[]
for b,c in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
 left=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(b).GetNeighbors() if a.GetIdx()!=c and a.GetAtomicNum()>1)
 right=sorted(a.GetIdx() for a in mol.GetAtomWithIdx(c).GetNeighbors() if a.GetIdx()!=b and a.GetAtomicNum()>1)
 if not left or not right:continue
 ids=[left[0],b,c,right[0]];values={}
 for s,p in paths.items():
  mapping=smiles_to_xyz_map(row['smiles'],p);assert mapping is not None
  xyz_indices=[mapping[i] for i in ids];atoms=read(p);values[s]=atoms.get_dihedral(*xyz_indices)
 diff=(values['ox1']-values['neu']+180)%360-180
 angles.append(dict(smiles_indices=ids,elements=[mol.GetAtomWithIdx(i).GetSymbol() for i in ids],angles_degrees=values,circular_difference_degrees=diff))
report=dict(id=gid,tether_torsions=angles,input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths.values()},basin_association_verified=False,interpretation='Symmetry and charge-dependent core deformation matter; reciprocal optimization is required before assigning the same local conformational branch.')
atomic_write_json(base/'endpoint_geometry.json',report);print(angles)
