"""Independently recompute graph, energy, matrix and conformer-selection evidence."""
import hashlib
import json
import math
import re
from pathlib import Path
import time

import numpy as np
from ase import Atoms
from ase.io import read
from ase.build.rotate import minimize_rotation_and_translation
from rdkit import Chem
from rdkit.Chem import Lipinski
from redox.core.common import ROOT,read_manifest
from redox.build.geometry import smiles_to_xyz_map
from redox.validation.protonation import hydrogen_topology_status

base=Path(__file__).resolve().parent
plan=json.loads((base/'plan.json').read_text());result=json.loads((base/'review.json').read_text())
checks=[]
def check(name,passed):checks.append(dict(check=name,passed=bool(passed)))
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for p,h in plan['source_sha256'].items():check('frozen source '+p,sha(p)==h)
canonical=ROOT/'calcs/dft/exp_ndi_hexyl_pentyloxy/red1/opt.xyz'
row=next(r for r in read_manifest() if r['id']==plan['id'] and r['state']=='red1')
mol=Chem.MolFromSmiles(row['smiles']);mapping=smiles_to_xyz_map(row['smiles'],canonical)
torsions=[]
for i,j in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
    a=sorted(n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx()!=j)
    d=sorted(n.GetIdx() for n in mol.GetAtomWithIdx(j).GetNeighbors() if n.GetIdx()!=i)
    if a and d:torsions.append(tuple(mapping[k] for k in (a[0],i,j,d[0])))


def comparison(path,reference):
    actual,ref=read(path),read(reference)
    check('preserved atom order '+str(path),np.array_equal(actual.numbers,ref.numbers))
    mask=actual.numbers>1
    a,b=Atoms(numbers=actual.numbers[mask],positions=actual.positions[mask]),Atoms(numbers=ref.numbers[mask],positions=ref.positions[mask])
    minimize_rotation_and_translation(b,a)
    rmsd=float(np.linalg.norm(a.positions-b.positions)/np.sqrt(len(a)))
    def angle(x,t):
        p,q,r,s=x[list(t)];v=q-p;w=r-q;z=s-r;n1=np.cross(v,w);n2=np.cross(w,z)
        return math.degrees(math.atan2(float(np.dot(np.cross(n1,n2),w/np.linalg.norm(w))),float(np.dot(n1,n2))))
    delta=np.array([abs((angle(actual.positions,t)-angle(ref.positions,t)+180)%360-180) for t in torsions])
    return rmsd,int(np.sum(delta>=60))


eligible=[];records=[]
for record in result['records']:
    key=Path(record['seed']).stem;folder=base/key
    check(key+' successful ALPB',record['returncode']==0 and record['status']=='ALPB_geometry_and_CDS_screen_completed')
    if record['status']!='ALPB_geometry_and_CDS_screen_completed':continue
    path=Path(record['geometry']);text=(folder/'xtb.log').read_text()
    energy=float(re.findall(r'TOTAL ENERGY\s+([-+0-9.EeDd]+)',text)[-1].replace('D','E').replace('d','e'))
    check(key+' ALPB energy/convergence',abs(energy-record['ALPB_energy_Ha'])<1e-12 and 'GEOMETRY OPTIMIZATION CONVERGED' in text)
    check(key+' unchanged declared H graph',hydrogen_topology_status(row['smiles'],path)=='ok')
    check(key+' geometry identity',sha(path)==record['geometry_sha256'])
    coarse=np.load(folder/'coarse_CDS_Ha_Bohr2.npy');fd=np.load(folder/'FD_CDS_Ha_Bohr2.npy')
    check(key+' finite CDS matrices',np.isfinite(coarse).all() and np.isfinite(fd).all())
    asym=float(np.max(abs(coarse-coarse.T)));difference=float(np.max(abs(coarse-fd)))
    check(key+' independent matrix metrics',abs(asym-record['CDS_asymmetry_Ha_Bohr2'])<1e-14 and abs(difference-record['CDS_step_difference_Ha_Bohr2'])<1e-14)
    smooth=asym<=1e-6 and difference<=1e-5;distinct=True
    for name,previous in record['comparisons'].items():
        reference=canonical if name=='canonical' else ROOT/f"results/validation/hexyl_finer_derivatives_2026-10-19/frozen_{name.removeprefix('failed_')}/opt.xyz"
        rmsd,changed=comparison(path,reference)
        check(key+' independent '+name+' conformer metrics',abs(rmsd-previous['heavy_RMSD_A'])<1e-8 and changed==previous['changed_torsions_60deg'])
        distinct &= rmsd>=.5 and changed>=2
    check(key+' predeclared heuristic classification',smooth==record['CDS_smoothness_heuristic_passed'] and distinct==record['distinct_from_preserved_basins'])
    records.append(dict(seed=key,geometry=str(path),ALPB_energy_Ha=energy,CDS_smoothness=smooth,distinct_basin=distinct))
    if smooth and distinct:eligible.append(records[-1])
selected=[]
for r in sorted(eligible,key=lambda r:r['ALPB_energy_Ha']):
    if any((c:=comparison(r['geometry'],s['geometry']))[0]<.5 or c[1]<2 for s in selected):continue
    selected.append(r)
    if len(selected)==2:break
check('independent selection agrees',[r['geometry'] for r in selected]==result['selected_for_future_DFT'])
failed=sum(not c['passed'] for c in checks)
report=dict(status='CPU_seed_preparation_independently_reviewed' if not failed else 'failed',
    cases=len(records),checks=len(checks),failed_checks=failed,details=checks,selected=selected,
    all_records=records,DFT_minima_or_native_partitions_validated=False,
    lowest_DFT_or_free_energy_conformer_selected=False,
    rule='A higher-energy conformer must not replace an unqualified lower minimum merely because it passes derivatives; full native/active-energy comparison and conformer scope remain necessary.',checked=time.time())
(base/'independent_review.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('details','all_records')}));assert not failed
