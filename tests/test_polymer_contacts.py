"""Physical geometry checks independent of the energy-producing implementation."""
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from redox.build.oligomers import decorate_sites
from redox.build.polymer_contacts import contact_metrics,pi_core_groups,enantiomer_representatives,starting_ensemble

ROOT=Path(__file__).resolve().parents[1]


def test_parallel_geometry_and_rigid_motion_invariance():
    angles=np.arange(6)*np.pi/3
    first=np.column_stack([1.4*np.cos(angles),1.4*np.sin(angles),np.zeros(6)])
    xyz=np.concatenate([first,first+[1.2,0,3.5]])
    groups=[list(range(6)),list(range(6,12))]
    result=contact_metrics(xyz,groups)[0]
    assert result['stack_geometry']
    assert result['plane_separation_A']==pytest.approx(3.5)
    assert result['lateral_offset_A']==pytest.approx(1.2)
    rotation,_=np.linalg.qr(np.random.default_rng(4).normal(size=(3,3)))
    transformed=contact_metrics(xyz@rotation+[10.,-4.,8.],groups)[0]
    for field in ('centroid_distance_A','plane_separation_A','lateral_offset_A','nearest_core_atom_A'):
        assert result[field]==pytest.approx(transformed[field],abs=1e-10)
    assert transformed['stack_geometry']


def test_perpendicular_cores_are_not_parallel_stacks():
    t=np.arange(6)*np.pi/3
    first=np.column_stack([np.cos(t),np.sin(t),np.zeros(6)])
    second=np.column_stack([np.cos(t),np.zeros(6),np.sin(t)])+[0,0,4]
    result=contact_metrics(np.concatenate([first,second]),[list(range(6)),list(range(6,12))])[0]
    assert result['interplane_angle_deg']==pytest.approx(90)
    assert not result['stack_geometry']


@pytest.mark.parametrize('candidate',['design_ndi_ethyl','design_pmdi_ethyl','exp_phthalimide_control','exp_bq_oxy_methoxy'])
def test_pi_cores_belong_to_pendants_not_backbone(candidate):
    feed=json.loads((ROOT/'tests/fixtures/motifs.json').read_text())
    row=next(r for r in feed['original39']+feed['additional_candidates'] if r['id']==candidate)
    variants=decorate_sites(row['frag'],3,(1,2,3))
    mol=Chem.AddHs(Chem.MolFromSmiles(variants[0]));groups=pi_core_groups(mol,row['frag'])
    assert len(groups)==3
    assert not any(set(a)&set(b) for i,a in enumerate(groups) for b in groups[i+1:])
    assert all(len(g)>=6 for g in groups)
    assert len(enantiomer_representatives(variants))==2
    if candidate=='exp_bq_oxy_methoxy':
        assert all(len(g)==6 for g in groups)
        assert all(not mol.GetAtomWithIdx(i).GetIsAromatic() for g in groups for i in g)


def test_trimer_has_a_valid_random_coordinate_fallback():
    fragment='[*:1]N1C(=O)c2ccccc2C1=O'
    smi=decorate_sites(fragment,3,(1,2,3))[0]
    mol,groups,rank=starting_ensemble(smi,fragment,'unbiased',count=1)
    assert len(groups)==3 and len(rank)==1
    assert mol.GetConformer(rank[0][1]).Is3D()
