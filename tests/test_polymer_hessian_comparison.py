"""A matched protocol and all numerical dimensions must agree separately."""
import copy
import importlib.util
import json
import hashlib
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def panel():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('tested_full_Hessian_comparison',root/'scripts/validation/compare_screening_polymer_hessians.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    limits=json.loads((root/'config/polymer_hessian_reproduction.json').read_text())
    left=dict(namespace='left',plan=dict(original_namespace='same',original_namespace_sha256='sha',
        numerical_settings=dict(solvent='acetonitrile',basis='def2-svp'),symbols=['O','H','H']),
        center=dict(coordinates_Bohr=[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],gradient_Ha_Bohr=np.zeros((3,3)).tolist(),E_Ha=-10.),
        Hessian=np.eye(9,dtype=np.float64),nuclear=dict(frequencies_cm=[30.,40.,50.],nuclear_by_cutoff_eV={'25':1.,'50':2.,'100':3.}))
    right=copy.deepcopy(left);right['namespace']='right'
    return module,left,right,limits


def test_matching_controls_pass_all_numerical_comparisons(panel):
    module,left,right,limits=panel;result=module.pair(left,right,limits)
    assert result['all_within_declared_tolerances']
    assert not any(result['differences'].values())


@pytest.mark.parametrize('axis',['energy','force','Hessian','frequency','nuclear'])
def test_disagreement_is_preserved_as_failed_diagnostic_not_silently_accepted(panel,axis):
    module,left,right,limits=panel
    if axis=='energy':right['center']['E_Ha']+=1e-5
    if axis=='force':right['center']['gradient_Ha_Bohr'][0][0]=1e-4
    if axis=='Hessian':right['Hessian'][0,0]+=.001
    if axis=='frequency':right['nuclear']['frequencies_cm'][0]+=.1
    if axis=='nuclear':right['nuclear']['nuclear_by_cutoff_eV']['50']+=1e-4
    result=module.pair(left,right,limits)
    assert not result['all_within_declared_tolerances'] and not all(result['passed'].values())


@pytest.mark.parametrize('change',['model','source','geometry','nonfinite','matrix_shape','precision'])
def test_incomparable_or_incomplete_controls_are_rejected(panel,change):
    module,left,right,limits=panel
    if change=='model':right['plan']['numerical_settings']['solvent']='water'
    if change=='source':right['plan']['original_namespace']='other'
    if change=='geometry':right['center']['coordinates_Bohr'][0][0]=.01
    if change=='nonfinite':right['Hessian'][0,0]=np.nan
    if change=='matrix_shape':right['Hessian']=np.eye(6)
    if change=='precision':right['Hessian']=right['Hessian'].astype(np.float32)
    with pytest.raises(ValueError):module.pair(left,right,limits)


@pytest.fixture
def records(panel,tmp_path,monkeypatch):
    module,left,right,limits=panel;reviews=[]
    for row in (left,right):
        root=tmp_path/row['namespace'];case=root/'case';(case/'points/0000').mkdir(parents=True)
        namespace=root/'namespace.json';namespace.write_text(json.dumps(dict(case=str(case))))
        row['plan']['screening_variant']='corrected'
        for filename,data in [('plan.json',row['plan']),('points/0000/evaluation.json',row['center']),('nuclear.json',row['nuclear'])]:
            (case/filename).write_text(json.dumps(data))
        np.save(case/'Hessian.npy',row['Hessian'])
        bound=[namespace,case/'plan.json',case/'points/0000/evaluation.json',case/'nuclear.json',case/'Hessian.npy']
        review=root/'review.json';review.write_text(json.dumps(dict(status='corrected_full_polymer_Hessian_independently_reconstructed',
            native_G_qualified=False,namespace=str(namespace),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in bound})))
        reviews.append(review)
    config=tmp_path/'limits.json';config.write_text(json.dumps(limits))
    monkeypatch.setattr(module,'verify',lambda path:json.loads(path.read_text()))
    return module,reviews,config


def test_comparison_binds_actual_numerical_files_and_preserves_source_hashes(records):
    module,reviews,config=records;result=module.compare(reviews,config)
    assert result['controls']==2 and result['all_pairs_within_declared_tolerances']
    assert len(result['source_sha256'])>=14


def test_changed_backing_cannot_hide_behind_a_passing_review_status(records):
    module,reviews,config=records
    np.save(reviews[0].parent/'case/Hessian.npy',np.eye(9)*2)
    with pytest.raises(ValueError,match='backing'):module.compare(reviews,config)


def test_duplicate_review_does_not_count_as_an_independent_reproduction(records):
    module,reviews,config=records
    with pytest.raises(ValueError,match='Distinct'):module.compare([reviews[0],reviews[0]],config)
