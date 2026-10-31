import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from redox.core.protocol import ACTIVE_SP,geom_sha1,record_path
from redox.workflows.molecular_recovery import recovered_metadata,select_verified_repair


@pytest.fixture
def closed(tmp_path):
    case=tmp_path/'case';case.mkdir()
    (case/'plan.json').write_text(json.dumps(dict(charge=0,mult=1,product_smiles='CCO')))
    (case/'status.json').write_text(json.dumps(dict(status='molecular_native_failed_preserved',
        error="OSError(18, 'Invalid cross-device link')",finished=1)))
    (case/'residency_guard').mkdir();(case/'residency_guard/owner.json').write_text(json.dumps(dict(guard_stopped=True)))
    base=case/'repairs/id/neu/run';base.mkdir(parents=True);branches=[]
    for sign,energy in [(-1,-10.),(1,-10.00001)]:
        directory=base/f'branch_{sign}';directory.mkdir();mol=Chem.AddHs(Chem.MolFromSmiles('CCO'))
        AllChem.EmbedMolecule(mol,randomSeed=17);Chem.MolToXYZFile(mol,str(directory/'opt.xyz'))
        path=record_path(directory,0,1);path.parent.mkdir()
        path.write_text(json.dumps(dict(protocol=ACTIVE_SP,protocol_hash=path.stem,geom_sha1=geom_sha1(directory/'opt.xyz'),
            charge=0,mult=1,converged_smd=True,converged_gas=True,e_smd_eV=energy,e_gas_eV=energy+1.,source='computed:unit_fixture')))
        branches.append(dict(sign=sign,status='verified',geometry_optimization_converged=True,optimization_SCF_converged=True,
            geometry='/old/results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs/id/neu/run/'+directory.name+'/opt.xyz',e_smd_eV=energy))
    status=base/'status.json';status.write_text(json.dumps(dict(branches=branches)))
    return dict(case=str(case),prior_root='/old',source_id='id',state='neu'),status


def test_both_terminal_branches_select_the_same_lowest_active_energy_rule(closed):
    spec,status=closed
    assert select_verified_repair(spec)[3].name=='branch_1'


@pytest.mark.parametrize('change',['pending','unmatched_energy','other_failure','live_guard'])
def test_incomplete_or_different_failure_cannot_be_recovered(closed,change):
    spec,status=closed;case=Path(spec['case'])
    data=json.loads(status.read_text())
    if change=='pending':data['branches'][1]['status']='optimizing'
    if change=='unmatched_energy':data['branches'][1]['e_smd_eV']=-11.
    status.write_text(json.dumps(data))
    if change=='other_failure':(case/'status.json').write_text(json.dumps(dict(status='molecular_native_failed_preserved',error='SCF failed',finished=1)))
    if change=='live_guard':(case/'residency_guard/owner.json').write_text(json.dumps(dict(guard_stopped=False)))
    with pytest.raises(ValueError):select_verified_repair(spec)


def test_source_metadata_relocates_provenance_without_fabricating_energy_proof(tmp_path):
    source=tmp_path/'source.json';source.write_text('{}')
    raw=dict(charge=0,mult=1,energy_protocol={'stale':True},energy_geom_sha1='old',e_smd_eV=-10.)
    original=tmp_path/'original';target=tmp_path/'target'
    result=recovered_metadata(raw,original,source,target)
    assert result['e_smd_eV']==-10. and result['charge']==0 and result['mult']==1
    assert 'energy_protocol' not in result and 'energy_geom_sha1' not in result
    assert result['minimum_repair_source']==str(original/'opt.xyz')
    assert result['opt_out']==str(target/'opt.xyz')
