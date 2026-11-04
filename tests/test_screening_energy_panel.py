"""Reject missing endpoints and changed source geometry/protocol in the panel."""
import copy
import hashlib
import json

import pytest

from redox.core.common import HARTREE_EV
from redox.core.protocol import ACTIVE_SP
from redox.workflows.screening_energy_panel import benchmark_inventory,output_budget


@pytest.fixture
def panel(tmp_path):
    xyz=tmp_path/'geometry.xyz';xyz.write_text('2\nfixture\nH 0 0 0\nH 0 0 0.75\n')
    data=dict(states=[],files_sha256={})
    for key,aliases,charge,mult in [('a',['benchmark:ferrocene:ox'],1,2),
            ('b',['benchmark:ferrocene:neu','benchmark:toy:neu'],0,1),
            ('c',['benchmark:toy:red1'],-1,2)]:
        energy=tmp_path/(key+'.json')
        energy.write_text(json.dumps(dict(protocol=ACTIVE_SP,charge=charge,mult=mult,converged_smd=True,
            geom_sha1=hashlib.sha1(xyz.read_bytes()).hexdigest()[:16],e_smd_Ha=-1.,e_smd_eV=-HARTREE_EV)))
        files=dict(geometry=xyz.name,energy=energy.name)
        data['states'].append(dict(id=key,aliases=aliases,smiles='[H][H]',charge=charge,mult=mult,
            original_state_directory='/source/'+key,files=files))
        for path in [xyz,energy]:data['files_sha256'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    benchmark=dict(admitted_events=1,outliers=[dict(id='toy',event='neu->red1')])
    return benchmark,data,tmp_path


def test_all_required_states_and_reference_are_kept(panel):
    benchmark,data,root=panel
    rows=benchmark_inventory(benchmark,data,root)
    assert {r['key'] for r in rows}=={'a','b','c'}
    assert all(r['reference_energy_Ha']==-1. for r in rows)
    assert all(r['output_budget_bytes']>0 for r in rows)


@pytest.mark.parametrize('change',['missing_alias','duplicate_alias','geometry','changed_protocol','stale_geometry_parent'])
def test_incomplete_or_changed_backing_is_rejected(panel,change):
    benchmark,data,root=panel
    data=copy.deepcopy(data)
    if change=='missing_alias':data['states'].pop()
    if change=='duplicate_alias':data['states'][0]['aliases'].append('benchmark:toy:neu')
    if change=='geometry':(root/'geometry.xyz').write_text('2\nchanged\nH 0 0 0\nH 0 0 1.5\n')
    if change in ('changed_protocol','stale_geometry_parent'):
        p=root/'b.json';e=json.loads(p.read_text())
        if change=='changed_protocol':e['protocol']['basis']='def2-svp'
        else:e['geom_sha1']='0000000000000000'
        p.write_text(json.dumps(e));data['files_sha256'][p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    with pytest.raises(ValueError):benchmark_inventory(benchmark,data,root)


@pytest.mark.parametrize('args',[(True,1,False),(0,1,False),(10,0,False),(10,2,'Fe')])
def test_unbounded_or_ambiguous_checkpoint_dimensions_are_rejected(args):
    with pytest.raises(ValueError):output_budget(*args)
