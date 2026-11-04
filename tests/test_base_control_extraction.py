import hashlib
import importlib.util
import json
from pathlib import Path
import pytest


def preparer():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('base_extraction_test',root/'scripts/validation/prepare_base_reaction_controls.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def inputs(root):
    root.mkdir();geometry=root/'geometry.xyz';geometry.write_text('1\nTest-only data fixture\nCl 0 0 0\n')
    energy=root/'energy.json';energy.write_text('{}')
    states=[dict(id=k,aliases=a,smiles='Cl',charge=0,mult=1,files=dict(geometry=geometry.name,energy=energy.name),original_state_directory='fixture')
            for k,a in [('a',[]),('b',[]),('fc_n',['benchmark:ferrocene:neu']),('fc_o',['benchmark:ferrocene:ox'])]]
    reactions=[]
    for i in range(5):
        reactions.append(dict(name='reference:c'+str(i),coefficients={'a':-1,'b':1},delta_n_solutes=0))
        for base in ['dbu','dipea']:
            for model in ['free_ions_kJ_mol','contact_pair_kJ_mol']:
                reactions.append(dict(name=base+':'+model+':c'+str(i),coefficients={'a':-1,'b':1},delta_n_solutes=0))
    data=dict(states=states,reactions=reactions,files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [geometry,energy]})
    (root/'manifest.json').write_text(json.dumps(data));return data


@pytest.mark.parametrize('bad',['scope','state_id','external_input'])
def test_extraction_rejects_incomplete_scope_or_path_escape(tmp_path,bad):
    root=tmp_path/'source';data=inputs(root)
    if bad=='scope':data['reactions'].pop()
    if bad=='state_id':
        data['states'][0]['id']='../a'
        for r in data['reactions']:r['coefficients']={'../a':-1,'b':1}
    if bad=='external_input':data['states'][0]['files']['geometry']='../outside.xyz'
    (root/'manifest.json').write_text(json.dumps(data))
    with pytest.raises(ValueError):preparer().prepare(root,tmp_path/'output')
