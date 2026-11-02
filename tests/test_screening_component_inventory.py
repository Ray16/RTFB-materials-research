from copy import deepcopy
import hashlib
import json
from pathlib import Path
import shutil
import pytest
pytest.importorskip('pyscf')
from redox.workflows.screening_component_panel import inventory


def package(tmp_path):
    root=Path(__file__).resolve().parents[1]/'examples/exact_viologen_replay/data'
    data=json.loads((root/'manifest.json').read_text());row=next(r for r in data['states'] if r['id']=='chloride')
    files={}
    for name in ('geometry','energy'):
        target=tmp_path/(name+Path(row['files'][name]).suffix)
        shutil.copyfile(root/row['files'][name],target);files[name]=target.name
    row=dict(row,files=files)
    manifest=dict(states=[row],files_sha256={f:hashlib.sha256((tmp_path/f).read_bytes()).hexdigest() for f in files.values()})
    (tmp_path/'manifest.json').write_text(json.dumps(manifest));return manifest


def test_bound_monatomic_component_is_explicit(tmp_path):
    package(tmp_path);_,rows=inventory(tmp_path)
    assert len(rows)==1 and rows[0]['atoms']==1 and rows[0]['source_state']['charge']==-1
    assert rows[0]['source_state']['aliases']==['exact_viologen:chloride']


@pytest.mark.parametrize('bad',['file_hash','charge','basis','nonfinite_eV'])
def test_rehashed_wrong_state_or_nonfinite_unit_cannot_enter_panel(tmp_path,bad):
    data=package(tmp_path)
    if bad=='file_hash':data['files_sha256']['energy.json']='0'*64
    elif bad=='charge':data['states'][0]['charge']=0
    else:
        p=tmp_path/'energy.json';record=json.loads(p.read_text())
        if bad=='basis':record['protocol']['basis']='def2-svp'
        else:record['e_smd_eV']=float('nan')
        p.write_text(json.dumps(record));data['files_sha256']['energy.json']=hashlib.sha256(p.read_bytes()).hexdigest()
    (tmp_path/'manifest.json').write_text(json.dumps(data))
    with pytest.raises(ValueError):inventory(tmp_path)
