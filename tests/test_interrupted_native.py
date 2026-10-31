import json
from pathlib import Path

import pytest

from redox.core.hashing import sha256_file
from redox.workflows import interrupted_native as continuation


@pytest.fixture
def source(tmp_path,monkeypatch):
    case=tmp_path/'original_case';case.mkdir();raw=tmp_path/'original_raw';raw.mkdir();archive=tmp_path/'archive';archive.mkdir()
    (case/'input.xyz').write_text('exact geometry')
    plan=dict(tasks=[dict(charge=0,mult=1)],producer_sha256={},input_sha256={},maximum_host_RSS_GiB=80,
        maximum_peak_GPU_MiB=31744,experimental_model={'solvent':'acetonitrile'},
        execution_driver='/prior/scripts/pipeline/smooth_polymer_native.py',scored_geometry_sha256=sha256_file(case/'input.xyz'))
    (case/'plan.json').write_text(json.dumps(plan))
    owner=dict(host='node',pid=42,process_start='100');guard=case/'residency_guard';guard.mkdir();native=case/'experimental_native';native.mkdir()
    (native/'status.json').write_text(json.dumps(dict(status='experimental_native_failed_preserved',finished=1,owner=owner,error='KeyboardInterrupt from own guard')))
    (guard/'owner.json').write_text(json.dumps(dict(owner,guard_stopped=True)))
    (guard/'interruption.json').write_text(json.dumps(dict(pid=42,foreign_resident_PIDs=[99],foreign_processes_signaled=False,signal_target='own reserved process only')))
    spec=dict(attempt=str(case),raw_parent=str(raw),archive=str(archive),prior_root='/prior',key='molecule',maximum_host_RSS_GiB=80)
    path=tmp_path/'namespace.json';path.write_text(json.dumps(spec))
    monkeypatch.setattr(continuation,'verify_spec',lambda path:json.loads(path.read_text()))
    class Reader:
        def __init__(self,*args):pass
        def verify_all(self,*args):pass
        def recheck(self):pass
    monkeypatch.setattr(continuation,'NamespaceSources',Reader)
    monkeypatch.setattr(continuation,'producer_inventory',lambda *args:{'node':[]})
    return path,tmp_path/'new',lambda owner:False,['node']


def test_foreign_interruption_preserves_geometry_and_every_old_partial(source):
    path,out,alive,hosts=source;old=json.loads(path.read_text());oldbytes=(Path(old['attempt'])/'plan.json').read_bytes()
    spec=continuation.prepare(path,out,alive,hosts)
    assert (Path(spec['attempt'])/'plan.json').read_bytes()==oldbytes
    assert (Path(old['attempt'])/'plan.json').read_bytes()==oldbytes
    assert not spec['scientific_settings_changed'] and spec['original_partial_outputs_preserved']


@pytest.mark.parametrize('change',['live','unknown','producer','no_foreign','wrong_medium'])
def test_live_unknown_or_other_failure_cannot_be_restarted_as_residency_failure(source,monkeypatch,change):
    path,out,alive,hosts=source;spec=json.loads(path.read_text());case=Path(spec['attempt'])
    if change=='live':alive=lambda _:True
    if change=='unknown':alive=lambda _:None
    if change=='producer':monkeypatch.setattr(continuation,'producer_inventory',lambda *args:{'node':[{'pid':1}]})
    if change=='no_foreign':
        p=case/'residency_guard/interruption.json';d=json.loads(p.read_text());d['foreign_resident_PIDs']=[];p.write_text(json.dumps(d))
    if change=='wrong_medium':
        p=case/'plan.json';d=json.loads(p.read_text());d['experimental_model']['solvent']='water';p.write_text(json.dumps(d))
    with pytest.raises(ValueError):continuation.prepare(path,out,alive,hosts)
    assert not out.exists()
