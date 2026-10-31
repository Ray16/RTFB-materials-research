import json
from pathlib import Path

import pytest

from redox.core.hashing import sha256_file
from redox.workflows import residency_refusal as recovery


@pytest.fixture
def refused(tmp_path,monkeypatch):
    output=tmp_path/'queue';case=tmp_path/'attempt';archive=tmp_path/'archive';raw=tmp_path/'raw';raw.mkdir()
    source=archive/'scripts/fleet/exclusive_python_worker.py';source.parent.mkdir(parents=True)
    source.write_text('def main():\n    if check():\n        raise RuntimeError('+repr(recovery.MESSAGE)+')\n    runpy.run_path("program")\n')
    spec=dict(prior_root='/original',archive=str(archive),key='key',attempt=str(case),raw_parent=str(raw),
        original_campaign='/original/results/polymer_models/campaigns/test/key',original_raw_parent='/original/calcs/polymer_models/test',original_raw_root='/original/calcs/polymer_models/test/key')
    case.mkdir();(case/'plan.json').write_text('{}');(case/'input.xyz').write_text('geometry')
    namespace=tmp_path/'namespace.json';namespace.write_text(json.dumps(spec));digest=sha256_file(namespace)
    claim=case/'native_dispatch_claim';claim.mkdir();guard=case/'residency_guard';guard.mkdir()
    placement=output/'placements/key/attempt/placement';placement.mkdir(parents=True)
    logical='/original/scripts/fleet/exclusive_python_worker.py'
    owner=dict(host='node',pid=10,process_start='100')
    (claim/'owner.json').write_text(json.dumps(dict(queue=str(output))))
    (claim/'launch.json').write_text(json.dumps(dict(status='controller_spawned_scientific_start_unverified',namespace=str(namespace),namespace_sha256=digest,argv=['launch','--out',str(placement)],owner=owner)))
    (guard/'owner.json').write_text(json.dumps(dict(owner,GPU=0)))
    (guard/'state.json').write_text(json.dumps(dict(GPU=0,foreign_resident_PIDs=[99],own_residency_observed=False)))
    terminal=dict(host='node',pid=11,process_start='101',GPU=0,namespace=str(namespace),namespace_sha256=digest,returncode=1,argv=['--ro-bind',str(archive),'/original',logical])
    (placement/'completion.json').write_text(json.dumps(terminal))
    (placement/'command.log').write_text('Traceback (most recent call last):\n  File "'+logical+'", line 3, in main\nRuntimeError: '+recovery.MESSAGE+'\n')
    monkeypatch.setattr(recovery,'producer_inventory',lambda *args:{'node':[]})
    return spec,output,lambda owner:False,['node']


def test_closed_pre_science_refusal_releases_claim_but_preserves_evidence(refused):
    spec,output,alive,hosts=refused;case=Path(spec['attempt'])
    assert recovery.release_residency_refusal(spec,output,alive,hosts)
    assert not (case/'residency_guard').exists() and not (case/'native_dispatch_claim').exists()
    resolution=next((output/'residency_refusal_history/key').glob('*/resolution.json'));proof=json.loads(resolution.read_text())
    assert not proof['scientific_attempt_consumed'] and not proof['foreign_processes_signalled']
    assert all(sha256_file(p)==h for p,h in proof['source_sha256'].items())


@pytest.mark.parametrize('change',['live','unknown','producer','unreachable','raw','host_guard','wrong_namespace','no_archive_bind','late_branch','no_foreign','own_GPU'])
def test_uncertainty_or_any_scientific_start_preserves_ownership(refused,monkeypatch,change):
    spec,output,alive,hosts=refused;case=Path(spec['attempt']);placement=output/'placements/key/attempt/placement'
    if change=='live':alive=lambda _:True
    if change=='unknown':alive=lambda _:None
    if change=='producer':monkeypatch.setattr(recovery,'producer_inventory',lambda *args:{'node':[{'pid':123}]})
    if change=='unreachable':
        def missing(*args):raise OSError('node unavailable')
        monkeypatch.setattr(recovery,'producer_inventory',missing)
    if change=='raw':(Path(spec['raw_parent'])/'key').mkdir()
    if change=='host_guard':(case/'host_guard').mkdir()
    if change in ('wrong_namespace','no_archive_bind'):
        p=placement/'completion.json';d=json.loads(p.read_text());d['namespace_sha256']='wrong' if change=='wrong_namespace' else d['namespace_sha256']
        if change=='no_archive_bind':d['argv']=[]
        p.write_text(json.dumps(d))
    if change=='late_branch':
        p=Path(spec['archive'])/'scripts/fleet/exclusive_python_worker.py';p.write_text('def main():\n    runpy.run_path("program")\n    if check():\n        raise RuntimeError('+repr(recovery.MESSAGE)+')\n')
    if change in ('no_foreign','own_GPU'):
        p=case/'residency_guard/state.json';d=json.loads(p.read_text())
        if change=='no_foreign':d['foreign_resident_PIDs']=[]
        else:d['own_residency_observed']=True
        p.write_text(json.dumps(d))
    assert not recovery.release_residency_refusal(spec,output,alive,hosts)
    assert (case/'native_dispatch_claim').exists() and (case/'residency_guard').exists()
