"""Recovery must preserve admitted records and refuse ambiguous ownership."""
import importlib.util
import json
from pathlib import Path
import subprocess

import pytest


@pytest.fixture
def recovery(tmp_path, monkeypatch):
    path=Path(__file__).resolve().parents[1]/'results/validation/aryl_solution_2026-10-19/recover_memory.py'
    spec=importlib.util.spec_from_file_location('aryl_memory_recovery_test',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module,'BASE',tmp_path/'campaign')
    monkeypatch.setattr(module,'DFT',tmp_path/'dft')
    monkeypatch.setattr(module,'RECOVERY',tmp_path/'recovery')
    for name in ('done','tasks','failures'):(module.BASE/name).mkdir(parents=True)
    return module


def test_completed_state_cannot_be_archived(recovery,tmp_path):
    task=dict(key='m__ox1',id='m',state='ox1')
    (recovery.BASE/'done/m__ox1.json').write_text('{}')
    with pytest.raises(RuntimeError,match='Completed state'):
        recovery.archive_attempt(task,dict(status='failed_preserved'),tmp_path/'archive')


def test_unreachable_owner_keeps_all_paths(recovery,tmp_path,monkeypatch):
    task=dict(key='m__ox1',id='m',state='ox1')
    original=recovery.BASE/'tasks/m__ox1';original.mkdir();(original/'status.json').write_text('{}')
    def fail(*args,**kwargs):raise subprocess.TimeoutExpired('ssh',20)
    monkeypatch.setattr(recovery.subprocess,'run',fail)
    owner=dict(status='failed_preserved',host='unreachable',pid=100,process_start='1')
    with pytest.raises(subprocess.TimeoutExpired):
        recovery.archive_attempt(task,owner,tmp_path/'archive')
    assert original.is_dir() and not (tmp_path/'archive').exists()


def test_published_geometry_is_not_replayed(recovery,tmp_path,monkeypatch):
    task=dict(key='m__ox1',id='m',state='ox1')
    sd=recovery.DFT/'m/ox1';sd.mkdir(parents=True);(sd/'result.json').write_text('{}')
    monkeypatch.setattr(recovery,'old_owner_absent',lambda *a:None)
    with pytest.raises(RuntimeError,match='Published state result'):
        recovery.archive_attempt(task,dict(status='failed_preserved'),tmp_path/'archive')
    assert (sd/'result.json').read_text()=='{}'


def test_partial_geometries_and_failure_evidence_survive(recovery,tmp_path,monkeypatch):
    task=dict(key='m__ox1',id='m',state='ox1')
    state=recovery.BASE/'tasks/m__ox1';state.mkdir();(state/'status.json').write_text('{"status":"failed_preserved"}')
    failure=recovery.BASE/'failures/m__ox1.json';failure.write_text('{"reason":"OOM"}')
    sd=recovery.DFT/'m/ox1';sd.mkdir(parents=True);(sd/'_opt_00.xyz').write_bytes(b'preserved geometry\n')
    monkeypatch.setattr(recovery,'old_owner_absent',lambda *a:None)
    archive=tmp_path/'archive'
    recovery.archive_attempt(task,dict(status='failed_preserved'),archive)
    assert (archive/'calculation/_opt_00.xyz').read_bytes()==b'preserved geometry\n'
    assert (archive/'failure.json').read_text()=='{"reason":"OOM"}'
    assert len(json.loads((archive/'archive.json').read_text())['source_sha256'])==3
    assert not state.exists() and not sd.exists() and not failure.exists()


def test_live_original_attempt_is_not_terminal(recovery,tmp_path):
    with pytest.raises(RuntimeError,match='not a preserved terminal failure'):
        recovery.archive_attempt(dict(key='m__ox1',id='m',state='ox1'),
            dict(status='computing'),tmp_path/'archive')


def test_native_resume_retains_geometry_energy_and_partial_evidence(recovery,tmp_path,monkeypatch):
    import redox.core.protocol as protocol
    import redox.qm.sp as sp
    import redox.qm.native_thermal as native
    import redox.validation.protonation as protonation
    import redox.validation.readiness as readiness
    task=dict(key='m__ox2',id='m',state='ox2',manifest_row=dict(charge='2',smiles='CC'))
    old=dict(status='failed_preserved',task=task,finished=1)
    state=recovery.BASE/'tasks/m__ox2/status.json';state.parent.mkdir()
    state.write_text(json.dumps(old))
    sd=recovery.DFT/'m/ox2';sd.mkdir(parents=True)
    geometry=b'geometry remains byte identical\n';(sd/'opt.xyz').write_bytes(geometry)
    result=json.dumps(dict(charge=2,mult=1,n_seeds=3,optimized=True,converged_smd=True,e_smd_eV=-100))
    (sd/'result.json').write_text(result)
    address=sd/'native_thermal/test.json';partial=address.with_suffix('');partial.mkdir(parents=True)
    (partial/'producer_source.py').write_text('preserved partial producer')
    monkeypatch.setattr(recovery,'old_owner_absent',lambda *a:None)
    monkeypatch.setattr(protocol,'load_record',lambda *a:dict(e_smd_eV=-100))
    monkeypatch.setattr(sp,'ensure',lambda *a,**k:('present',dict(e_smd_eV=-100)))
    monkeypatch.setattr(native,'cache_path',lambda *a:(address,{}))
    monkeypatch.setattr(native,'select_policy',lambda *a:None)
    monkeypatch.setattr(protonation,'hydrogen_topology_status',lambda *a:'ok')
    monkeypatch.setattr(readiness,'inspect_states',lambda *a:[dict(status='ready')])
    def finish(*a,**k):
        assert not partial.exists()
        address.write_text('{}')
        return 'computed',dict(g_thermal_total_eV=1)
    monkeypatch.setattr(native,'ensure_native_thermal',finish)
    recovery.resume_published_state(task,old,state)
    assert (sd/'opt.xyz').read_bytes()==geometry and (sd/'result.json').read_text()==result
    assert json.loads(state.read_text())['status']=='state_records_verified'
    archived=list((recovery.RECOVERY/'archives/m__ox2').glob('*/incomplete_native_work/producer_source.py'))
    assert len(archived)==1 and archived[0].read_text()=='preserved partial producer'


@pytest.mark.parametrize('changed',['proof','state'])
def test_control_node_release_proof_is_bound_to_owner_and_state(recovery,monkeypatch,changed):
    import hashlib
    recovery.RECOVERY.mkdir()
    owner=dict(host='lambda11',pid=100,process_start='old-start')
    state=recovery.BASE/'tasks/m__ox2/status.json';state.parent.mkdir();state.write_text('{"status":"failed_preserved"}')
    proof=recovery.RECOVERY/'ownership_release.json'
    proof.write_text(json.dumps(dict(records={'m__ox2':dict(owner=owner,
        original_owner_and_predictor_verified_absent=True,
        status_report_sha256=hashlib.sha256(state.read_bytes()).hexdigest())})))
    (recovery.RECOVERY/'recovery_plan.json').write_text(json.dumps(dict(
        ownership_release_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())))
    def unexpected_ssh(*a,**k):raise AssertionError('Sealed control-node release should not require node-to-node SSH')
    monkeypatch.setattr(recovery.subprocess,'run',unexpected_ssh)
    recovery.old_owner_absent(owner,'m__ox2')
    if changed=='proof':proof.write_text(proof.read_text()+' ')
    else:state.write_text('{"status":"computing"}')
    with pytest.raises(AssertionError):recovery.old_owner_absent(owner,'m__ox2')
