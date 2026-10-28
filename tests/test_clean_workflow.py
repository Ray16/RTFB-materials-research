"""Regression contracts for the clean CLI, exports and bounded parallel scheduler."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import threading

import pytest

from redox.core.hashing import sha256_file
from redox.workflows.environment import THREAD_KEYS, worker_environment
from redox.workflows.stages import execute_graph, validate_graph


@dataclass(frozen=True)
class Stage:
    name: str
    dependencies: tuple = ()


def test_parallel_stages_overlap_but_join_waits_for_both_dependencies():
    barrier = threading.Barrier(2, timeout=5)
    completed = set()
    stages = [Stage('a'), Stage('b'), Stage('join', ('a', 'b'))]
    def run(stage, slot):
        if stage.name in ('a', 'b'):
            barrier.wait()
        else:
            assert completed == {'a', 'b'}
        return 0
    result = execute_graph(stages, run, workers=2,
                           finished=lambda stage, code: completed.add(stage.name))
    assert result == dict(completed=['a', 'b', 'join'], failed=[], not_started=[])


def test_failure_keeps_running_results_but_never_starts_dependents():
    barrier = threading.Barrier(2, timeout=5)
    called = []
    stages = [Stage('failure'), Stage('peer'), Stage('dependent', ('failure',)), Stage('later', ('peer',))]
    def run(stage, slot):
        called.append(stage.name)
        if stage.name in ('failure', 'peer'):
            barrier.wait()
        if stage.name == 'peer':
            # Make the failure observable before the other branch becomes ready.
            failed.wait(5)
        return int(stage.name == 'failure')
    failed = threading.Event()
    def finished(stage, code):
        if code: failed.set()
    result = execute_graph(stages, run, workers=2, finished=finished)
    assert set(called) == {'failure', 'peer'}
    assert result['failed'] == ['failure']
    assert result['not_started'] == ['dependent', 'later']


@pytest.mark.parametrize('stages', [[Stage('x'), Stage('x')], [Stage('x', ('missing',))],
                                    [Stage('x', ('y',)), Stage('y', ('x',))]])
def test_invalid_graph_fails_before_any_work(stages):
    with pytest.raises(ValueError):
        execute_graph(stages, lambda *args: pytest.fail('unexpected work'), workers=2)


def test_worker_environment_isolated_and_thread_caps_complete(tmp_path):
    base = {'PYTHONPATH': '/old/src', 'OMP_NUM_THREADS': '64', 'PATH': '/usr/bin'}
    env = worker_environment(tmp_path, threads=2, base=base, python='/new/bin/python')
    assert env['PYTHONPATH'] == str(tmp_path / 'src')
    assert all(env[key] == '2' for key in THREAD_KEYS)
    assert env['PATH'] == '/new/bin:/usr/bin'
    assert base['OMP_NUM_THREADS'] == '64'


def test_streaming_hash_matches_sha256_across_chunk_boundaries(tmp_path):
    data = bytes(range(256)) * 1031
    path = tmp_path / 'hessian.npy'; path.write_bytes(data)
    assert sha256_file(path, chunk_size=997) == hashlib.sha256(data).hexdigest()
    with pytest.raises(ValueError): sha256_file(path, chunk_size=0)


def test_outer_sphere_v2_priority_and_legacy_metadata(tmp_path):
    from redox.reports import outer_sphere
    header = 'id,couple,born_lambda_o_eV,pcm_lambda_o_eV,pcm_over_born,sasa_radius_A,solvent,xc,basis,scf_converged\n'
    v2 = tmp_path / 'v2.csv'; v1 = tmp_path / 'v1.csv'
    v2.write_text(header + 'm_ox_red,ox/red,1,2,2,3,acetonitrile,wb97m-v,def2-tzvpd,True\n')
    v1.write_text(header + 'm_ox_red,ox/red,9,9,1,3,acetonitrile,,,\n' + 'n_ox_red,ox/red,4,5,1.25,3,acetonitrile,,,\n')
    out = tmp_path / 'new/results.csv'
    outer_sphere.main(sources=[v2, v1], output=out)
    import csv
    rows = list(csv.DictReader(out.open()))
    assert [(r['id'], r['lambda_o_pcm_eV']) for r in rows] == [('m', '2'), ('n', '5')]
    assert rows[0]['scf_checked'] == 'True'
    assert rows[1]['method'] == 'b3lyp/6-31g(d,p) (legacy)'


def test_qualification_boolean_alone_cannot_admit_nitro(tmp_path, monkeypatch):
    from redox.qm import uma_policy
    monkeypatch.setattr(uma_policy, 'ROOT', tmp_path)
    (tmp_path / 'gate.json').write_text(json.dumps({'production_qualified': True}))
    with pytest.raises(RuntimeError, match='Incomplete'):
        uma_policy.require_production_qualification({'qualification_file': 'gate.json',
            'checkpoint_sha256': 'a', 'reference_sha256': {}})


def test_read_manifest_cache_never_serves_stale_rows(tmp_path, monkeypatch):
    from redox.core import common
    library = tmp_path / 'library'; library.mkdir()
    monkeypatch.setattr(common, 'LIBRARY', library)
    manifest = library / 'manifest.csv'
    manifest.write_text('id,state\nm,neu\n')
    assert common.read_manifest() == [{'id': 'm', 'state': 'neu'}]
    assert common.read_manifest() == [{'id': 'm', 'state': 'neu'}]  # cached, still correct
    # A rewrite (different size/mtime) must invalidate the cache, not serve stale rows.
    manifest.write_text('id,state\nm,neu\nm,red1\n')
    assert common.read_manifest() == [{'id': 'm', 'state': 'neu'}, {'id': 'm', 'state': 'red1'}]
    manifest.unlink()
    assert common.read_manifest() == []


def test_config_check_failure_cannot_claim_runtime_ready(monkeypatch):
    from redox.validation import uma_environment
    monkeypatch.setattr(uma_environment.subprocess, 'run', lambda *a, **kw: (_ for _ in ()).throw(FileNotFoundError('python')))
    result = uma_environment.inspect_policy(dict(nitro_python='/missing',legacy_python='/missing',
        fairchem_cache='/cache',nitro_library_path='/lib'))
    assert not result['configuration_ok'] and len(result['errors']) == 2
    assert not result['production_qualified'] and not result['cuda_context_created']
