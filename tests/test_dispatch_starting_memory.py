"""Exercise asynchronous launch promises and exact stopped prelaunch retries."""
import json
from pathlib import Path
from types import SimpleNamespace

from redox.workflows import charged_native_followup as flow


def setup(tmp_path, monkeypatch, bounds):
    output = tmp_path / 'run'; output.mkdir(); rows = []; specs = {}; launched = []
    for i, bound in enumerate(bounds):
        root = output / 'states' / str(i); case = root / 'case'; case.mkdir(parents=True)
        (case / 'plan.json').write_text(json.dumps(dict(maximum_host_RSS_GiB=bound,
            maximum_durable_output_bytes=256 * 2**20, minimum_GPU_memory_MiB=31744)))
        ns = root / 'namespace.json'; ns.write_text(json.dumps(dict(case=str(case)))); specs[str(ns)] = dict(case=str(case))
        rows.append(dict(key=str(i), namespace=str(ns)))
    policy = tmp_path / 'policy.env'; policy.write_text('ALLOWED_HOSTS="lambdaA lambdaB"\n')
    monkeypatch.setenv('GPU_RESERVE_CONFIG', str(policy))
    monkeypatch.setattr(flow, 'verify', lambda path: specs[str(path)])
    # These tests isolate RAM promises and claim ownership. Real campaign
    # priority declarations must not trigger SSH against the synthetic hosts.
    monkeypatch.setattr(flow, 'priority_requests', lambda path: [])
    capacity = SimpleNamespace(discover_slots=lambda hosts: ([('lambdaA', i) for i in range(3)] + [('lambdaB', 0)], []),
        capacity=lambda *args: dict(GPU_MiB=32768, available=503 * 2**30, promised=0, load=0, cpus=80),
        fits=lambda plan, data: data['available'] - data['promised'] >= (plan['maximum_host_RSS_GiB'] + 32) * 2**30)
    def inspect(row, claim, alive):
        return ('prelaunch_refusal_no_scientific_attempt', {}) if claim.exists() else ('pending', {})
    def launch(row, host, gpu, interpreter, claim, number):
        launched.append((row['key'], host, number))
        receipt = dict(placement=str(Path(row['namespace']).parent / f'placement_{number:03d}'))
        (claim / f'try_{number:03d}.json').write_text(json.dumps(receipt))
        return receipt
    modules = dict(charged_native_capacity=capacity, charged_native_owners=SimpleNamespace(alive=lambda owner: False),
        charged_native_shared_launcher=SimpleNamespace(inspect=inspect, launch=launch),
        charged_native_shared_storage=SimpleNamespace(storage_headroom=lambda *args: dict(sufficient=True)))
    # Dispatch freezes every loaded module's source. Test doubles must expose
    # that same contract so the tests reach the RAM/ownership decisions.
    for name, module in modules.items():
        source = tmp_path / f'{name}.py'
        source.write_text(f'# Isolated dispatcher test double: {name}\n')
        module.__file__ = str(source)
    monkeypatch.setattr(flow, 'load', lambda path, name: modules[name])
    return output, rows, launched, modules


def test_not_yet_started_workers_cannot_consume_the_same_host_RAM(tmp_path, monkeypatch):
    output, rows, launched, _ = setup(tmp_path, monkeypatch, [240, 160, 160])
    result = flow.dispatch_rows(rows, output, 'test-python', [])
    assert launched == [('0', 'lambdaA', 1), ('1', 'lambdaA', 1), ('2', 'lambdaB', 1)]
    assert len(result['launches']) == 3 and not result['refusals']
    assert all(Path(p).is_file() for p in result['source_sha256'])


def test_large_queue_head_does_not_hide_a_fitting_ready_input(tmp_path, monkeypatch):
    output, rows, launched, _ = setup(tmp_path, monkeypatch, [240, 240, 160])
    flow.dispatch_rows(rows, output, 'test-python', [])
    assert launched == [('0', 'lambdaA', 1), ('2', 'lambdaA', 1), ('1', 'lambdaB', 1)]


def test_stopped_prelaunch_retry_preserves_first_receipt_and_number(tmp_path, monkeypatch):
    output, rows, launched, _ = setup(tmp_path, monkeypatch, [160])
    claim = Path(rows[0]['namespace']).parent / 'dispatch_claim'; claim.mkdir()
    (claim / 'owner.json').write_text(json.dumps(dict(queue=str(output))))
    first = claim / 'try_001.json'; first.write_text('{"old":true}')
    flow.dispatch_rows(rows, output, 'test-python', [])
    assert launched == [('0', 'lambdaA', 2)] and first.read_text() == '{"old":true}'


def test_unverified_or_foreign_claims_never_restart(tmp_path, monkeypatch):
    output, rows, launched, modules = setup(tmp_path, monkeypatch, [160])
    claim = Path(rows[0]['namespace']).parent / 'dispatch_claim'; claim.mkdir()
    (claim / 'owner.json').write_text(json.dumps(dict(queue='different-run')))
    result = flow.dispatch_rows(rows, output, 'test-python', [])
    assert not launched and 'foreign dispatch claim' in result['refusals'][0]['error']
    modules['charged_native_shared_launcher'].inspect = lambda *args: ('ownership_unverified', {})
    flow.dispatch_rows(rows, output, 'test-python', [])
    assert not launched
