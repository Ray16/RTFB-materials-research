#!/usr/bin/env python
"""Bound-aware electronic dispatch with safe gate-contention recovery.

Numerical producers and plans remain untouched. Uncertain SSH placement errors
are persisted while the original fail-closed ownership behavior is retained.
"""
import inspect
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def recover_gate_contention(case, handle, owner_confirmed_dead):
    """Archive and reopen only an explicit atomic gate refusal with no worker output."""
    if not owner_confirmed_dead or handle.get('status') != 'gate_started_scientific_start_pending':
        return None
    case=Path(case).resolve();launch=case/'launch.json';intent=case/'dispatch_intent.json'
    if not launch.is_file() or not intent.is_file() or json.loads(launch.read_text()) != handle:
        return None
    attempt=Path(handle.get('attempt','')).resolve()
    if not attempt.is_relative_to(case/'dispatch_attempts') or attempt.parent != case/'dispatch_attempts':
        return None
    log=attempt/'worker.log'
    if not log.is_file() or not re.fullmatch(r'gpu_reserve: gpu'+str(handle.get('GPU'))+r' BUSY/reserved -> not running\s*',log.read_text()):
        return None
    plan=json.loads((case/'plan.json').read_text())
    if (Path(plan['raw_root']).exists() or any((case/name).exists() for name in
            ('experimental_energy','residency_guard','host_guard','tensor_runtime'))):
        return None
    launch_bytes,intent_bytes=launch.read_bytes(),intent.read_bytes()
    history=case/'gate_contention_history'/str(time.time_ns());history.mkdir(parents=True,exist_ok=False)
    (history/'launch.json').write_bytes(launch_bytes)
    (history/'intent.json').write_bytes(intent_bytes)
    shutil.copyfile(log,history/'worker.log')
    audit=dict(status='atomic_gate_refusal_no_scientific_attempt',checked=time.time(),
        owner=handle['owner'],GPU=handle['GPU'],attempt=str(attempt),
        launch_sha256=api.digest(history/'launch.json'),intent_sha256=api.digest(history/'intent.json'),
        worker_log_sha256=api.digest(history/'worker.log'),scientific_attempts=0)
    atomic_write_json(history/'audit.json',audit)
    if launch.read_bytes()!=launch_bytes or intent.read_bytes()!=intent_bytes:
        raise ValueError('Gate contention records changed during recovery')
    if Path(plan['raw_root']).exists() or any((case/name).exists() for name in
            ('experimental_energy','residency_guard','host_guard','tensor_runtime')):
        raise ValueError('Scientific output appeared during gate contention recovery')
    launch.unlink()
    prior=json.loads(intent_bytes);prior.update(outcome='not_started',contention_audit=str(history/'audit.json'))
    atomic_write_json(intent,prior)
    return str(history/'audit.json')


def coordinate_with_durable_errors(module):
    body = inspect.getsource(module.coordinate)
    old = """                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            observations[entry['key']]=dict(owner=dict(host=host,pid=-1,process_start='unknown'),error=repr(error))
                            pending.pop(0);continue
"""
    new = """                        except (OSError,ValueError,subprocess.SubprocessError) as error:
                            detail=dict(owner=dict(host=host,pid=-1,process_start='unknown'),error=repr(error),
                                stderr=str(getattr(error,'stderr','') or '')[-4000:],checked=time.time(),
                                status='placement_response_unverified_preserved')
                            observations[entry['key']]=detail
                            atomic_write_json(case/'placement_error.json',detail)
                            pending.pop(0);continue
"""
    stopped = """                        else:
                            item['status']='stopped_attempt_preserved_review_required'
                            if key not in processed:
                                host=handle['owner']['host'];failures[host]=failures.get(host,0)+1;processed.add(key)
"""
    recovered = """                        else:
                            gate_audit=_recover_gate_contention(case,handle,owner_confirmed_dead=True)
                            if gate_audit:
                                observations.pop(key,None)
                                item.update(status='pending',last_gate_contention=gate_audit)
                                pending.append(entry)
                            else:
                                item['status']='stopped_attempt_preserved_review_required'
                                if key not in processed:
                                    host=handle['owner']['host'];failures[host]=failures.get(host,0)+1;processed.add(key)
"""
    if body.count(old) != 1 or body.count(stopped) != 1:
        raise ValueError('Electronic controller error-handling template changed')
    scope = {}
    namespace=dict(vars(module),_recover_gate_contention=recover_gate_contention)
    exec(compile(body.replace(old,new,1).replace(stopped,recovered,1),str(Path(__file__).resolve()),'exec'),namespace,scope)
    return scope['coordinate']


def main():
    # SSH starts in the remote account's home. Some frozen evidence paths are
    # repository-relative; resolve them identically on every node.
    os.chdir(ROOT)
    upstream = ROOT / 'scripts/pipeline/coordinate_smooth_electronic_batch.py'
    adapter_path = ROOT / 'scripts/pipeline/dispatch_smooth_energy_capacity.py'
    proof = {str(p): api.digest(p) for p in (Path(__file__).resolve(), upstream, adapter_path)}
    original_load, original_verify = api.load_module, api.verify_sources
    adapter = original_load('bound_aware_energy_capacity_adapter_v2', adapter_path)

    def verify(mapping):
        original_verify(proof)
        original_verify(mapping)

    def load(name, path):
        verify({})
        loaded = original_load(name, path)
        if Path(path).resolve() == ROOT / 'scripts/pipeline/dispatch_smooth_energy_pilot.py':
            loaded.placement_module = adapter.capacity_placement
            loaded.preflight = adapter.capacity_preflight(loaded)
        return loaded

    api.load_module, api.verify_sources = load, verify
    module = load('electronic_controller_capacity_workdir_v2', upstream)
    module.SCRIPT = Path(__file__).resolve()
    module.coordinate = coordinate_with_durable_errors(module)
    sys.argv = [str(Path(__file__).resolve()), *sys.argv[1:]]
    module.main()


if __name__ == '__main__':
    main()
