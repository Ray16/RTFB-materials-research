#!/usr/bin/env python
"""Bound-aware electronic dispatch with a fixed remote working directory.

Numerical producers and plans remain untouched. Uncertain SSH placement errors
are persisted while the original fail-closed ownership behavior is retained.
"""
import inspect
import os
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


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
    if body.count(old) != 1:
        raise ValueError('Electronic controller error-handling template changed')
    scope = {}
    exec(compile(body.replace(old, new, 1), str(Path(__file__).resolve()), 'exec'), vars(module), scope)
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
