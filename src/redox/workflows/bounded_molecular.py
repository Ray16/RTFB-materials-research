"""Shared execution/budget adapter for source-frozen molecular namespaces.

Callers must first enforce their scientific pilot and workload admission rules.
The existing atomic reservation, node probe, CPU/RAM and process owner checks
remain in the underlying molecular launcher and budget observer.
"""
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def launch(namespace, gpu, output, *, caller_source):
    from redox.workflows.molecular_native_namespace import verify, launch as execute
    from redox.workflows.screening_control import probe
    from redox.workflows.response_diagnostic import observe_budget
    import json
    spec = verify(namespace)
    plan = json.loads((Path(spec['case'])/'plan.json').read_text())
    probe(namespace)
    child, result = execute(namespace, gpu, output)
    source_paths = [(Path(observe_budget.__code__.co_filename), 'budget_observer_source.py'),
                    (Path(__file__), 'bounded_molecular_source.py'),
                    (Path(caller_source), 'scientific_launcher_source.py')]
    for source, name in source_paths:
        shutil.copyfile(source, output/name)
    result.update(maximum_wall_seconds=plan['maximum_wall_seconds'],
        maximum_durable_output_bytes=plan['maximum_durable_output_bytes'],
        budget_observer_owner=dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19]),
        budget_observer_source_sha256=sha256_file(output/'budget_observer_source.py'),
        bounded_launcher_sources={str(output/name): sha256_file(output/name) for _, name in source_paths})
    atomic_write_json(output/'launch.json', result)
    return observe_budget(child, spec, result, output)
