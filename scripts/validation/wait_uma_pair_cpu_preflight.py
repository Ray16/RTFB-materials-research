#!/usr/bin/env python
"""Resume the fixed-device UMA comparison using CPU-only dependency preflight.

The paired scientific producer is unchanged. The original waiter is executed
with a process-local substitution of the qualified CPU-only preflight helper.
"""
from pathlib import Path
import runpy
import subprocess

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


UPSTREAM = ROOT / 'scripts/validation/wait_uma_pair.py'
ORIGINAL = ROOT.parent / '.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
CPU_ONLY = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'


def preflight_command(argv):
    values = list(argv)
    indices = [i for i, value in enumerate(values) if str(value) == str(ORIGINAL)]
    if len(indices) > 1:
        raise ValueError('Ambiguous dependency preflight command')
    if indices:
        values[indices[0]] = str(CPU_ONLY)
    return values


def main():
    paths = (Path(__file__).resolve(), UPSTREAM, ORIGINAL, CPU_ONLY)
    proof = {str(path): api.digest(path) for path in paths}
    original_verify, original_run = api.verify_sources, subprocess.run

    def verify(mapping):
        original_verify(proof)
        return original_verify(mapping)

    def run(argv, **kwargs):
        verify({})
        return original_run(preflight_command(argv), **kwargs)

    api.verify_sources, subprocess.run = verify, run
    try:
        runpy.run_path(str(UPSTREAM), run_name='__main__')
    finally:
        api.verify_sources, subprocess.run = original_verify, original_run


if __name__ == '__main__':
    main()
