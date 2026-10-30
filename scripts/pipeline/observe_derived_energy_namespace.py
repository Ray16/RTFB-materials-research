#!/usr/bin/env python3
"""Collect measured-class namespace energies with the original physical audit."""
import importlib.util
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.derived_energy_namespace import mount


def main():
    source = ROOT/'scripts/pipeline/observe_archived_energy_namespace.py'
    paths = (Path(__file__).resolve(), source,
             ROOT/'src/redox/workflows/derived_energy_namespace.py',
             ROOT/'src/redox/workflows/derived_native_namespace.py',
             ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap')
    hashes = {str(path): sha256_file(path) for path in paths}
    output = Path(sys.argv[sys.argv.index('--out')+1]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output.with_name(output.name+'_adapter.json'),
                      dict(status='filesystem_adapter_only_independent_energy_gates_unchanged',
                           source_sha256=hashes, native_G_computed=False))
    specification = importlib.util.spec_from_file_location('measured_class_original_energy_observer', source)
    observer = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(observer)

    def command(spec):
        if any(sha256_file(path) != digest for path, digest in hashes.items()):
            raise ValueError('Namespace adapter changed before independent energy review')
        return mount(spec)

    observer.energy_mount = command
    observer.main()


if __name__ == '__main__':
    main()
