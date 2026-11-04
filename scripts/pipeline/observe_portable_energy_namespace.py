#!/usr/bin/env python3
"""Use the unchanged energy observer with a bound portable namespace executable."""
import importlib.util
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.archived_energy_namespace import energy_mount


def main():
    original = ROOT/'scripts/pipeline/observe_archived_energy_namespace.py'
    portable = ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    if not portable.is_file():
        raise FileNotFoundError('Qualified shared namespace executable required')
    sources = {str(path): sha256_file(path) for path in
               (Path(__file__).resolve(), original, portable, ROOT/'src/redox/workflows/archived_energy_namespace.py')}
    output = Path(sys.argv[sys.argv.index('--out')+1]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output.with_name(output.name+'_adapter.json'),
                      dict(status='namespace_executable_adapter_only', source_sha256=sources,
                           scientific_review_or_thresholds_changed=False))
    specification = importlib.util.spec_from_file_location('unchanged_energy_observer_with_portable_binary', original)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    def build(spec):
        if any(sha256_file(path) != digest for path, digest in sources.items()):
            raise ValueError('Energy namespace executable/adapter changed before review')
        command = energy_mount(spec)
        command[0] = str(portable)
        return command

    module.energy_mount = build
    module.main()


if __name__ == '__main__':
    main()
