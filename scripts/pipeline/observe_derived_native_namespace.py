#!/usr/bin/env python3
"""Use the existing physical native observer with derived-batch mountpoints.

The adapter changes only filesystem construction. Scientific review, closed
ownership, raw-source verification and missing-result handling remain in the
maintained observer. Its code is not edited while earlier collectors are live.
"""
import importlib.util
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.derived_native_namespace import namespace_command


def main():
    original = ROOT/'scripts/pipeline/observe_archived_native_namespace.py'
    sources = {str(path): sha256_file(path) for path in
               (Path(__file__).resolve(), original, ROOT/'src/redox/workflows/derived_native_namespace.py')}
    index = sys.argv.index('--out')+1
    output = Path(sys.argv[index]).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(output.with_name(output.name+'_adapter.json'),
                      dict(status='filesystem_adapter_only_scientific_gates_unchanged', source_sha256=sources))
    specification = importlib.util.spec_from_file_location('unchanged_native_observer_for_derived_batches', original)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    def build(spec):
        if any(sha256_file(path) != digest for path, digest in sources.items()):
            raise ValueError('Derived namespace/observer source changed before review')
        command = namespace_command(spec)
        # The identical shared binary is available even on hosts without bwrap.
        portable = ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
        if portable.is_file():
            command[0] = str(portable)
        return command

    module.mount_command = build
    module.main()


if __name__ == '__main__':
    main()
