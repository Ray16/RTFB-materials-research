#!/usr/bin/env python3
"""CPU-safe fleet preflight; GPU imports occur only inside reserved workers.

Retains the canonical skill's CPU/RAM/path/SSH checks. GPU packages are checked
through distribution metadata, with a process-wide import guard against direct
and transitive GPU-library imports. The shared skill installation is untouched.
"""
import importlib
from importlib import metadata
import importlib.util
import inspect
from pathlib import Path

GPU_DISTRIBUTIONS = {'gpu4pyscf':'gpu4pyscf-cuda12x', 'cupy':'cupy-cuda12x',
                     'torch':'torch', 'fairchem':'fairchem-core'}


def module_version(name, importer=importlib.import_module, version=metadata.version):
    root = name.split('.')[0]
    if root in GPU_DISTRIBUTIONS:
        return version(GPU_DISTRIBUTIONS[root]), True
    value = importer(name)
    return str(getattr(value, '__version__', 'import OK')), False


def safe_remote(original):
    anchor = 'import importlib, json, os, pathlib, sys'
    if original.count(anchor) != 1:
        raise ValueError('Canonical preflight import template changed')
    setup = '''
from importlib import metadata
from importlib.abc import MetaPathFinder
GPU_DISTRIBUTIONS = '''+repr(GPU_DISTRIBUTIONS)+'''
class _CPUOnlyImports(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in GPU_DISTRIBUTIONS:
            raise ImportError('GPU imports are deferred until inside a gpu_reserve worker: '+fullname)
sys.meta_path.insert(0,_CPUOnlyImports())
'''+inspect.getsource(module_version)+'''
'''
    original = original.replace(anchor, anchor+'\n'+setup, 1)
    anchor = "report={'host':os.uname().nodename,'python':sys.executable,'ok':False,'errors':[]}"
    if original.count(anchor) != 1:
        raise ValueError('Canonical preflight report template changed')
    original = original.replace(anchor, anchor+"\nreport.update(pid=os.getpid(),GPU_imports_deferred=[],GPU_imports_blocked_in_CPU_preflight=True)", 1)
    anchor = "            module=importlib.import_module(name)\n            versions[name]=str(getattr(module,'__version__','import OK'))"
    replacement = """            value,deferred=module_version(name)
            versions[name]=value
            if deferred:report['GPU_imports_deferred'].append(name)"""
    if original.count(anchor) != 1:
        raise ValueError('Canonical preflight dependency template changed')
    original = original.replace(anchor, replacement, 1)
    anchor = "    report['modules']=versions"
    if original.count(anchor) != 1:
        raise ValueError('Canonical preflight result template changed')
    return original.replace(anchor, anchor+"\n    report['GPU_modules_loaded']=[n for n in sys.modules if n.split('.')[0] in GPU_DISTRIBUTIONS]\n    if report['GPU_modules_loaded']:report['errors'].append('GPU module loaded in CPU-only preflight')", 1)


def main():
    root = Path(__file__).resolve().parents[2]
    source = root.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'
    if not source.is_file():
        raise RuntimeError('The required canonical lambda-fleet preflight could not be found')
    spec = importlib.util.spec_from_file_location('canonical_CPU_dependency_preflight',source)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.REMOTE = safe_remote(module.REMOTE)
    module.main()


if __name__ == '__main__': main()
