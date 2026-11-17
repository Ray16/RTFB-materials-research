"""Discover the authoritative reservation tool without a checkout-specific path."""
import os
import json
from pathlib import Path
import shutil
import importlib.util


def site_defaults():
    path=Path(__file__).resolve().parents[3]/'config/runtime_paths.json'
    if not path.exists():return {}
    data=json.loads(path.read_text())
    if data.get('schema')!=1 or not isinstance(data.get('site_defaults'),dict):
        raise ValueError('Unsupported runtime-path configuration')
    return data['site_defaults']


def reservation_gate(*, required=True):
    explicit = os.environ.get('REDOX_GPU_RESERVE')
    candidates = [explicit] if explicit else [site_defaults().get('gpu_reserve'), shutil.which('gpu_reserve'), Path.home()/'bin/gpu_reserve']
    for candidate in candidates:
        if candidate:
            path = Path(candidate).expanduser().resolve()
            if path.is_file() and os.access(path, os.X_OK):
                return path
    if not required and not explicit:
        # CPU-only imports and documentation must work without a GPU service.
        # Actual execution of this absent path fails; it never bypasses the gate.
        return (Path.home()/'bin/gpu_reserve').resolve()
    raise RuntimeError('Authoritative gpu_reserve not found; set REDOX_GPU_RESERVE to the installed tool')


def reservation_lock_root():
    # This environment name is also consumed by gpu_reserve itself. Using it
    # keeps ownership checks and the gate in the same reservation namespace.
    if 'GPU_LOCK_ROOT' in os.environ:
        return Path(os.environ['GPU_LOCK_ROOT']).expanduser().resolve()
    configured=site_defaults();gate=configured.get('gpu_reserve');locks=configured.get('gpu_lock_root')
    if gate and locks and reservation_gate(required=False)==Path(gate).expanduser().resolve():
        return Path(locks).expanduser().resolve()
    return (Path.home()/'.gpu_locks').resolve()


def gpu4pyscf_source_root():
    """Prefer declared shared sources only when their actual bytes match locally.

    Preparation may run in a node-local environment while the worker executes
    from a shared pinned environment. Identical paths are not assumed; the
    three implementation bindings must be byte-identical before relocation.
    """
    from redox.core.hashing import sha256_file
    spec=importlib.util.find_spec('gpu4pyscf')
    if spec is None or not spec.origin:raise RuntimeError('Installed GPU4PySCF package source not found')
    installed=Path(spec.origin).resolve().parent
    declared=os.environ.get('REDOX_GPU4PYSCF_ROOT') or site_defaults().get('gpu4pyscf_source_root')
    if not declared:return installed
    shared=Path(declared).expanduser().resolve()
    for name in ('dft/gen_grid.py','dft/numint.py','lib/libgdft.so'):
        if not (shared/name).is_file() or sha256_file(shared/name)!=sha256_file(installed/name):
            raise ValueError('Declared shared GPU4PySCF implementation differs: '+name)
    return shared
