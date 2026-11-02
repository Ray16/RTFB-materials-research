"""Discover the authoritative reservation tool without a checkout-specific path."""
import os
from pathlib import Path
import shutil


def reservation_gate(*, required=True):
    explicit = os.environ.get('REDOX_GPU_RESERVE')
    candidates = [explicit] if explicit else [shutil.which('gpu_reserve'), Path.home()/'bin/gpu_reserve']
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
    return Path(os.environ.get('GPU_LOCK_ROOT', str(Path.home()/'.gpu_locks'))).expanduser().resolve()
