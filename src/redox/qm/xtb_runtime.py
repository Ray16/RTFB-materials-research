"""Explicit, hash-qualified xTB path binding for thermochemistry entry points.

This does not alter producer source files, certificates or native physical gates.
It permits only the previously qualified library and executable bytes. A caller
must bind once before starting a qRRHO worker, within a single Python thread.
"""
from contextlib import contextmanager
import hashlib
import os
from pathlib import Path
import threading

QUALIFIED_XTB_BINARY_SHA256='271dd6cdd5da127a8e39919e7e32d924790fb346f0057b4183769ba8b3411887'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def validate_prefix(prefix,library_sha256,binary_sha256=QUALIFIED_XTB_BINARY_SHA256):
    prefix=Path(prefix).expanduser().resolve()
    library=(prefix/'lib/libxtb.so.6.7.1').resolve();binary=(prefix/'bin/xtb').resolve()
    if (not library.is_relative_to(prefix) or not binary.is_relative_to(prefix)
            or not library.is_file() or not binary.is_file() or not os.access(binary,os.X_OK)):
        raise ValueError('xTB library and executable must exist inside the declared prefix')
    if digest(library)!=library_sha256 or digest(binary)!=binary_sha256:
        raise ValueError('xTB files differ from the qualified bytes')
    return dict(prefix=str(prefix),library=str(library),library_sha256=library_sha256,
                binary=str(binary),binary_sha256=binary_sha256)


@contextmanager
def qualified_xtb_prefix(prefix):
    """Bind the exact qualified files process-locally and restore them afterwards."""
    if threading.current_thread() is not threading.main_thread() or threading.active_count()!=1:
        raise RuntimeError('Bind xTB before starting Python worker threads')
    from redox.qm import qrrho, native_thermal, projected_thermal
    if qrrho._WORKER is not None or qrrho._LAST_WORKER is not None:
        raise RuntimeError('Bind xTB before the first thermodynamic worker request')
    runtime=validate_prefix(prefix,qrrho.LIBRARY_SHA256)
    previous=(qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB,os.environ.get('LD_LIBRARY_PATH'))
    qrrho.LIBRARY=Path(runtime['library']);qrrho._checked_library.cache_clear()
    native_thermal.XTB=projected_thermal.XTB=Path(runtime['binary'])
    os.environ['LD_LIBRARY_PATH']=str(Path(runtime['library']).parent)+(os.pathsep+previous[3] if previous[3] else '')
    try:
        qrrho._checked_library()
        yield runtime
    finally:
        try:
            qrrho._close_worker()
        finally:
            qrrho._LAST_WORKER=None
            qrrho.LIBRARY,native_thermal.XTB,projected_thermal.XTB=previous[:3]
            qrrho._checked_library.cache_clear()
            if previous[3] is None:os.environ.pop('LD_LIBRARY_PATH',None)
            else:os.environ['LD_LIBRARY_PATH']=previous[3]
