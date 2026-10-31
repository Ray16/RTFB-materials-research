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
import sys

QUALIFIED_XTB_BINARY_SHA256='271dd6cdd5da127a8e39919e7e32d924790fb346f0057b4183769ba8b3411887'
LEGACY_PREFIX=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox')


def _absolute(value,name):
    if not value or not Path(value).expanduser().is_absolute():
        raise ValueError(name+' must be a nonempty absolute path')
    return Path(value).expanduser().resolve()


def configured_prefix(environ=None):
    """Find a complete xTB installation without importing numerical libraries.

    An explicit prefix wins. Normal conda/venv interpreter locations follow;
    the historical cluster location is only a compatibility fallback. Hash
    qualification remains the responsibility of the thermodynamic entry point.
    """
    environ=os.environ if environ is None else environ
    if 'REDOX_XTB_PREFIX' in environ:return _absolute(environ['REDOX_XTB_PREFIX'],'REDOX_XTB_PREFIX')
    for prefix in dict.fromkeys((Path(sys.prefix),Path(sys.base_prefix),LEGACY_PREFIX)):
        if (prefix/'lib/libxtb.so.6.7.1').is_file() and (prefix/'bin/xtb').is_file():return prefix.resolve()
    return Path(sys.prefix).resolve()


def configured_library(environ=None):
    environ=os.environ if environ is None else environ
    if 'REDOX_XTB_LIBRARY' in environ:return _absolute(environ['REDOX_XTB_LIBRARY'],'REDOX_XTB_LIBRARY')
    return configured_prefix(environ)/'lib/libxtb.so.6.7.1'


def configured_binary(environ=None):
    environ=os.environ if environ is None else environ
    if 'REDOX_XTB_EXECUTABLE' in environ:return _absolute(environ['REDOX_XTB_EXECUTABLE'],'REDOX_XTB_EXECUTABLE')
    return configured_prefix(environ)/'bin/xtb'


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
