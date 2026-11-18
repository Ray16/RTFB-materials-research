"""Shared byte-preserving readers for native orbital checkpoint arrays.

This module validates storage layout only; callers own state/frame and physical QC.
"""
import hashlib
import math


def array_digest(array):
    """SHA-256 over an array's raw bytes, without copying it into a new buffer."""
    return hashlib.sha256(memoryview(array).cast('B')).hexdigest()


def to_host(value):
    """Return a host NumPy view of a possibly device-resident array, unchanged.

    Device backends expose a ``.get()`` host transfer; host arrays pass through.
    Values are copied verbatim, so checkpoint bytes and digests stay exact.
    """
    import numpy as np
    return np.asarray(value.get() if hasattr(value, 'get') else value)


def contiguous_payload(handle,name,content):
    """Read only native FP64 contiguous payloads with verified in-file bounds.

    The installed HDF5 high-level read stalls on the real source checkpoint.
    Metadata offsets avoid that path; orbital SHA-256 checks below remain exact.
    Chunked, compressed, external or unallocated datasets are unsupported.
    """
    import h5py
    import numpy as np
    dataset=handle[name];creation=dataset.id.get_create_plist();offset=dataset.id.get_offset()
    if (dataset.dtype!=np.dtype('float64') or creation.get_layout()!=h5py.h5d.CONTIGUOUS
            or creation.get_nfilters()!=0 or creation.get_external_count()!=0
            or offset is None or offset<0):raise ValueError('Only allocated native FP64 contiguous checkpoint datasets are supported')
    count=math.prod(dataset.shape);size=count*8
    if offset+size>len(content):raise ValueError('Contiguous dataset extends beyond the source checkpoint')
    return np.frombuffer(content,dtype=np.float64,count=count,offset=offset).reshape(dataset.shape)
